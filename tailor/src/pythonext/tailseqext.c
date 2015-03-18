/**
 * Assorted accelerator functions for TAIL-seq from Narry Kim Lab.
 *
 * Copyright (c) 2014-5 Institute of Basic Sciences
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Hyeshik Chang <hyeshik@snu.ac.kr>
 */

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <numpy/arrayobject.h>

#define NCHANNELS                   4
#define INTENSITY_CODING_BASE       33
#define INTENSITY_CODING_RADIX      91
#define INTENSITY_CODING_WIDTH      8192
#define INTENSITY_BOTTOM_SHIFT      255

static int
decode_intensity(const char *encoded, int length, int16_t *signals)
{
    int chan, cycle;

    for (cycle = 0; cycle < length; cycle++) {
        for (chan = 0; chan < NCHANNELS; chan++) {
            int16_t high, low;

            high = (int16_t)*(encoded++) - INTENSITY_CODING_BASE;
            low = (int16_t)*(encoded++) - INTENSITY_CODING_BASE;

            if (high < 0 || high >= INTENSITY_CODING_RADIX ||
                    low < 0 || low >= INTENSITY_CODING_RADIX) {
                PyErr_SetString(PyExc_ValueError, "Illegal character found.");
                return -1;
            }

            *(signals++) = (high * INTENSITY_CODING_RADIX + low) - INTENSITY_BOTTOM_SHIFT;
        }
    }

    return 0;
}

static PyObject *
Py_decode_intensity(PyObject *self, PyObject *args)
{
    PyArrayObject *r;
    char *encoded;
    Py_ssize_t encoded_length;
    npy_intp dim[2];

    if (!PyArg_ParseTuple(args, "s#:decode_intensity", &encoded, &encoded_length))
        return NULL;

    if (encoded_length % 8 != 0) {
        PyErr_Format(PyExc_ValueError, "Irregular length of encoded intensity: \"%s\"",
                     encoded);
        return NULL;
    }

    dim[0] = encoded_length / 8;
    dim[1] = NCHANNELS;

    if ((r = (PyArrayObject *)PyArray_SimpleNew(2, dim, NPY_INT16)) == NULL)
        return NULL;

    if (decode_intensity(encoded, dim[0], PyArray_DATA(r)) == -1)
        Py_XDECREF(r);

    return (PyObject *)r;
}


/* ======================================================== */
/* Poly(A) locator */

typedef struct {
    PyObject_HEAD
    int weight[256];
} PolyALocatorObject;

static PyTypeObject PolyALocator_Type;

#define PolyALocator_Type(v)        (Py_TYPE(v) == &PolyALocator_Type)

static PolyALocatorObject *
PolyALocator_create(PyObject *args)
{
    PolyALocatorObject *self;
    PyObject *weightmap=NULL;
    char c;

    if (!PyArg_ParseTuple(args, "O:PolyALocator", &weightmap))
        return NULL;

    if (!PyDict_Check(weightmap)) {
        PyErr_SetString(PyExc_TypeError, "Weights must be given in a dict.");
        return NULL;
    }

    self = PyObject_New(PolyALocatorObject, &PolyALocator_Type);
    if (self == NULL)
        return NULL;

    memset(self->weight, 0, sizeof(int) * 256);
    for (c = 'A'; c <= 'Z'; c++) {
        char key[2];
        PyObject *value;

        key[0] = c; key[1] = 0;

        value = PyDict_GetItemString(weightmap, key);
        if (value == NULL)
            continue;

        if (!PyLong_Check(value)) {
            PyErr_SetString(PyExc_TypeError, "Values of weight map must be an integer");
            Py_DECREF(self);
            return NULL;
        }

        self->weight[(int)c] = (int)PyLong_AsLong(value);
    }

    return self;
}

static void
PolyALocator_dealloc(PolyALocatorObject *self)
{
    PyObject_Del(self);
}

static PyObject *
PolyALocator_call(PolyALocatorObject *self, PyObject *args, PyObject *kw)
{
    char *inseq;
    Py_ssize_t inseqlen;
    int max_term_mod, i, j;
    int longest_i, longest_j, longest_length;

    if (!PyArg_ParseTuple(args, "s#i:decode", &inseq, &inseqlen, &max_term_mod))
        return NULL;

    if (inseqlen < max_term_mod)
        max_term_mod = inseqlen;

    longest_i = longest_j = inseqlen + 1;
    longest_length = -1;

    /* calculate match scores for all possible [i, j] */
    for (i = 0; i < max_term_mod; i++) {
        int curlength, scoresum;

        scoresum = self->weight[(int)inseq[i]];
        if (longest_length < 1 && scoresum > 0) {
            longest_length = 1;
            longest_i = longest_j = i;
        }

        for (j = i + 1, curlength = 2; j < inseqlen; j++, curlength++) {
            scoresum += self->weight[(int)inseq[j]];

            if (scoresum > 0 && curlength > longest_length) {
                longest_i = i;
                longest_j = j;
                longest_length = curlength;
            }
        }
    }

    if (longest_length < 0)
        return Py_BuildValue("(ii)", 0, 0);

    while (inseq[longest_i] != 'T' && longest_i <= longest_j)
        longest_i++;

    while (inseq[longest_j] != 'T' && longest_j >= longest_i)
        longest_j--;

    return Py_BuildValue("(ii)", longest_i, longest_j + 1);
}


static PyTypeObject PolyALocator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tailseqext.PolyALocator",           /*tp_name*/
    sizeof(PolyALocatorObject),          /*tp_basicsize*/
    0,                      /*tp_itemsize*/
    /* methods */
    (destructor)PolyALocator_dealloc, /*tp_dealloc*/
    0,                      /*tp_print*/
    0,                      /*tp_getattr*/
    0,                      /*tp_setattr*/
    0,                      /*tp_compare*/
    0,                      /*tp_repr*/
    0,                      /*tp_as_number*/
    0,                      /*tp_as_sequence*/
    0,                      /*tp_as_mapping*/
    0,                      /*tp_hash*/
    (ternaryfunc)PolyALocator_call, /*tp_call*/
    0,                      /*tp_str*/
    0,                      /*tp_getattro*/
    0,                      /*tp_setattro*/
    0,                      /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,     /*tp_flags*/
    0,                      /*tp_doc*/
    0,                      /*tp_traverse*/
    0,                      /*tp_clear*/
    0,                      /*tp_richcompare*/
    0,                      /*tp_weaklistoffset*/
    0,                      /*tp_iter*/
    0,                      /*tp_iternext*/
    0,                      /*tp_methods*/
    0,                      /*tp_members*/
    0,                      /*tp_getset*/
    0,                      /*tp_base*/
    0,                      /*tp_dict*/
    0,                      /*tp_descr_get*/
    0,                      /*tp_descr_set*/
    0,                      /*tp_dictoffset*/
    0,                      /*tp_init*/
    0,                      /*tp_alloc*/
    0,                      /*tp_new*/
    0,                      /*tp_free*/
    0,                      /*tp_is_gc*/
};

static PyObject *
PolyALocator_new(PyObject *self, PyObject *args)
{
    return (PyObject *)PolyALocator_create(args);
}


/* List of functions defined in the module */

static PyMethodDef tailseqext_methods[] = {
    {"decode_intensity", Py_decode_intensity, METH_VARARGS,
        PyDoc_STR("decode_intensity(str) -> array")},
    {"PolyALocator", PolyALocator_new, METH_VARARGS,
        PyDoc_STR("PolyALocator(weights) -> object")},
    {NULL, NULL}           /* sentinel */
};

PyDoc_STRVAR(tailseqext_doc,
"This is a template module just for instruction.");

struct tailseqext_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct tailseqext_state*)PyModule_GetState(m))

static int tailseqext_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int tailseqext_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef tailseqext_def = {
    PyModuleDef_HEAD_INIT,
    "tailseqext",
    tailseqext_doc,
    sizeof(struct tailseqext_state),
    tailseqext_methods,
    NULL,
    tailseqext_traverse,
    tailseqext_clear,
    NULL
};

PyObject *
PyInit_tailseqext(void)
{
    PyObject *m;

    import_array();

    /* Create the module and add the functions */
    m = PyModule_Create(&tailseqext_def);
    if (m == NULL)
        return NULL;

    /* Add some symbolic constants to the module */
    struct tailseqext_state *st = GETSTATE(m);

    st->error = PyErr_NewException("tailseqext.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
