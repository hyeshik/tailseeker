/**
 * Assorted accelerator functions for Tail-seq from Narry Kim Lab.
 *
 * Hyeshik Chang <hyeshik@snu.ac.kr>
 *
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <numpy/arrayobject.h>

static PyObject *ErrorObject;

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
    PyObject *r;
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

    if ((r = PyArray_SimpleNew(2, dim, NPY_INT16)) == NULL)
        return NULL;

    if (decode_intensity(encoded, dim[0], PyArray_DATA(r)) == -1)
        Py_XDECREF(r);

    return r;
}


/* List of functions defined in the module */

static PyMethodDef tailseqext_methods[] = {
    {"decode_intensity", Py_decode_intensity, METH_VARARGS,
        PyDoc_STR("decode_intensity(str) -> array")},
    {NULL, NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
"This is a template module just for instruction.");

/* Initialization function for the module (*must* be called initxx) */

PyMODINIT_FUNC
inittailseqext2(void)
{
    PyObject *m;

    import_array();

    /* Create the module and add the functions */
    m = Py_InitModule3("tailseqext2", tailseqext_methods, module_doc);
    if (m == NULL)
        return;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("tailseqext2.error", NULL, NULL);
        if (ErrorObject == NULL)
            return;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
}
