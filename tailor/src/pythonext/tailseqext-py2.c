/**
 * Assorted accelerator functions for Tail-seq from Narry Kim Lab.
 *
 * Hyeshik Chang <hyeshik@snu.ac.kr>
 *
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"

static PyObject *ErrorObject;

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

        if (!PyInt_Check(value)) {
            PyErr_SetString(PyExc_TypeError, "Values of weight map must be an integer");
            Py_DECREF(self);
            return NULL;
        }

        self->weight[(int)c] = (int)PyInt_AsLong(value);
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
    {"PolyALocator", PolyALocator_new, METH_VARARGS,
        PyDoc_STR("PolyALocator(weights) -> object")},
    {NULL, NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
"This is a template module just for instruction.");

/* Initialization function for the module (*must* be called initxx) */

PyMODINIT_FUNC
inittailseqext(void)
{
    PyObject *m;

    /* Create the module and add the functions */
    m = Py_InitModule3("tailseqext", tailseqext_methods, module_doc);
    if (m == NULL)
        return;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("tailseqext.error", NULL, NULL);
        if (ErrorObject == NULL)
            return;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
}
