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

#define NCHANNELS   4

#define SIGNAL_MIN  -255
#define CODED_MAX   4095

static const char *base64_enc_table="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"
                                    "ghijklmnopqrstuvwxyz0123456789+/";

#define NE          64
static const int16_t base64_dec_table[256]={
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, 62, NE, NE, NE, 63,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, NE, NE, NE, NE, NE, NE,
    NE,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, NE, NE, NE, NE, NE,
    NE, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE,
    NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE, NE
};


static int
decode_intensity(const char *encoded, int length, int16_t *signals)
{
    int chan, cycle;

    for (cycle = 0; cycle < length; cycle++) {
        for (chan = 0; chan < NCHANNELS; chan++) {
            int16_t high, low;

            high = base64_dec_table[(unsigned char)*(encoded++)];
            low = base64_dec_table[(unsigned char)*(encoded++)];

            if (high == NE || low == NE) {
                PyErr_SetString(PyExc_ValueError, "Illegal character found.");
                return -1;
            }

            *(signals++) = ((high << 6) | low) + SIGNAL_MIN;
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

static void
encode_intensity(int length, const int16_t *signals, char *outbuf)
{
    int chan, cycle;
    int16_t v;

    for (cycle = 0; cycle < length; cycle++) {
        for (chan = 0; chan < NCHANNELS; chan++) {
            v = *(signals++) - SIGNAL_MIN;
            v = (v <= CODED_MAX ? v : CODED_MAX);
            v = (v >= 0 ? v : 0);
            *(outbuf++) = base64_enc_table[v >> 6];
            *(outbuf++) = base64_enc_table[v & 63];
        }
    }
}

static PyObject *
Py_encode_intensity(PyObject *self, PyObject *args)
{
    PyObject *signals, *r;
    int length;

    if (!PyArg_ParseTuple(args, "O:encode_intensity", &signals))
        return NULL;

    if (!PyArray_Check(signals)) {
        PyErr_SetString(PyExc_TypeError, "needs a numpy array.");
        return NULL;
    }

    if (PyArray_NDIM(signals) != 2 || PyArray_DIMS(signals)[1] != NCHANNELS ||
            PyArray_DESCR(signals) != PyArray_DescrFromType(NPY_INT16)) {
        PyErr_Format(PyExc_TypeError, "array must be N x %d int16.", NCHANNELS);
        return NULL;
    }

    length = PyArray_DIMS(signals)[0];
    if ((r = PyString_FromStringAndSize(NULL, length * 8)) == NULL)
        return NULL;

    encode_intensity(length, PyArray_DATA(signals), PyString_AS_STRING(r));
    return r;
}

static int
encode_intensity_from_string(const char *intxt, int ncycles, char *outbuf)
{
    int chan, cycle;
    int16_t v;

    for (cycle = 0; cycle < ncycles; cycle++) {
        int readchars, intensity[NCHANNELS];

        /* The next line is hard coded with fixed number of channels. Change it when needed */
#if NCHANNELS == 4
        if (sscanf(intxt, "%d %d %d %d%n", intensity, intensity+1, intensity+2, intensity+3,
                    &readchars) < NCHANNELS)
#endif
            return -1;

        intxt += readchars;

        for (chan = 0; chan < NCHANNELS; chan++) {
            v = intensity[chan] - SIGNAL_MIN;
            v = (v <= CODED_MAX ? v : CODED_MAX);
            v = (v >= 0 ? v : 0);
            *(outbuf++) = base64_enc_table[v >> 6];
            *(outbuf++) = base64_enc_table[v & 63];
        }
    }

    return 0;
}

static PyObject *
Py_encode_intensity_from_string(PyObject *self, PyObject *args)
{
    PyObject *r;
    int length;
    char *instr, *outstr;

    if (!PyArg_ParseTuple(args, "si:encode_intensity_from_string", &instr, &length))
        return NULL;

    r = PyString_FromStringAndSize(NULL, length * NCHANNELS * 2);
    if (r == NULL)
        return NULL;

    outstr = PyString_AS_STRING(r);
    outstr[length * NCHANNELS * 2] = 0;

    if (encode_intensity_from_string(instr, length, outstr)) {
        PyErr_SetString(PyExc_ValueError, "Couldn't parse the string.");
        return NULL;
    }

    return r;
}

static PyObject *
Py_phred_64to33(PyObject *self, PyObject *args)
{
    PyObject *r;
    Py_ssize_t length;
    char *instr, *outstr;

    if (!PyArg_ParseTuple(args, "s#:phred_64to33", &instr, &length))
        return NULL;

    r = PyString_FromStringAndSize(NULL, length);
    if (r == NULL)
        return NULL;

    outstr = PyString_AS_STRING(r);
    outstr[length] = 0;

    while (*instr != 0)
        *(outstr++) = *(instr++) - 31;

    return r;
}


/* ======================================================== */
/* decode_index */

#define MAX_INDICES         128
#define MAX_INDEX_LENGTH    8

typedef struct {
    PyObject_HEAD
    int indexlen;
    int numindices;
    char indices[MAX_INDICES][MAX_INDEX_LENGTH];
} IndexDecoderObject;

static PyTypeObject IndexDecoder_Type;

#define IndexDecoder_Type(v)        (Py_TYPE(v) == &IndexDecoder_Type)

static IndexDecoderObject *
IndexDecoder_create(PyObject *args)
{
    IndexDecoderObject *self;
    PyObject *indexlist=NULL;
    Py_ssize_t i;

    if (!PyArg_ParseTuple(args, "O:IndexDecoder", &indexlist))
        return NULL;

    if (!PyList_Check(indexlist)) {
        PyErr_SetString(PyExc_TypeError, "List of index must be given in a list.");
        return NULL;
    }

    self = PyObject_New(IndexDecoderObject, &IndexDecoder_Type);
    if (self == NULL)
        return NULL;

    self->numindices = PyList_Size(indexlist);
    self->indexlen = 0;

    if (self->numindices > MAX_INDICES) {
        PyErr_SetString(PyExc_TypeError, "too many indices.");
        Py_DECREF(self);
        return NULL;
    }

    for (i = 0; i < self->numindices; i++) {
        PyObject *item;

        item = PyList_GetItem(indexlist, i);
        if (!PyString_Check(item)) {
            PyErr_Format(PyExc_TypeError, "%d-th index is not a string.", (int)i + 1);
            Py_DECREF(self);
            return NULL;
        }

        if (i == 0) {
            self->indexlen = PyString_Size(item);

            if (self->indexlen > MAX_INDEX_LENGTH) {
                PyErr_SetString(PyExc_TypeError, "index is too long.");
                Py_DECREF(self);
                return NULL;
            }
        }
        else if (PyString_Size(item) != self->indexlen) {
            PyErr_Format(PyExc_TypeError, "length of %d-th index is not %d.", (int)i + 1,
                         (int)self->indexlen);
            Py_DECREF(self);
            return NULL;
        }

        memcpy(self->indices[i], PyString_AS_STRING(item), self->indexlen);
    }

    return self;
}

static void
IndexDecoder_dealloc(IndexDecoderObject *self)
{
    PyObject_Del(self);
}

static PyObject *
IndexDecoder_call(IndexDecoderObject *self, PyObject *args, PyObject *kw)
{
    char *inseq;
    Py_ssize_t inseqlen;
    int bestidx, bestmismatches, i, secondbestfound;

    if (!PyArg_ParseTuple(args, "s#:decode", &inseq, &inseqlen))
        return NULL;

    if (inseqlen != self->indexlen) {
        PyErr_SetString(PyExc_TypeError, "index size mismatch.");
        return NULL;
    }

    bestidx = -1;
    bestmismatches = self->indexlen + 1;
    secondbestfound = 0;

    for (i = 0; i < self->numindices; i++) {
        int j, mismatches=0;

        for (j = 0; j < self->indexlen; j++)
            mismatches += (inseq[j] != self->indices[i][j]);

        if (mismatches < bestmismatches) {
            bestidx = i;
            bestmismatches = mismatches;
            secondbestfound = 0;
        }
        else if (mismatches == bestmismatches)
            secondbestfound = 1;
    }

    if (bestidx == -1) {
        PyErr_SetString(PyExc_ValueError, "no feasible index matched.");
        return NULL;
    }

    if (secondbestfound == 1) {
        PyErr_SetString(PyExc_ValueError, "multiple indices have same distance.");
        return NULL;
    }

    return Py_BuildValue("(is#)", bestmismatches, self->indices[bestidx], self->indexlen);
}

static PyTypeObject IndexDecoder_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tailseqext.IndexDecoder",           /*tp_name*/
    sizeof(IndexDecoderObject),          /*tp_basicsize*/
    0,                      /*tp_itemsize*/
    /* methods */
    (destructor)IndexDecoder_dealloc, /*tp_dealloc*/
    0,                      /*tp_print*/
    0,                      /*tp_getattr*/
    0,                      /*tp_setattr*/
    0,                      /*tp_compare*/
    0,                      /*tp_repr*/
    0,                      /*tp_as_number*/
    0,                      /*tp_as_sequence*/
    0,                      /*tp_as_mapping*/
    0,                      /*tp_hash*/
    (ternaryfunc)IndexDecoder_call, /*tp_call*/
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
IndexDecoder_new(PyObject *self, PyObject *args)
{
    return (PyObject *)IndexDecoder_create(args);
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
    {"decode_intensity", Py_decode_intensity, METH_VARARGS,
        PyDoc_STR("decode_intensity(str) -> array")},
    {"encode_intensity", Py_encode_intensity, METH_VARARGS,
        PyDoc_STR("encode_intensity(array) -> str")},
    {"encode_intensity_from_string", Py_encode_intensity_from_string, METH_VARARGS,
        PyDoc_STR("encode_intensity_from_string(str) -> str")},
    {"phred_64to33", Py_phred_64to33, METH_VARARGS,
        PyDoc_STR("phred_64to33(str) -> str")},
    {"IndexDecoder", IndexDecoder_new, METH_VARARGS,
        PyDoc_STR("IndexDecoder(indices) -> object")},
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

    import_array();

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
