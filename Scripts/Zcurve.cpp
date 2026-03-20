#include <Python.h>
#include <numpy/arrayobject.h>
#include <thread>
#include <vector>
// memory allocation
#define  NEW    new (std::nothrow)
// offset constants
#define  X      0
#define  Y      1
#define  Z      2
#define  PHASE  3
// float constants
#define  V      1.0F/2
#define  W      1.0F/3
#define  Q      1.0F/4
// Z-curve transformation dimension
#define  DIM_A  765
/* 
 * Map for converting ASCII chars into one-hot vectors
 *
 * A = [1, 0, 0, 0] G = [0, 1, 0, 0]
 * C = [0, 0, 1, 0] T = [0, 0, 0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float ONE_HOT[][4] = 
{   
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0},  
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0},
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0}, 
    {0, 0, V, V}, {1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0}, 
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0},
    {0, 0, V, V}, {1, 0, 0, 0}
};
/* 
 * Map for converting ASCII chars into Z-curve coordinates
 *
 * A = [+1, +1, +1]  G = [+1, -1, -1]
 * C = [-1, +1, -1]  T = [-1, -1, +1]
 * 
 * Degenerate symbols are calculated by weighted vector sum
 */
static float Z_COORD[][3] = {
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}
};
/**
 * @brief           Calculate 1-mer Z-curve params for a given sequence.
 * @param seq       The input sequence.
 * @param len       The length of the input sequence.
 * @param params    The output Z-curve parameters.
 */
static void mono_trans(const char *seq, int len, float *params) {
    float counts[PHASE][3] = {{0.0f}};
        
    for (int i = 0; i < len; i ++) {
        int p = i % PHASE;

        counts[p][X] += Z_COORD[(int)seq[i]][X];
        counts[p][Y] += Z_COORD[(int)seq[i]][Y];
        counts[p][Z] += Z_COORD[(int)seq[i]][Z];
    }
        
    for (int p = 0; p < PHASE; p ++, params += 3) {
        params[X] = counts[p][X] / len * PHASE;
        params[Y] = counts[p][Y] / len * PHASE;
        params[Z] = counts[p][Z] / len * PHASE;
    }
}
/**
 * @brief           Calculate 2-mer Z-curve params for a given sequence.
 * @param seq       The input sequence.
 * @param len       The length of the input sequence.
 * @param params    The output Z-curve parameters.
 */
static void di_trans(const char *seq, int len, float *params) {
    float counts[PHASE][4][3] = {{{0.0f}}};
    int sublen = len - 1;

    for (int i = 0; i < sublen; i++) {
        int p = i % PHASE;
        for (int b = 0; b < 4; b ++) {
            counts[p][b][X] += ONE_HOT[(int)seq[i]][b] * Z_COORD[(int)seq[i + 1]][X];
            counts[p][b][Y] += ONE_HOT[(int)seq[i]][b] * Z_COORD[(int)seq[i + 1]][Y];
            counts[p][b][Z] += ONE_HOT[(int)seq[i]][b] * Z_COORD[(int)seq[i + 1]][Z];
        }
    }
        
    for (int p = 0; p < PHASE; p ++)
        for (int b = 0; b < 4; b ++, params += 3) {
            params[X] = counts[p][b][X] / len * PHASE;
            params[Y] = counts[p][b][Y] / len * PHASE;
            params[Z] = counts[p][b][Z] / len * PHASE;
        }
}
/**
 * @brief           Calculate 3-mer Z-curve params for a given sequence.
 * @param seq       The input sequence.
 * @param len       The length of the input sequence.
 * @param params    The output Z-curve parameters.
 */
static void tri_trans(const char *seq, int len, float *params) {
    float counts[PHASE][4][4][3] = {{{{0.0f}}}};
    int sublen = len - 2;

    for (int i = 0; i < sublen; i ++) {
        int p = i % PHASE;
        for (int s = 0; s < 4; s ++)
        for (int b = 0; b < 4; b ++) {
            counts[p][s][b][X] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i + 1]][b] * 
                                  Z_COORD[(int)seq[i + 2]][X];
            counts[p][s][b][Y] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i + 1]][b] * 
                                  Z_COORD[(int)seq[i + 2]][Y];
            counts[p][s][b][Z] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i + 1]][b] * 
                                  Z_COORD[(int)seq[i + 2]][Z];
        }
    }

    for (int p = 0; p < PHASE; p ++)
    for (int s = 0; s < 4; s ++)
    for (int b = 0; b < 4; b ++, params += 3) {
        params[X] = counts[p][s][b][X] / len * PHASE;
        params[Y] = counts[p][s][b][Y] / len * PHASE;
        params[Z] = counts[p][s][b][Z] / len * PHASE;
    }
}
/**
 * @brief           Calculate 3-mer Z-curve params for a given sequence.
 * @param seq       The input sequence.
 * @param len       The length of the input sequence.
 * @param params    The output Z-curve parameters.
 */
static void quart_trans(const char *seq, int len, float *params) {
    float counts[PHASE][4][4][4][3] = {{{{{0.0f}}}}};
    int sublen = len - 2;

    for (int i = 0; i < sublen; i ++) {
        int p = i % PHASE;
        for (int s = 0; s < 4; s ++)
        for (int t = 0; t < 4; t ++)
        for (int b = 0; b < 4; b ++) {
            counts[p][s][t][b][X] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i]][t] * 
                                     ONE_HOT[(int)seq[i + 1]][b] * Z_COORD[(int)seq[i + 2]][X];
            counts[p][s][t][b][Y] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i]][t] * 
                                     ONE_HOT[(int)seq[i + 1]][b] * Z_COORD[(int)seq[i + 2]][Y];
            counts[p][s][t][b][Z] += ONE_HOT[(int)seq[i]][s] * ONE_HOT[(int)seq[i]][t] * 
                                     ONE_HOT[(int)seq[i + 1]][b] * Z_COORD[(int)seq[i + 2]][Z];
        }
    }

    for (int p = 0; p < PHASE; p ++)
    for (int s = 0; s < 4; s ++)
    for (int t = 0; t < 4; t ++)
    for (int b = 0; b < 4; b ++, params += 3) {
        params[X] = counts[p][s][t][b][X] / len * PHASE;
        params[Y] = counts[p][s][t][b][Y] / len * PHASE;
        params[Z] = counts[p][s][t][b][Z] / len * PHASE;
    }
}

PyObject *encode(PyObject *self, PyObject *args, PyObject *kw) {
    import_array();
    static char *kwlist[] = {(char *)"records", (char *)"n_jobs", NULL};
    PyObject *records = NULL;
    int n_jobs = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "O|i", kwlist, &records, &n_jobs))
        return NULL;
    
    // Check if n_jobs is valid
    if (n_jobs <= 0) n_jobs = std::thread::hardware_concurrency();

    // Check if records is List[str]
    if (!PyList_Check(records)) {
        PyErr_SetString(PyExc_TypeError, "records must be List[str]");
        return NULL;
    }
    int n_records = (int) PyList_Size(records);
    
    // Prepare input sequences
    std::vector<char *> seqs(n_records);
    std::vector<int> lens(n_records);
    for (int i = 0; i < n_records; i++) {
        PyObject *item = PyList_GetItem(records, i);
        if (!PyUnicode_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "records must be List[str]");
            return NULL;
        }
        lens[i] = (int) PyUnicode_GetLength(item);
        seqs[i] = (char *)PyUnicode_AsUTF8(item);
    }

    float *params = new float[n_records * DIM_A];
    if (params == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }
    std::thread **threads = new std::thread *[n_jobs];
    if (threads == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }
    // convert sequence to Z-curve params
    for (int i = 0; i < n_jobs; i++) {
        threads[i] = NEW std::thread([i, n_jobs, n_records, params, &seqs, &lens]() {
            for (int j = i; j < n_records; j += n_jobs) {
                float *head = params + (j * DIM_A);
                mono_trans(seqs[j], lens[j], head);
                di_trans(seqs[j], lens[j], head+9);
                tri_trans(seqs[j], lens[j], head+45);
                quart_trans(seqs[j], lens[j], head+189);
            }
        });
    }
    for (int i = 0; i < n_jobs; i ++) {
        if (threads[i] == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
            return NULL;
        }
        threads[i]->join();
        delete threads[i];
    }
    delete[] threads;

    npy_intp dims[] = { n_records, DIM_A };
    PyObject* np_array = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, params);

    if (!np_array) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create NumPy array");
        delete[] params;
        return NULL;
    }

    PyArray_ENABLEFLAGS((PyArrayObject*)np_array, NPY_ARRAY_OWNDATA);
    
    return np_array;
}

PyMODINIT_FUNC PyInit_Zcurve(void) {
    import_array();
    static PyMethodDef methods[] = {
        {"encode", (PyCFunction)encode, METH_VARARGS | METH_KEYWORDS, NULL},
        {NULL, NULL, 0, NULL}
    };
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "Zcurve",
        NULL,
        -1,
        methods
    };
    return PyModule_Create(&moduledef);
}