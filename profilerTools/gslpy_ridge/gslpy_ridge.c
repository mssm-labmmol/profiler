#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <Python.h>
#include <numpy/arrayobject.h>

void mark_debug (void)
{
    static int number;
    printf("---> Mark %d <---\n", number);
    fflush(stdout);
    number++;
}


void label_debug (const char * label)
{
    printf("---> Mark %s <---\n", label);
    fflush(stdout);
}

/* =========================== Pure-C source code =========================== */

/* Regularization workspace */
struct rglz_wksp {
    double lamb; /* Lambda value */
    gsl_vector *Ldiag; /* Pointer to vector with diagonal values of the
			      * regularization matrix L. */
    const gsl_vector *Y; /* Pointer to original Y */
    gsl_vector *Ys; /* Ponter to Y vector of the standard form. */
    const gsl_matrix *X; /* Pointer to original X matrix  */
    gsl_matrix *Xs; /* Pointer to X matrix of the standard form. */
    gsl_multifit_linear_workspace *lwksp; /* Pointer to GSL workspace where SVDs are stored. */
    /* L-curve data */
    struct lcurve {
	gsl_vector *reg_param;
	gsl_vector *rho;
	gsl_vector *eta;
    } lcurve;
};
typedef struct rglz_wksp rglz_wksp;


/* Allocate and initialize inner variables of a regularization workspace. */
rglz_wksp * rglz_init(const double *Ldiag, const int p, 
		      const gsl_vector *Y, const gsl_matrix *X, const gsl_vector *W,
		      gsl_multifit_linear_workspace *lwksp)
{
    rglz_wksp *rw;

    /* Allocate space for the workspace */
    rw = (rglz_wksp*) malloc(sizeof(rglz_wksp));

    rw->lamb = 0.0; /* Initial value */
    rw->lcurve.reg_param = NULL;
    rw->lcurve.rho = NULL;
    rw->lcurve.eta = NULL;

    /* Allocate space for the standard-form matrixes */
    rw->Xs = gsl_matrix_alloc(X->size1, X->size2);
    rw->Ys = gsl_vector_alloc(Y->size);

    /* Initialize the regularization weights Ldiag */
    rw->Ldiag = gsl_vector_alloc(p);
    for (int i = 0; i < p; i++)
    {
	gsl_vector_set(rw->Ldiag, i, Ldiag[i]);
    }

    /* Convert the system to standard form */
    gsl_multifit_linear_wstdform1(rw->Ldiag, X, W, Y, rw->Xs, rw->Ys, lwksp);

    /* Calculate SVD of standard-form matrix Xs */
    gsl_multifit_linear_svd(rw->Xs, lwksp);

    rw->lwksp = lwksp;
    rw->X = X;
    rw->Y = Y;
    return rw;
}

void rglz_estimate_lambda(rglz_wksp *rw, const int lcsize)
{
    size_t idx; /* Temporary index to store corner of the L-curve */
    /* Allocate vectors */
    rw->lcurve.reg_param = gsl_vector_alloc(lcsize);
    rw->lcurve.rho = gsl_vector_alloc(lcsize);
    rw->lcurve.eta = gsl_vector_alloc(lcsize);
    /* Create Lcurve */
    gsl_multifit_linear_lcurve(rw->Y, rw->lcurve.reg_param, rw->lcurve.rho, rw->lcurve.eta, rw->lwksp);
    /* Estimate corner */
    gsl_multifit_linear_lcorner(rw->lcurve.rho, rw->lcurve.eta, &idx);
    /* Set lambda estimate */
    rw->lamb = gsl_vector_get(rw->lcurve.reg_param, idx);
}

void rglz_solve(gsl_vector *c, double *rnorm, double *snorm, rglz_wksp *rw)
{
    /* Solve in standard form */
    gsl_multifit_linear_solve(rw->lamb, rw->Xs, rw->Ys, c, rnorm, snorm, rw->lwksp);
    /* Convert back to original parameters */
    gsl_multifit_linear_genform1(rw->Ldiag, c, c, rw->lwksp);
    return;
}

inline void rglz_set_lambda(rglz_wksp *rw, const int lamb)
{
    rw->lamb = lamb;
}

void rglz_close(rglz_wksp *rw)
{
    gsl_vector_free(rw->Ldiag);
    gsl_matrix_free(rw->Xs);
    gsl_vector_free(rw->Ys);
    if (rw->lcurve.reg_param != NULL)
	gsl_vector_free(rw->lcurve.reg_param);
    if (rw->lcurve.rho != NULL)
	gsl_vector_free(rw->lcurve.rho);
    if (rw->lcurve.eta != NULL)
	gsl_vector_free(rw->lcurve.eta);
    free(rw);
}

int
gslpy_ridge_regression(const gsl_vector * y,
		       const gsl_matrix * X,
		       const gsl_vector * w,
		       const double lambda, 
		       const double * L,
		       const unsigned long Lcurve_size,
		       gsl_vector * coefs)
{

    
    
    unsigned long n_samples = X->size1;
    unsigned long n_features = X->size2;

    double snorm;
    double chisq;

    gsl_multifit_linear_workspace
	* gmlw = gsl_multifit_linear_alloc(n_samples, n_features);
    
    rglz_wksp * rw = rglz_init(L, n_features, y, X, w, gmlw);

    /* If the user specified gamma, do not attempt to optimize it. */
    if (lambda >= 0)
	rglz_set_lambda(rw, lambda);
    else
	rglz_estimate_lambda(rw, Lcurve_size);

    rglz_solve(coefs, &chisq, &snorm, rw);

    rglz_close(rw);

    gsl_multifit_linear_free(gmlw);
    return 0;
}

/* ============================= Python wrapper ============================= */

static PyObject*
gslpy_ridge_fit(PyObject* self, PyObject* args, PyObject* kwargs)
{
    static char * keywords[] = {
	"X",
	"y",
	"sample_weight",
	"lamb",
	"ldiag",
	"lcurve_size",
	NULL,
    };

    PyArrayObject
	*X=NULL,
	*y=NULL,
	*w=NULL,
	*Ldiag=NULL,
	*coefs=NULL;

    /* Negative lambda means that lambda is optimized */
    double lambda = -1;
    unsigned long Lcurve_size = 10000;

    /* Flag to check if things are set */
    int wei_set = 0;
    int L_set = 0;

    /* GSL stuff */
    gsl_vector * g_y;
    gsl_matrix * g_X;
    gsl_vector * g_w;
    gsl_vector * g_coefs;
    
    if (!PyArg_ParseTupleAndKeywords(args,
                                     kwargs,
                                     "O!O!|O!dO!k",
                                     keywords,
                                     &PyArray_Type, &X,
				     &PyArray_Type, &y,
				     &PyArray_Type, &w,
				     &lambda,
				     &PyArray_Type, &Ldiag,
				     &Lcurve_size)) {
        return NULL;
    }

    
    /* If weights are not given, set them as 1 */
    if (w == NULL)
    {
	wei_set = 1;
	w = (PyArrayObject*) PyArray_ZEROS(PyArray_NDIM(y),
					   PyArray_DIMS(y),
					   NPY_DOUBLE,
					   0);
	
	npy_double * w_data = PyArray_DATA(w);
	
	for (npy_intp i = 0; i < PyArray_SIZE(w); ++i)
	    w_data[i] = 1.0;
    }

    /* If Ldiag is not given, set it to identity. */
    if (Ldiag == NULL)
    {
	npy_intp Lshape[] = {PyArray_DIM(X,1)};
	L_set = 1;
	Ldiag = (PyArrayObject*) PyArray_ZEROS(1,
					       Lshape,
					       NPY_DOUBLE,
					       0);

	for (npy_intp i = 0; i < PyArray_DIM(X, 1); i++)
	    ((double *) PyArray_DATA(Ldiag))[i] = 1.0;
    }
    

    /* Convert PyArray objects to gsl types */
    g_y = gsl_vector_calloc(PyArray_DIM(y, 0));
    g_w = gsl_vector_calloc(PyArray_DIM(w, 0));
    g_X = gsl_matrix_calloc(PyArray_DIMS(X)[0], PyArray_DIMS(X)[1]);
    g_coefs = gsl_vector_calloc(PyArray_DIM(X, 1));
    
    for (npy_intp i = 0; i < (npy_intp) g_y->size; i++)
    {
	gsl_vector_set(g_y, i, ((double *) PyArray_DATA(y))[i]);
	gsl_vector_set(g_w, i, ((double *) PyArray_DATA(w))[i]);
    }
    
    for (npy_intp i = 0; i < PyArray_SIZE(X); i++)
	g_X->data[i] = ((double *) PyArray_DATA(X))[i];

    /* Run regression */
    gslpy_ridge_regression(g_y, g_X,
			   g_w, lambda,
			   (double *) PyArray_DATA(Ldiag),
			   Lcurve_size, g_coefs);
    

    /* Convert g_coefs back to a Python object */
    /* Local scope just to delete ``coefs_shape``after. */
    {
	npy_intp coefs_shape[] = {g_coefs->size, 1};
	coefs = (PyArrayObject *) PyArray_ZEROS(2,
						coefs_shape,
						NPY_DOUBLE,
						0);
	for (npy_intp i = 0; i < (npy_intp) g_coefs->size; i++)
	    ((double *) PyArray_DATA(coefs))[i] = gsl_vector_get(g_coefs, i);

    }
    
    if (wei_set)
	Py_XDECREF(w);

    if (L_set)
	Py_XDECREF(Ldiag);

    gsl_matrix_free(g_X);
    gsl_vector_free(g_y);
    gsl_vector_free(g_w);
    gsl_vector_free(g_coefs);

    return (PyObject *) coefs;
}

/* =========================== Python boilerplate =========================== */

PyDoc_STRVAR(module_doc,
	     "Ridge regression model methods using the GNU Scientific Library.");

PyDoc_STRVAR(fit_doc,
	     "Fits a ridge regression model to data.\n"
	     "fit(X, y, sample_weights=None, lambda=None, ldiag=None, Lcurve_size=10000).\n");


PyMethodDef methods[] = {
    {
	"fit",
	(PyCFunction) gslpy_ridge_fit,
	METH_VARARGS | METH_KEYWORDS,
	fit_doc,
    },
    {NULL},
};


PyModuleDef gslpyridge_module = {
    PyModuleDef_HEAD_INIT,
    "gslpyridge",
    module_doc,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_gslpyridge(void)
{
    import_array();
    return PyModule_Create(&gslpyridge_module);
}
