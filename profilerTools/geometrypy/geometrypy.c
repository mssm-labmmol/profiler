#define PY_SSIZE_T_CLEAN
#define GPY_DIMENSIONS 3

#include <Python.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <numpy/arrayobject.h>

inline int nan_check (const double x)
{
    return (x != x);
}

void * NanTracker (const char *label, const int n, ...)
{
    va_list ap;
    va_start (ap, n);
    for (int i = 0; i < n; i++)
    {
	double number;
	number = va_arg (ap, double);
	// check NaN and raise errror
	if (nan_check(number))
	{
	    char msg[64];
	    sprintf(msg, "NaN at label %s, argument %d.\n", label, i+1);
	    PyErr_SetString (PyExc_ValueError, msg);
	    return NULL;
	}
    }
    va_end (ap);
    return (void *) label;
}

PyArrayObject *PyArray_Einsum_i_j_ij (PyArrayObject *a1, PyArrayObject *a2) {
  PyArrayObject *output=NULL;
  npy_intp ndims, ndims_other, *dims=NULL, *dims_other=NULL, dims_out[]={0, 0};
  npy_double *a1_c=NULL, *a2_c=NULL, *output_c=NULL;
  ndims = PyArray_NDIM(a1);
  ndims_other = PyArray_NDIM(a2);
  dims = PyArray_DIMS(a1);
  dims_other = PyArray_DIMS(a2);
  dims_out[0] = dims[0];
  dims_out[1] = dims_other[0];
  if ((ndims != ndims_other) ||
      (ndims != 1))
    goto badargument;
  output = (PyArrayObject*) PyArray_EMPTY(2, dims_out, NPY_DOUBLE, 0);
  output_c = PyArray_DATA(output);
  a1_c = PyArray_DATA(a1);
  a2_c = PyArray_DATA(a2);
  for (int i = 0; i < dims_out[0]; i++) {
    for (int j = 0; j < dims_out[1]; j++) {
      output_c[dims_out[1]*i + j] = a1_c[i] * a2_c[j];
      //NanTracker("Einsum_i_j_ij", 1, output_c[dims_out[1]*i + j]);
    }
  }
  return output;
 badargument:
  PyErr_SetString(PyExc_ValueError, "Einsum i,j->ij expects two 1-D PyArrays.");
  return NULL;
}

PyArrayObject *PyArray_Einsum_ij_ij_i (PyArrayObject *a1, PyArrayObject *a2) {
  PyArrayObject *output=NULL;
  npy_intp ndims, ndims_other, *dims=NULL, *dims_other=NULL, dims_out[]={0};
  npy_double *a1_c=NULL, *a2_c=NULL, *output_c=NULL;
  ndims = PyArray_NDIM(a1);
  ndims_other = PyArray_NDIM(a2);
  dims = PyArray_DIMS(a1);
  dims_other = PyArray_DIMS(a2);
  dims_out[0] = dims[0];
  if ((ndims != ndims_other) ||
      (ndims != 2))
    goto badargument;
  for (int i = 0; i < ndims; i++) {
    if (dims[i] != dims_other[i]) {
      goto badargument;
    }
  }
  output = (PyArrayObject*) PyArray_SimpleNew(1, dims_out, NPY_DOUBLE);
  output_c = PyArray_DATA(output);
  a1_c = PyArray_DATA(a1);
  a2_c = PyArray_DATA(a2);
  for (int i = 0; i < dims[0]; i++) {
    output_c[i] = 0.0;
    for (int j = 0; j < dims[1]; j++) {
      output_c[i] += a1_c[dims[1]*i + j] * a2_c[dims[1]*i + j];
      //NanTracker("Einsum_ij_ij_i", 1, output_c[i]);
    }
  }
  return output;
 badargument:
  PyErr_SetString(PyExc_ValueError,
     "Einsum ij,ij->i expects two 2-D PyArrays with the same shape.");
  return NULL;
}

PyArrayObject *PyArray_Einsum_ij_ij_j (PyArrayObject *a1, PyArrayObject *a2) {
  PyArrayObject *output=NULL;
  npy_intp ndims, ndims_other, *dims=NULL, *dims_other=NULL, dims_out[]={0};
  npy_double *a1_c=NULL, *a2_c=NULL, *output_c=NULL;
  ndims = PyArray_NDIM(a1);
  ndims_other = PyArray_NDIM(a2);
  dims = PyArray_DIMS(a1);
  dims_other = PyArray_DIMS(a2);
  dims_out[0] = dims[1];
  if ((ndims != ndims_other) ||
      (ndims != 2))
    goto badargument;
  for (int i = 0; i < ndims; i++) {
    if (dims[i] != dims_other[i]) {
      goto badargument;
    }
  }
  output = (PyArrayObject*) PyArray_SimpleNew(1, dims_out, NPY_DOUBLE);
  output_c = PyArray_DATA(output);
  a1_c = PyArray_DATA(a1);
  a2_c = PyArray_DATA(a2);
  for (int i = 0; i < dims[1]; i++) {
    output_c[i] = 0.0;
    for (int j = 0; j < dims[0]; j++) {
      output_c[i] += a1_c[dims[1]*j + i] * a2_c[dims[1]*j + i];
      //NanTracker("Einsum_ij_ij_j", 1, output_c[i]);
    }
  }
  return output;
 badargument:
  PyErr_SetString(PyExc_ValueError,
    "Einsum ij,ij->j expects two 2-D PyArrays with the same shape.");
  return NULL;
}

PyArrayObject *PyArray_Einsum_ij_kj_kj (PyArrayObject *a1, PyArrayObject *a2) {
  PyArrayObject *output=NULL;
  npy_intp ndims, ndims_other, *dims=NULL, *dims_other=NULL, *dims_out=NULL;
  npy_double *a1_c=NULL, *a2_c=NULL, *output_c=NULL;
  ndims = PyArray_NDIM(a1);
  ndims_other = PyArray_NDIM(a2);
  dims = PyArray_DIMS(a1);
  dims_other = PyArray_DIMS(a2);
  dims_out = dims_other;
  if ((ndims != ndims_other) ||
      (ndims != 2))
    goto badargument;
  /* only second dimensions need to match */
  if (dims[1] != dims_other[1])
    goto badargument;
  output = (PyArrayObject*) PyArray_SimpleNew(2, dims_out, NPY_DOUBLE);
  output_c = PyArray_DATA(output);
  a1_c = PyArray_DATA(a1);
  a2_c = PyArray_DATA(a2);
  for (int k = 0; k < dims_out[0]; k++) {
    for (int j = 0; j < dims_out[1]; j++) {
      output_c[dims_out[1]*k + j] = 0.0;
      for (int i = 0; i < dims[0]; i++) {
	output_c[dims_out[1]*k + j] +=
	    a1_c[dims[1]*i + j] * a2_c[dims_out[1]*k + j];
	//NanTracker("Einsum_ij_kj_kj", 1, output_c[dims_out[1]*k + j]);
      }
    }
  }
  return output;
 badargument:
  PyErr_SetString(PyExc_ValueError,
    "Einsum ij,kj->kj expects two 2-D PyArrays with the same number of columns.");
  return NULL;
}

typedef struct t_einsum {
  const char *key; /* einsum string */
  PyArrayObject (*(*einsum_func) (PyArrayObject*, PyArrayObject*)); /* pointer
								     * to
								     * einsum
								     * function */
} t_einsum;

static t_einsum einsum_lookup_table[] = {
					 {"ij,ij->i", PyArray_Einsum_ij_ij_i},
					 {"ij,ij->j", PyArray_Einsum_ij_ij_j},
					 {"ij,kj->kj", PyArray_Einsum_ij_kj_kj},
					 {"i,j->ij", PyArray_Einsum_i_j_ij},
					 {NULL, NULL} /* Sentinel */
};

int MyCArray_CrossProduct_InPlace(npy_double *v1, npy_double *v2, npy_double *out) {
  if (GPY_DIMENSIONS != 3) {
    PyErr_BadArgument();
    return NULL;
  }
  out[0] = v1[1]*v2[2] - v1[2]*v2[1];
  out[1] = v1[2]*v2[0] - v1[0]*v2[2];
  out[2] = v1[0]*v2[1] - v1[1]*v2[0];
  //NanTracker("CrossProduct_InPlace", 3, out[0], out[1], out[2]);
  return 1;
}

/* Remember to FREE out after using it. */
npy_double* MyCArray_CrossProduct(npy_double *v1, npy_double *v2) {
  npy_double *out;
  out = (npy_double*) malloc (GPY_DIMENSIONS*sizeof(npy_double));
  if (out == NULL)
    return (void*) PyErr_NoMemory();
  if (!MyCArray_CrossProduct_InPlace(v1, v2, out))
    return NULL;
  return out;
}

/* Remember to FREE after use. */
npy_double* MyCArray_Subtract(npy_double *v1, npy_double* v2) {
  npy_double *out;
  if (GPY_DIMENSIONS != 3) {
    PyErr_BadArgument();
    return NULL;
  }
  out = (npy_double*) malloc (GPY_DIMENSIONS*sizeof(npy_double));
  if (out == NULL)
    return (void*) PyErr_NoMemory();
  out[0] = v1[0] - v2[0];
  out[1] = v1[1] - v2[1];
  out[2] = v1[2] - v2[2];
  //NanTracker("MyCArray_Subtract", 3, out[0], out[1], out[2]);
  return out;
}

int MyCArray_DotProduct(npy_double *v1, npy_double *v2, npy_double *dot) {
  if (GPY_DIMENSIONS != 3) {
    PyErr_BadArgument();
    return NULL;
  }
  *dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  //NanTracker("DotProduct", 1, *dot);
  return 1;
}

PyArrayObject* MyPyArray_Subtract(PyArrayObject* arr1, PyArrayObject* arr2) {
  PyArrayObject *output;
  output = (PyArrayObject*) PyArray_EnsureArray(PyNumber_Subtract(arr1, arr2));
  return output;
}

static PyObject *
einsumWrapper(PyObject *dummy, PyObject *args)
{
  char *einsum_string=NULL;
  PyObject *obj_1=NULL, *obj_2=NULL;
  PyArrayObject *a_1=NULL, *a_2=NULL;
  int lookup_i=0;
  
  if (!PyArg_ParseTuple(args, "sO!O!",
			&einsum_string,
			&PyArray_Type, &obj_1,
			&PyArray_Type, &obj_2))
    return NULL;

  a_1 = PyArray_FROM_OTF(obj_1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (a_1 == NULL) goto fail;
  a_2 = PyArray_FROM_OTF(obj_2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (a_2 == NULL) goto fail;

  while (einsum_lookup_table[lookup_i].key != NULL) {
    if (!strcmp(einsum_lookup_table[lookup_i].key, einsum_string)) {
      PyArrayObject *out = einsum_lookup_table[lookup_i].einsum_func(a_1, a_2);
      Py_XDECREF(a_1);
      Py_XDECREF(a_2);
      return out;
    }
    lookup_i++;
  }
  goto fail;

 fail:
  Py_XDECREF(a_1);
  Py_XDECREF(a_2);
  return NULL;
}

static PyObject *
crossProduct(PyObject *dummy, PyObject *args)
{
  PyObject	 *a_v1=NULL,	 *a_v2=NULL;
  PyObject	 *v1=NULL,	 *v2=NULL;
  PyObject	 *output=NULL;
  npy_double	 *v1_c=NULL,	 *v2_c=NULL,	 *output_c=NULL;
  npy_intp       *dims=NULL;

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &a_v1,
			&PyArray_Type, &a_v2))
    return NULL;

  /* Convert generic PyArrays to their appropriate types */
  v1 = PyArray_FROM_OTF(a_v1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (v1 == NULL) goto fail;
  v2 = PyArray_FROM_OTF(a_v2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (v2 == NULL) goto fail;

  if ((PyArray_NDIM(v1) != 2) ||
      (PyArray_NDIM(v2) != 2) ||
      (PyArray_DIMS(v1)[1] != GPY_DIMENSIONS) ||
      (PyArray_DIMS(v2)[1] != GPY_DIMENSIONS))
    goto badargument;

  /* Perform calculations */
  dims = (npy_intp*) PyArray_DIMS(v1);
  v1_c = (npy_double*)  PyArray_DATA(v1);
  v2_c = (npy_double*)  PyArray_DATA(v2);
  output = PyArray_EMPTY (2, dims, NPY_DOUBLE, 0);
  output_c = (npy_double*) PyArray_DATA(output);
  for (npy_intp i = 0; i < dims[0]; ++i)
    MyCArray_CrossProduct_InPlace(v1_c + GPY_DIMENSIONS*i,
				  v2_c + GPY_DIMENSIONS*i,
				  output_c + GPY_DIMENSIONS*i);

  /* DECREF's */
  Py_XDECREF(v1);
  Py_XDECREF(v2);
  return output;

 fail:
  Py_XDECREF(v1);
  Py_XDECREF(v2);
  return NULL;
  
 badargument:
  Py_XDECREF(v1);
  Py_XDECREF(v2);
  PyErr_BadArgument();
  return NULL;
}

static PyObject *
calculateDisplacements(PyObject *dummy, PyObject *args)
{
    PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL;
    PyObject *conf=NULL, *ii=NULL, *ij=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL, oshape[2];
    npy_int  *iic=NULL, *ijc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
                &PyArray_Type, &a_ij))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    oshape[0] = nidx;
    oshape[1] = GPY_DIMENSIONS;
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(2, oshape, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (npy_intp i = 0; i < nidx; ++i)
    {
        for (npy_intp j = 0; j < GPY_DIMENSIONS; ++j)
	{
	  outputc[GPY_DIMENSIONS*i + j] =
	      (npy_double) confc[GPY_DIMENSIONS*iic[i]+j]
	        - confc[GPY_DIMENSIONS*ijc[i]+j];
	  //NanTracker("calculateDisplacements", 1, outputc[GPY_DIMENSIONS*i + j]);
	}
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    PyErr_BadArgument();
    return NULL;
}

static PyObject *
calculateDistances(PyObject *dummy, PyObject *args)
{
    PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL;
    PyObject *conf=NULL, *ii=NULL, *ij=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
                &PyArray_Type, &a_ij))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
        outputc[i] = (npy_double) 0.00;
        for (int j = 0; j < GPY_DIMENSIONS; ++j)
	{
            outputc[i] += (npy_double) pow(confc[GPY_DIMENSIONS*iic[i]+j]
					   - confc[GPY_DIMENSIONS*ijc[i]+j], 2);
	    //NanTracker("calculateDistances", 1, outputc[i]);
	}

        outputc[i] = (npy_double) sqrt(outputc[i]);
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    PyErr_BadArgument();
    return NULL;
}

static PyObject *
calculateDistances2(PyObject *dummy, PyObject *args)
{
    PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL;
    PyObject *conf=NULL, *ii=NULL, *ij=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
                &PyArray_Type, &a_ij))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
        outputc[i] = (npy_double) 0.00;
        for (int j = 0; j < GPY_DIMENSIONS; ++j)
	{
            outputc[i] += (npy_double) pow(confc[GPY_DIMENSIONS*iic[i]+j]
					   - confc[GPY_DIMENSIONS*ijc[i]+j], 2);
	    //NanTracker("calculateDistances2", 1, outputc[i]);
	}
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    PyErr_BadArgument();
    return NULL;
}

static PyObject *
calculateCosines(PyObject *dummy, PyObject *args)
{
  PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL, *a_ik=NULL;
  PyObject *conf=NULL, *ii=NULL, *ij=NULL, *ik=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL, *ikc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
		&PyArray_Type, &a_ij,
			  &PyArray_Type, &a_ik	  ))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;
    ik = PyArray_FROM_OTF(a_ik, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ik == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(ik) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    ikc     = (npy_int*) PyArray_DATA(ik);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
      npy_double v1_norm = 0.0;
      npy_double v2_norm = 0.0;
      outputc[i] = (npy_double) 0.00;
      /* dot product */
      for (int j = 0; j < GPY_DIMENSIONS; ++j) {
	npy_double v1 = confc[GPY_DIMENSIONS*iic[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	npy_double v2 = confc[GPY_DIMENSIONS*ikc[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	v1_norm += v1 * v1;
	v2_norm += v2 * v2;
	outputc[i] += v1 * v2;
	//NanTracker("calculateCosines", 5, v1, v2, v1_norm, v2_norm, outputc[i]);
      }
      v1_norm = sqrt(v1_norm);
      v2_norm = sqrt(v2_norm);
      outputc[i] /= (v1_norm * v2_norm);
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    PyErr_BadArgument();
    return NULL;
}

PyObject *
calculateSines(PyObject *dummy, PyObject *args)
{
  PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL, *a_ik=NULL;
  PyObject *conf=NULL, *ii=NULL, *ij=NULL, *ik=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL, *ikc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
		&PyArray_Type, &a_ij,
			  &PyArray_Type, &a_ik	  ))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;
    ik = PyArray_FROM_OTF(a_ik, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ik == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(ik) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    ikc     = (npy_int*) PyArray_DATA(ik);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
      npy_double v1_norm = 0.0;
      npy_double v2_norm = 0.0;
      outputc[i] = (npy_double) 0.00;
      /* dot product */
      for (int j = 0; j < GPY_DIMENSIONS; ++j) {
	npy_double v1 = confc[GPY_DIMENSIONS*iic[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	npy_double v2 = confc[GPY_DIMENSIONS*ikc[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	v1_norm += v1 * v1;
	v2_norm += v2 * v2;
	outputc[i] += v1 * v2;
	//NanTracker("calculateSines", 5, v1, v2, v1_norm, v2_norm, outputc[i]);
      }
      v1_norm = sqrt(v1_norm);
      v2_norm = sqrt(v2_norm);
      outputc[i] /= (v1_norm * v2_norm);
      outputc[i] = sqrt(1 - pow(outputc[i],2));
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    PyErr_BadArgument();
    return NULL;
}

static PyObject *
calculateAngles(PyObject *dummy, PyObject *args)
{
  PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL, *a_ik=NULL;
  PyObject *conf=NULL, *ii=NULL, *ij=NULL, *ik=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL, *ikc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                &PyArray_Type, &a_conf,
                &PyArray_Type, &a_ii,
		&PyArray_Type, &a_ij,
			  &PyArray_Type, &a_ik	  ))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;
    ik = PyArray_FROM_OTF(a_ik, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ik == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(ik) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    ikc     = (npy_int*) PyArray_DATA(ik);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
      npy_double v1_norm = 0.0;
      npy_double v2_norm = 0.0;
      outputc[i] = (npy_double) 0.00;
      /* dot product */
      for (int j = 0; j < GPY_DIMENSIONS; ++j) {
	npy_double v1 = confc[GPY_DIMENSIONS*iic[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	npy_double v2 = confc[GPY_DIMENSIONS*ikc[i]+j]
	    - confc[GPY_DIMENSIONS*ijc[i]+j];
	v1_norm += v1 * v1;
	v2_norm += v2 * v2;
	outputc[i] += v1 * v2;
	//NanTracker("calculateAngles", 5, v1, v2, v1_norm, v2_norm, outputc[i]);
      }
      v1_norm = sqrt(v1_norm);
      v2_norm = sqrt(v2_norm);
      outputc[i] /= (v1_norm * v2_norm);
      outputc[i] = acos(outputc[i]);
      //NanTracker("calculateAngles::acos", 1, outputc[i]);
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    PyErr_BadArgument();
    return NULL;
}

static PyObject *
calculateDihedrals(PyObject *dummy, PyObject *args)
{
    PyObject *a_conf=NULL, *a_ii=NULL, *a_ij=NULL, *a_ik=NULL, *a_il=NULL;
    PyObject *conf=NULL, *ii=NULL, *ij=NULL, *ik=NULL, *il=NULL;
    PyObject *output=NULL;
    npy_intp *iishape=NULL;
    npy_int  *iic=NULL, *ijc=NULL, *ikc=NULL, *ilc=NULL;
    npy_double   *confc=NULL, *outputc=NULL;
    npy_int       nidx;

    /* Parse arguments */
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!",
			  &PyArray_Type, &a_conf,
			  &PyArray_Type, &a_ii,
			  &PyArray_Type, &a_ij,
			  &PyArray_Type, &a_ik,
			  &PyArray_Type, &a_il))
        return NULL;

    /* Convert generic PyArrays to their appropriate types */
    conf = PyArray_FROM_OTF(a_conf, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (conf == NULL) goto fail;
    ii = PyArray_FROM_OTF(a_ii, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ii == NULL) goto fail;
    ij = PyArray_FROM_OTF(a_ij, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ij == NULL) goto fail;
    ik = PyArray_FROM_OTF(a_ik, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (ik == NULL) goto fail;
    il = PyArray_FROM_OTF(a_il, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (il == NULL) goto fail;

    if (PyArray_NDIM(ii) != 1) goto badargument;
    if (PyArray_NDIM(ij) != 1) goto badargument;
    if (PyArray_NDIM(ik) != 1) goto badargument;
    if (PyArray_NDIM(il) != 1) goto badargument;
    if (PyArray_NDIM(conf) != 2) goto badargument;
    if (PyArray_DIMS(conf)[1] != GPY_DIMENSIONS) goto badargument;

    /* Perform calculations */
    iishape = PyArray_DIMS(ii);
    nidx    = iishape[0];
    iic     = (npy_int*) PyArray_DATA(ii);
    ijc     = (npy_int*) PyArray_DATA(ij);
    ikc     = (npy_int*) PyArray_DATA(ik);
    ilc     = (npy_int*) PyArray_DATA(il);
    confc   = (npy_double*)  PyArray_DATA(conf);
    output  = PyArray_SimpleNew(1, (npy_intp[]) {nidx}, NPY_DOUBLE);
    outputc = (npy_double*) PyArray_DATA(output);
    for (int i = 0; i < nidx; ++i)
    {
      npy_double  *ri=&(confc[GPY_DIMENSIONS*iic[i]]),
	*rj=&(confc[GPY_DIMENSIONS*ijc[i]]),
	*rk=&(confc[GPY_DIMENSIONS*ikc[i]]),
	*rl=&(confc[GPY_DIMENSIONS*ilc[i]]);
      npy_double *rij=NULL, *rkj=NULL, *rkl=NULL, *nijk=NULL, *njkl=NULL;
      npy_double dp, mod_nijk, mod_njkl, dsign;
      
      rij = MyCArray_Subtract(ri, rj);
      rkj = MyCArray_Subtract(rk, rj);
      rkl = MyCArray_Subtract(rk, rl);
      nijk = MyCArray_CrossProduct(rij, rkj);
      njkl = MyCArray_CrossProduct(rkj, rkl);
      if (!MyCArray_DotProduct(nijk, njkl, &dp) ||
	  !MyCArray_DotProduct(nijk, nijk, &mod_nijk) ||
	  !MyCArray_DotProduct(njkl, njkl, &mod_njkl) ||
	  !MyCArray_DotProduct(rij, njkl, &dsign)) {
	free(rij);
	free(rkj);
	free(rkl);
	free(nijk);
	free(njkl);
	goto badargument;
      }
//      if (!NanTracker("calculateDihedrals", 1, dp))
//	  goto fail;
      dp /= sqrt(mod_nijk * mod_njkl);
      if (dp >= 1.0) {
	      dp = 0.0;
      }
      else if (dp <= -1.0)
      {
	      dp = 3.141592653589793238;
      }
      else
      {
	  dp = acos(dp);
      }
      //if (!NanTracker("calculateDihedrals::acos", 1, dp))
      //	  goto fail;
      if (signbit(dsign))
	dp *= -1.0;

      outputc[i] = dp;
      free(rij);
      free(rkj);
      free(rkl);
      free(nijk);
      free(njkl);
    }

    /* DECREF's */
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    Py_XDECREF(il);
    return output;

fail:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    Py_XDECREF(il);
    return NULL;
badargument:
    Py_XDECREF(conf);
    Py_XDECREF(ii);
    Py_XDECREF(ij);
    Py_XDECREF(ik);
    Py_XDECREF(il);
    PyErr_BadArgument();
    return NULL;
}

static PyMethodDef mymethods[] = {
			  { "einsum", einsumWrapper,
			    METH_VARARGS | METH_KEYWORDS,
			    ""},
			  { "crossProduct", crossProduct,
			    METH_VARARGS | METH_KEYWORDS,
			    ""},
			  { "calculateDisplacements", calculateDisplacements,
			    METH_VARARGS | METH_KEYWORDS,
			    ""},
			  { "calculateDistances", calculateDistances,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  { "calculateDistances2", calculateDistances2,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  { "calculateCosines", calculateCosines,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  { "calculateSines", calculateSines,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  { "calculateAngles", calculateAngles,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  { "calculateDihedrals", calculateDihedrals,
			    METH_VARARGS | METH_KEYWORDS, 
			    ""},
			  {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef geometrypymodule = {
    PyModuleDef_HEAD_INIT,
    "geometrypy",   /* name of module */
    NULL,               /* module documentation, may be NULL */
    -1,                 /* size of per-interpreter state of the module,
                           or -1 if the module keeps state in global variables. */
    mymethods
};

PyMODINIT_FUNC
PyInit_geometrypy(void)
{
   import_array();
   PyObject *m;
   m = PyModule_Create(&geometrypymodule);
    if (m == NULL)
        return NULL;
    return m;
}

