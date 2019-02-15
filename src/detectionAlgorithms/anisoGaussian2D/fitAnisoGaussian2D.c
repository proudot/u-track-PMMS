/* [prmVect prmStd covarianceMatrix residuals Jacobian] = fitAnisoGaussian2D(prmVect, initValues, mode);
 *
 * (c) Sylvain Berlemont, 2011 (last modified Jul 29, 2011)
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a fitAnisoGaussian2D.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.14" -I"..\..\mex\include" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output fitAnisoGaussian2D fitAnisoGaussian2D.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h> /* tolower() */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "mex.h"
#include "matrix.h"
#include "stats.h"

#define SIGN(x)	(x > 0 ? 1 : (x < 0 ? -1.0 : 0.0))

#define NPARAMS	7
#define REFMODE	"xyarstc"	/* r = sigma_x (along the feature),
				   s = sigma_y (aside the feature) */

typedef struct argStruct
{
  double xi, yi;
  double A;
  double g;
  double ct, st;
  double c2t, s2t;
  double sx2, sy2;
  double sx3, sy3;
  double a, b, c;
} argStruct_t;

typedef int(*pfunc_t)(gsl_matrix*, int, int, argStruct_t*);

typedef struct dataStruct
{
  int nx, ny, np;
  double *pixels;
  int *estIdx;
  int *idx;
  int nValid; /* number of non-NaN pixels */
  double *x_init;
  double prmVect[NPARAMS];
  pfunc_t *dfunc;
  gsl_vector *residuals;
  gsl_matrix *J;
} dataStruct_t;

/* 2 A gg (aa xi + bb yi) */
static int df_dx(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  double xi = argStruct->xi;
  double yi = argStruct->yi;
  double A = argStruct->A;
  double g = argStruct->g;
  double a = argStruct->a;
  double b = argStruct->b;
  
  gsl_matrix_set(J, i, k, 2 * A * g * (a * xi + b * yi));

  return 0;
}

/* 2 A g (bb xi + cc yi) */
static int df_dy(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  double xi = argStruct->xi;
  double yi = argStruct->yi;
  double A = argStruct->A;
  double g = argStruct->g;
  double b = argStruct->b;
  double c = argStruct->c;
	
  gsl_matrix_set(J, i, k, 2 * A * g * (b * xi + c * yi));
	
  return 0;
}

/* g */
static int df_dA(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  gsl_matrix_set(J, i, k, argStruct->g);
  return 0;
}

/* (A g (xi ct + yi st)^2)/sx^3 */
static int df_dsx(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  double xi = argStruct->xi;
  double yi = argStruct->yi;
  double A = argStruct->A;
  double g = argStruct->g;
  double ct = argStruct->ct;
  double st = argStruct->st;
  double sx3 = argStruct->sx3;
	
  double r = xi * ct + yi * st;
	
  gsl_matrix_set(J, i, k, A * g * r * r / sx3);
	
  return 0;
}

/* (A g (yi ct - xi st)^2)/sy^3 */
static int df_dsy(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  double xi = argStruct->xi;
  double yi = argStruct->yi;
  double A = argStruct->A;
  double g = argStruct->g;
  double ct = argStruct->ct;
  double st = argStruct->st;
  double sy3 = argStruct->sy3;
	
  double r = yi * ct - xi * st;
	
  gsl_matrix_set(J, i, k, A * g * r * r / sy3);
	
  return 0;
}

/* -((A g (sx^2 - sy^2) (-2 xi yi c2t + (xi^2 - yi^2) s2t)) / (2 sx^2 sy^2)) */
static int df_dt(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  double xi = argStruct->xi;
  double yi = argStruct->yi;
  double A = argStruct->A;
  double g = argStruct->g;
  double c2t = argStruct->c2t;
  double s2t = argStruct->s2t;
  double sx2 = argStruct->sx2;
  double sy2 = argStruct->sy2;
	
  gsl_matrix_set(J, i, k, -((A * g * (sx2 - sy2) *
			     (-2 * xi * yi * c2t +
			      (xi * xi - yi * yi) * s2t)) / (2 * sx2 * sy2)));
	
  return 0;
}

/* 1 */
static int df_dC(gsl_matrix *J, int i, int k, argStruct_t *argStruct)
{
  gsl_matrix_set(J, i, k, 1);
  return 0;
}

static int f(const gsl_vector *x, void *params, gsl_vector *f)
{
  dataStruct_t *dataStruct = (dataStruct_t *)params;
  int i,idx;

  int nx = dataStruct->nx;
  int ny = dataStruct->ny;
  int nx_div2 = (nx-1) >> 1;
  int ny_div2 = (ny-1) >> 1;
    
  double *pixels = dataStruct->pixels;
  
  /* update prmVect with new estimates */
  for (i=0; i<dataStruct->np; ++i) {
    dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
  }

  double xp = dataStruct->prmVect[0];
  double yp = dataStruct->prmVect[1];
  double A = dataStruct->prmVect[2];
  double sx = fabs(dataStruct->prmVect[3]);
  double sy = fabs(dataStruct->prmVect[4]);
  double t = dataStruct->prmVect[5];
  double C = dataStruct->prmVect[6];

  double ct = cos(t);
  double st = sin(t);
  double s2t = sin(2 * t);
  double sx2 = sx * sx;
  double sy2 = sy * sy;
	
  double xi, yi;
  double a = ct * ct / (2 * sx2) + st * st / (2 * sy2);
  double b = s2t / (4 * sx2) - s2t / (4 * sy2);
  double c = st * st / (2 * sx2) + ct * ct / (2 * sy2);
    
  div_t divRes;
	
  for (i=0; i < dataStruct->nValid; ++i)
    {
      idx = dataStruct->idx[i];
      divRes = div(idx, ny);
	  
      xi = divRes.quot-nx_div2-xp;
      yi = divRes.rem-ny_div2-yp;
	  
      gsl_vector_set(f, i, A * exp(-a * xi * xi - yi * (2 * b * xi + c * yi)) + C - pixels[idx]);
    }
	
  return GSL_SUCCESS;
}

static int df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  dataStruct_t *dataStruct = (dataStruct_t *)params;
  int i,idx,k;

  int nx = dataStruct->nx;
  int ny = dataStruct->ny;
  int nx_div2 = (nx-1) >> 1;
  int ny_div2 = (ny-1) >> 1;
    
  /* update prmVect with new estimates */
  for (i=0; i<dataStruct->np; ++i) {
    dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
  }
    
  double xp = dataStruct->prmVect[0];
  double yp = dataStruct->prmVect[1];
  double A = dataStruct->prmVect[2];
  double sx = fabs(dataStruct->prmVect[3]);
  double sy = fabs(dataStruct->prmVect[4]);
  double t = dataStruct->prmVect[5];
	
  double xi, yi;
  double ct = cos(t);
  double st = sin(t);
  double c2t = cos(2 * t);
  double s2t = sin(2 * t);
  double sx2 = sx * sx;
  double sy2 = sy * sy;
  double a = ct * ct / (2 * sx2) + st * st / (2 * sy2);
  double b = s2t / (4 * sx2) - s2t / (4 * sy2);
  double c = st * st / (2 * sx2) + ct * ct / (2 * sy2);

  argStruct_t argStruct;
  argStruct.A = A;
  argStruct.ct = ct;
  argStruct.st = st;
  argStruct.c2t = c2t;
  argStruct.s2t = s2t;
  argStruct.sx2 = sx2;
  argStruct.sy2 = sy2;
  argStruct.sx3 = sx2 * sx;
  argStruct.sy3 = sy2 * sy;
  argStruct.a = a;
  argStruct.b = b;
  argStruct.c = c;
    
  div_t divRes;

  for (i=0; i<dataStruct->nValid; ++i)
    {
      idx = dataStruct->idx[i];
      divRes = div(idx, ny);
		
      xi = divRes.quot-nx_div2-xp;
      yi = divRes.rem-ny_div2-yp;
		
      argStruct.xi = xi;
      argStruct.yi = yi;
      argStruct.g = exp(-a * xi * xi - yi * (2 * b * xi + c * yi));
		
      for (k=0; k<dataStruct->np; ++k)
	dataStruct->dfunc[k](J, i, k, &argStruct);
    }
  return GSL_SUCCESS;
}

static int fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
  dataStruct_t *dataStruct = (dataStruct_t *)params;
  int i, idx, k;

  int nx = dataStruct->nx;
  int ny = dataStruct->ny;
  int nx_div2 = (nx-1) >> 1;
  int ny_div2 = (ny-1) >> 1;
    
  double *pixels = dataStruct->pixels;
    
  /* update prmVect with new estimates */
  for (i=0; i<dataStruct->np; ++i)
    dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    
  double xp = dataStruct->prmVect[0];
  double yp = dataStruct->prmVect[1];
  double A = dataStruct->prmVect[2];
  double sx = fabs(dataStruct->prmVect[3]);
  double sy = fabs(dataStruct->prmVect[4]);
  double t = dataStruct->prmVect[5];
  double C = dataStruct->prmVect[6];
	
  double xi, yi;
  double ct = cos(t);
  double st = sin(t);
  double c2t = cos(2 * t);
  double s2t = sin(2 * t);
  double sx2 = sx * sx;
  double sy2 = sy * sy;
  double a = ct * ct / (2 * sx2) + st * st / (2 * sy2);
  double b = s2t / (4 * sx2) - s2t / (4 * sy2);
  double c = st * st / (2 * sx2) + ct * ct / (2 * sy2);
	
  argStruct_t argStruct;
  argStruct.A = A;
  argStruct.ct = ct;
  argStruct.st = st;
  argStruct.c2t = c2t;
  argStruct.s2t = s2t;
  argStruct.sx2 = sx2;
  argStruct.sy2 = sy2;
  argStruct.sx3 = sx2 * sx;
  argStruct.sy3 = sy2 * sy;
  argStruct.a = a;
  argStruct.b = b;
  argStruct.c = c;
    
  div_t divRes;

  for (i=0; i<dataStruct->nValid; ++i)
    {
      idx = dataStruct->idx[i];
      divRes = div(idx, ny);
		
      xi = divRes.quot-nx_div2-xp;
      yi = divRes.rem-ny_div2-yp;

      argStruct.xi = xi;
      argStruct.yi = yi;
      argStruct.g = exp(-a * xi * xi - yi * (2 * b * xi + c * yi));

      gsl_vector_set(f, i, A*argStruct.g + C - pixels[idx]);
        
      for (k=0; k<dataStruct->np; ++k)
	dataStruct->dfunc[k](J, i, k, &argStruct);
    }
  return GSL_SUCCESS;
}

static int MLalgo(struct dataStruct *data)
{
  int i;

  /* declare solvers */
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
    
  /* number of parameters to optimize */
  const size_t p = data->np;
    
  gsl_vector_view x = gsl_vector_view_array(data->x_init, p);
    
  const gsl_rng_type *type;
    
  gsl_rng_env_setup();
    
  type = gsl_rng_default;
    
  gsl_multifit_function_fdf fit;
  fit.f = &f;
  fit.df = &df;
  fit.fdf = &fdf;
  fit.n = data->nValid;
  fit.p = p;
  fit.params = data;
  
  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, data->nValid, data->np);
  gsl_multifit_fdfsolver_set(s, &fit, &x.vector);
    
  int status;
  int iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    if (status)
      break;

    status = gsl_multifit_test_delta(s->dx, s->x, 1e-8, 1e-8);
  }
  while (status == GSL_CONTINUE && iter < 500);
    
  for (i=0; i<data->np; ++i)
    data->prmVect[data->estIdx[i]] = gsl_vector_get(s->x, i);
    
  /* copy residual */
  data->residuals = gsl_vector_alloc(data->nValid);
  gsl_vector_memcpy(data->residuals, s->f);

  /* copy Jacobian */
  data->J = gsl_matrix_alloc(data->nValid, data->np);

  gsl_matrix_memcpy(data->J, s->J);

  gsl_multifit_fdfsolver_free(s);
	
  return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, k = 0, l = 0;
		
  /* check inputs */
  if (nrhs < 3) mexErrMsgTxt("Inputs should be: data, prmVect, mode.");
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Data input must be double array.");
  if (mxGetNumberOfElements(prhs[1])!=NPARAMS || !mxIsDouble(prhs[1])) mexErrMsgTxt("Incorrect parameter vector format.");
  if (!mxIsChar(prhs[2])) mexErrMsgTxt("Mode needs to be a string.");

  size_t nx = mxGetN(prhs[0]);
  size_t ny = mxGetM(prhs[0]);
  int N = nx*ny;

  /* read mode input */
  int np = mxGetNumberOfElements(prhs[2]);
  char *mode;
  mode = (char*)malloc(sizeof(char)*(np+1));
  mxGetString(prhs[2], mode, np+1);

  for (i = 0; i < strlen(mode); ++i)
    mode[i] = tolower(mode[i]);

  np = 0; /* number of parameters to fit */
  for (i = 0; i < NPARAMS; ++i)
    if (strchr(mode, REFMODE[i])!=NULL) { np++; }

  if (!np)
    mexErrMsgTxt("Unknown mode.");

  /* allocate */
  dataStruct_t data;
  data.nx = nx;
  data.ny = ny;
  data.np = np;
  data.pixels = mxGetPr(prhs[0]);
  data.estIdx = (int*)malloc(sizeof(int)*np);
  memcpy(data.prmVect, mxGetPr(prhs[1]), NPARAMS*sizeof(double));
  data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);

  /* read mask/pixels */
  data.nValid = N;
  for (i=0; i<N; ++i)
    data.nValid -= (int)mxIsNaN(data.pixels[i]);

  data.idx = (int*)malloc(sizeof(int)*data.nValid);
  int *nanIdx = (int*)malloc(sizeof(int)*(N-data.nValid));

  for (i=0; i<N; ++i)
    if (!mxIsNaN(data.pixels[i]))
      data.idx[k++] = i;
    else
      nanIdx[l++] = i;

  np = 0;
  if (strchr(mode, 'x')!=NULL) {data.estIdx[np] = 0; data.dfunc[np++] = df_dx;}
  if (strchr(mode, 'y')!=NULL) {data.estIdx[np] = 1; data.dfunc[np++] = df_dy;}
  if (strchr(mode, 'a')!=NULL) {data.estIdx[np] = 2; data.dfunc[np++] = df_dA;}
  if (strchr(mode, 'r')!=NULL) {data.estIdx[np] = 3; data.dfunc[np++] = df_dsx;}
  if (strchr(mode, 's')!=NULL) {data.estIdx[np] = 4; data.dfunc[np++] = df_dsy;}
  if (strchr(mode, 't')!=NULL) {data.estIdx[np] = 5; data.dfunc[np++] = df_dt;}
  if (strchr(mode, 'c')!=NULL) {data.estIdx[np] = 6; data.dfunc[np++] = df_dC;}
    
  data.x_init = (double*)malloc(sizeof(double)*np);
  for (i=0; i<np; ++i)
    data.x_init[i] = data.prmVect[data.estIdx[i]];
    
  MLalgo(&data);
    
  /* parameters */
  if (nlhs > 0)
    {
      /* Make sure the angle parameter lies in -pi/2...pi/2 */
      float tmp = data.prmVect[5];
      if (tmp > M_PI_2 || tmp < -M_PI_2)
	{
	  tmp = fmodl(tmp + M_PI_2,M_PI);
	  tmp -= SIGN(tmp) * M_PI_2;
	  data.prmVect[5] = tmp;
	}

      /* Make sure sigma_x and sigma_y are positive */
      data.prmVect[3] = fabs(data.prmVect[3]);
      data.prmVect[4] = fabs(data.prmVect[4]);

      plhs[0] = mxCreateDoubleMatrix(1, NPARAMS, mxREAL);
      memcpy(mxGetPr(plhs[0]), data.prmVect, NPARAMS * sizeof(double));
    }
    
  double RSS = 0.0;
  double* resValid = NULL;
  /* standard dev. of parameters & covariance matrix */

  if (nlhs > 1) {
      resValid = (double*)malloc(data.nValid*sizeof(double));
      for (i=0; i<data.nValid; ++i) {
          resValid[i] = gsl_vector_get(data.residuals, i);
          RSS += resValid[i]*resValid[i];
      }
      gsl_matrix *covar = gsl_matrix_alloc(np, np);
      gsl_multifit_covar(data.J, 0.0, covar);
      double iRSS = RSS/(data.nValid - data.np - 1);
      plhs[1] = mxCreateDoubleMatrix(1, data.np, mxREAL);
      double *prmStd = mxGetPr(plhs[1]);
      for (i=0; i<data.np; ++i)
          prmStd[i] = sqrt(iRSS*gsl_matrix_get(covar, i, i));
      if (nlhs > 2) {
          plhs[2] = mxCreateDoubleMatrix(np, np, mxREAL);
          /* cov. matrix is symmetric, no need to transpose */
          memcpy(mxGetPr(plhs[2]), covar->data, np*np*sizeof(double));
      }
      gsl_matrix_free(covar);
  }
  
  /* residuals */
  if(nlhs > 3){
      const char *fieldnames[]={"data", "pval", "mean", "std", "RSS"};
      mwSize dims[2]={1, 1};
      plhs[3]=mxCreateStructArray(2, dims, 5, fieldnames);
      mxArray *val=mxCreateDoubleMatrix(nx, ny, mxREAL);
      double *res=mxGetPr(val);
      
      double mean=0.0, std=0.0;
      for(i=0;i<data.nValid;++i){
          res[data.idx[i]]=resValid[i];
          mean+=resValid[i];
      }
      std=sqrt((RSS-mean*mean/data.nValid)/(data.nValid-1));
      mean/=data.nValid;
      
      for(i=0;i<N-data.nValid;++i) {
          res[nanIdx[i]]=mxGetNaN();
      }
      
      double pval=ksone(resValid, data.nValid, mean, std);
      mxSetFieldByNumber(plhs[3], 0, 0, val);
      mxSetFieldByNumber(plhs[3], 0, 1, mxCreateDoubleScalar(pval));
      mxSetFieldByNumber(plhs[3], 0, 2, mxCreateDoubleScalar(mean));
      mxSetFieldByNumber(plhs[3], 0, 3, mxCreateDoubleScalar(std));
      mxSetFieldByNumber(plhs[3], 0, 4, mxCreateDoubleScalar(RSS));
  }
    
  /* Jacobian */
  if (nlhs > 4)
    {
      /* convert row-major double* data.J to column-major double */
      plhs[4] = mxCreateDoubleMatrix(N, np, mxREAL);
      double *J = mxGetPr(plhs[4]);
      int k;
      for (k=0; k<np; ++k)
	{
	  for (i=0; i<data.nValid; ++i)
	    J[data.idx[i]+k*N] = gsl_matrix_get(data.J, i, k);
	  for (i=0; i<N-data.nValid; ++i)
	    J[nanIdx[i]+k*N] = mxGetNaN();
	}
    }
    
  gsl_vector_free(data.residuals);
  gsl_matrix_free(data.J);
  free(data.x_init);
  free(nanIdx);
  free(data.idx);
  free(data.dfunc);
  free(data.estIdx);
  free(mode);
  free(resValid);
}

/* Compile line for MAX OSX: */

/* export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2010b.app/bin/maci64
   && gcc -ansi -Wall -g -DARRAY_ACCESS_INLINING
   -I. -I/Applications/MATLAB_R2010b.app/extern/include
   -L/Applications/MATLAB_R2010b.app/bin/maci64 -lmx -lmex -lgsl
   -lgslcblas -lmat fitAnisoGaussian2D.c
*/

/* Compile line for Linux: */

/* export LD_RUN_PATH=/usr/local/Matlab/bin/glnxa64 && gcc -ansi -Wall
   -g -DARRAY_ACCESS_INLINING -I. -L/usr/local/Matlab/bin/glnxa64
   -I/usr/local/Matlab/extern/include fitAnisoGaussian2D.c -lmx -lmex
   -lgsl -lgslcblas
*/

int main(int argc, char** argv)
{ 
  int i, k = 0, l = 0;
  int nx = 15;
  int N = nx*nx;
    
  dataStruct_t data;
  data.nx = nx;
  data.np = NPARAMS;
  data.pixels = (double*)malloc(sizeof(double)*N);
	
  /* fill with noise */
  for (i=0; i<N; ++i)
    data.pixels[i] = rand();
    
  data.estIdx = (int*)malloc(sizeof(int) * NPARAMS);
  data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * NPARAMS);
	
  /* read mask/pixels */
  data.nValid = N;
  for (i=0; i<N; ++i) {
    data.nValid -= (int)mxIsNaN(data.pixels[i]);
  }
  data.idx = (int*)malloc(sizeof(int)*data.nValid);
  int *nanIdx = (int*)malloc(sizeof(int)*(N-data.nValid));

  for (i=0; i<N; ++i)
    if (!mxIsNaN(data.pixels[i]))
      data.idx[k++] = i;
    else
      nanIdx[l++] = i;
    
  data.prmVect[0] = 0;
  data.prmVect[1] = 0;
  data.prmVect[2] = 5;
  data.prmVect[3] = 1;
  data.prmVect[4] = 1;
  data.prmVect[5] = 0;
  data.prmVect[6] = 0;
	
  data.estIdx[0] = 0; data.dfunc[0] = df_dx;
  data.estIdx[1] = 1; data.dfunc[1] = df_dy;
  data.estIdx[2] = 2; data.dfunc[2] = df_dA;
  data.estIdx[3] = 3; data.dfunc[3] = df_dsx;
  data.estIdx[4] = 4; data.dfunc[4] = df_dsy;
  data.estIdx[5] = 5; data.dfunc[5] = df_dt;
  data.estIdx[6] = 6; data.dfunc[6] = df_dC;
    
  data.x_init = (double*)malloc(sizeof(double)*NPARAMS);
  for (i=0; i<NPARAMS; ++i)
    data.x_init[i] = data.prmVect[data.estIdx[i]];
    
  MLalgo(&data);
    
  gsl_vector_free(data.residuals);
  gsl_matrix_free(data.J);
  free(data.x_init);
  free(nanIdx);
  free(data.idx);
  free(data.dfunc);
  free(data.estIdx);
  free(data.pixels);
    
  return 0;
}

