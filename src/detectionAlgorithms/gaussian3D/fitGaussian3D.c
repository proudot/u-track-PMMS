/* [prmVect, prmStd, covarianceMatrix, residuals, Jacobian] = fitGaussiandD(prmVect, initValues, mode);
 *
 * Copyright (c) 2013 Francois Aguet
 *
 * Compilation:
 * mex -I/usr/local/include -I../../common/mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a fitGaussian3D.c
 * mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\extern\mex\include\gsl-1.15" -I"..\..\common\mex\include" "..\..\extern\mex\lib\gsl.lib" "..\..\extern\mex\lib\cblas.lib" -output fitGaussian3D fitGaussian3D.c
 * 
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a fitGaussian3D.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.15" -I"..\..\mex\include" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output fitGaussian3D fitGaussian3D.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "mex.h"
#include "matrix.h"
#include "stats.h"

#define NPARAMS 7
#define refMode "xyzasrc" // s = x,y sigma; r = z sigma = rho


typedef struct aStruct {
    double xi, yi, zi, A, g, sigma2, sigma3, rho2, rho3;
} argStruct_t;

typedef int(*pfunc_t)(gsl_matrix*, int, int, argStruct_t*);

typedef struct dataStruct {
    int nx, ny, nz, np;
    double *pixels;
    double *gx, *gy, *gz;
    int *estIdx;
    int *idx;
    int nValid; /* number of non-NaN pixels */
    double *x_init;
    double prmVect[NPARAMS];
    pfunc_t *dfunc;
    gsl_vector *residuals;
    gsl_matrix *J;
    double maxIter, eAbs, eRel;
} dataStruct_t;



// argStruct.xi contains (x-xp)
static int df_dx(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*xi*g);
    return 0;
}

static int df_dy(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*yi*g);
    return 0;
}

static int df_dz(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double zi = argStruct->zi;
    double g = argStruct->g;
    double A = argStruct->A;
    double r2 = argStruct->rho2;
    gsl_matrix_set(J, i, k, A/r2*zi*g);
    return 0;
}

static int df_dA(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, argStruct->g);
    return 0;
}

static int df_ds(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s3 = argStruct->sigma3;
    gsl_matrix_set(J, i, k, (xi*xi + yi*yi)*A/s3*g);
    return 0;
}

static int df_dr(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double zi = argStruct->zi;
    double g = argStruct->g;
    double A = argStruct->A;
    double r3 = argStruct->rho3;
    gsl_matrix_set(J, i, k, zi*zi*A/r3*g);
    return 0;
}

static int df_dc(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, 1);
    return 0;
}



static int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int ny = dataStruct->ny;
    int nz = dataStruct->nz;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    double *gz = dataStruct->gz;
    
    // update prmVect with new estimates
    int i;
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double zp = dataStruct->prmVect[2];
    double A = dataStruct->prmVect[3];
    double sigma = fabs(dataStruct->prmVect[4]);
    double rho = fabs(dataStruct->prmVect[5]);
    double c = dataStruct->prmVect[6];
    
    // calculate components of Gaussian (separable)
    double s = 2.0*sigma*sigma;
    double r = 2.0*rho*rho;
    double k;
    for (i=0; i<nx; ++i) {
        k = i-xp;
        gx[i] = exp(-k*k/s);
    }
    for (i=0; i<ny; ++i) {
        k = i-yp;
        gy[i] = exp(-k*k/s);
    }
    for (i=0; i<nz; ++i) {
        k = i-zp;
        gz[i] = exp(-k*k/r);
    }
    
    div_t divRes, divResZ;
    int idx;
    int nxy = nx*ny;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divResZ = div(idx, nxy);
        // Matlab indexing is column-major
        divRes = div(divResZ.rem, ny);
        gsl_vector_set(f, i, A*gx[divRes.quot]*gy[divRes.rem]*gz[divResZ.quot]+c - pixels[idx]);
    }
    return GSL_SUCCESS;
}



static int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int ny = dataStruct->ny;
    int nz = dataStruct->nz;
    
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    double *gz = dataStruct->gz;
    
    // update prmVect with new estimates
    int i;
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double zp = dataStruct->prmVect[2];
    double A = dataStruct->prmVect[3];
    double sigma = fabs(dataStruct->prmVect[4]);
    double rho = fabs(dataStruct->prmVect[5]);
    
    double sigma2 = sigma*sigma;
    double rho2 = rho*rho;
    
    double s = 2.0*sigma2;
    double r = 2.0*rho2;
    double k;
    for (i=0; i<nx; ++i) {
        k = i-xp;
        gx[i] = exp(-k*k/s);
    }
    for (i=0; i<ny; ++i) {
        k = i-yp;
        gy[i] = exp(-k*k/s);
    }
    for (i=0; i<nz; ++i) {
        k = i-zp;
        gz[i] = exp(-k*k/r);
    }
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.rho2 = rho2;
    argStruct.rho3 = rho2*rho;
    argStruct.A = A;
    
    div_t divRes, divResZ;
    int idx;
    int nxy = nx*ny;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divResZ = div(idx, nxy);
        divRes = div(divResZ.rem, ny);
        argStruct.xi = divRes.quot - xp;
        argStruct.yi = divRes.rem - yp;
        argStruct.zi = divResZ.quot - zp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem]*gz[divResZ.quot];
        
        for (int k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



static int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int ny = dataStruct->ny;
    int nz = dataStruct->nz;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    double *gz = dataStruct->gz;
    
    // update prmVect with new estimates
    int i;
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double zp = dataStruct->prmVect[2];
    double A = dataStruct->prmVect[3];
    double sigma = fabs(dataStruct->prmVect[4]);
    double rho = fabs(dataStruct->prmVect[5]);
    double c = dataStruct->prmVect[6];
    
    double sigma2 = sigma*sigma;
    double rho2 = rho*rho;
    double s = 2.0*sigma2;
    double r = 2.0*rho2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.rho2 = rho2;
    argStruct.rho3 = rho2*rho;
    argStruct.A = A;
    
    double k;
    for (i=0; i<nx; ++i) {
        k = i-xp;
        gx[i] = exp(-k*k/s);
    }
    for (i=0; i<ny; ++i) {
        k = i-yp;
        gy[i] = exp(-k*k/s);
    }
    for (i=0; i<nz; ++i) {
        k = i-zp;
        gz[i] = exp(-k*k/r);
    }
    
    div_t divRes, divResZ;
    int idx;
    int nxy = nx*ny;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divResZ = div(idx, nxy);
        divRes = div(divResZ.rem, ny);
        
        argStruct.xi = divRes.quot - xp;
        argStruct.yi = divRes.rem - yp;
        argStruct.zi = divResZ.quot - zp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem]*gz[divResZ.quot];
        gsl_vector_set(f, i, A*argStruct.g+c - pixels[idx]);
        
        for (int k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



static int MLalgo(struct dataStruct *data) {
    
    // declare solvers
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    
    gsl_vector_view x = gsl_vector_view_array(data->x_init, data->np);
    
    const gsl_rng_type *type;
    
    gsl_rng_env_setup();
    
    type = gsl_rng_default;
    
    gsl_multifit_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    f.n = data->nValid;
    f.p = data->np;
    f.params = data;
    
    
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, data->nValid, data->np);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);
    
    int status, status2;
    int iter = 0;
    gsl_vector *gradt = gsl_vector_alloc(data->np);
    
    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        
        status = gsl_multifit_test_delta(s->dx, s->x, data->eAbs, data->eRel);
        gsl_multifit_gradient(s->J, s->f, gradt);
        status2 = gsl_multifit_test_gradient(gradt, data->eAbs);
    }
    while ((status == GSL_CONTINUE || status2 == GSL_CONTINUE) && iter < data->maxIter);
    
    gsl_vector_free(gradt);
    
    int i;
    for (i=0; i<data->np; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(s->x, i);
    }
    // force sigma, rho to positive values
    data->prmVect[4] = fabs(data->prmVect[4]);
    data->prmVect[5] = fabs(data->prmVect[5]);
    
    // copy residuals
    data->residuals = gsl_vector_alloc(data->nValid);
    gsl_vector_memcpy(data->residuals, s->f);
    
    /* copy Jacobian */
    data->J = gsl_matrix_alloc(data->nValid, data->np);
    gsl_matrix_memcpy(data->J, s->J);
    
    gsl_multifit_fdfsolver_free(s);
    return 0;
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    dataStruct_t data;
    
    // check inputs
    if (nrhs < 3) mexErrMsgTxt("Inputs should be: data, prmVect, mode.");
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Data input must be double array.");
    int ndims = mxGetNumberOfDimensions(prhs[0]);
    if (ndims!=3) mexErrMsgTxt("Input data must be 3D.");
    const int* dims = mxGetDimensions(prhs[0]);
    int ny = dims[0];
    int nx = dims[1];
    int nz = dims[2];
    
    // verify initial values
    int nparam = mxGetNumberOfElements(prhs[1]);
    if (!(nparam==6||nparam==7) || !mxIsDouble(prhs[1])) mexErrMsgTxt("Incorrect parameter vector format.");
    
    // verify mode
    if (!mxIsChar(prhs[2])) mexErrMsgTxt("Mode needs to be a string.");
    
    if (nrhs<4) {
        data.maxIter = 500;
        data.eAbs = 1e-8;
        data.eRel = 1e-8;
    } else {
        if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=3) mexErrMsgTxt("Options must must be double array with 3 elements.");
        double *options = mxGetPr(prhs[3]);
        data.maxIter = options[0];
        data.eAbs = options[1];
        data.eRel = options[2];
    }
    
    // read mode input
    int np = (int)mxGetNumberOfElements(prhs[2]);
    char *mode;
    mode = (char*)malloc(sizeof(char)*(np+1));
    mxGetString(prhs[2], mode, np+1);
    
    int i;
    for (i=0; i<strlen(mode); ++i) {
        mode[i] = tolower(mode[i]);
    }
    
    np = 0; // number of parameters to fit
    for (i=0; i<NPARAMS; ++i) {
        if (strchr(mode, refMode[i])!=NULL) { np++; }
    }
    if (np==0) mexErrMsgTxt("Unknown mode.");
    
    // allocate
    data.nx = nx;
    data.ny = ny;
    data.nz = nz;
    data.np = np; // # params
    data.pixels = mxGetPr(prhs[0]);
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*ny);
    data.gz = (double*)malloc(sizeof(double)*nz);
    data.estIdx = (int*)malloc(sizeof(int)*np);
    double *prms = mxGetPr(prhs[1]);
    if (nparam==NPARAMS) {
        memcpy(data.prmVect, prms, NPARAMS*sizeof(double));
    } else {
        for (i=0;i<6;++i) {
            data.prmVect[i] = prms[i];
        }
        data.prmVect[5] = prms[4];
        data.prmVect[6] = prms[5];
    }
    data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);
    
    // read mask/pixels
    int N = nx*ny*nz;
    data.nValid = N;
    for (i=0; i<N; ++i) {
        data.nValid -= (int)mxIsNaN(data.pixels[i]);
    }
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int *nanIdx = (int*)malloc(sizeof(int)*(N-data.nValid));
    int k = 0, l = 0;
    for (i=0; i<N; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        } else {
            nanIdx[l++] = i;
        }
    }
    
    np = 0;
    if (strchr(mode, 'x')!=NULL) {data.estIdx[np] = 0; data.dfunc[np++] = df_dx;}
    if (strchr(mode, 'y')!=NULL) {data.estIdx[np] = 1; data.dfunc[np++] = df_dy;}
    if (strchr(mode, 'z')!=NULL) {data.estIdx[np] = 2; data.dfunc[np++] = df_dz;}
    if (strchr(mode, 'a')!=NULL) {data.estIdx[np] = 3; data.dfunc[np++] = df_dA;}
    if (strchr(mode, 's')!=NULL) {data.estIdx[np] = 4; data.dfunc[np++] = df_ds;}
    if (strchr(mode, 'r')!=NULL) {data.estIdx[np] = 5; data.dfunc[np++] = df_dr;}
    if (strchr(mode, 'c')!=NULL) {data.estIdx[np] = 6; data.dfunc[np++] = df_dc;}
    
    data.x_init = (double*)malloc(sizeof(double)*np);
    for (i=0; i<np; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    
    MLalgo(&data);
    
    // return fitted parameter values
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(1, NPARAMS, mxREAL);
        memcpy(mxGetPr(plhs[0]), data.prmVect, NPARAMS*sizeof(double));
    }
    
    // standard dev. of parameters & covariance matrix
    double RSS = 0.0;
    double* resValid = NULL;
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
        for (i=0; i<data.np; ++i) {
            prmStd[i] = sqrt(iRSS*gsl_matrix_get(covar, i, i));
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(np, np, mxREAL);
            // cov. matrix is symmetric, no need to transpose
            memcpy(mxGetPr(plhs[2]), covar->data, np*np*sizeof(double));
        }
        gsl_matrix_free(covar);
    }
    
    // residuals
    if (nlhs > 3) {
        const char *fieldnames[] = {"data", "hAD", "mean", "std", "RSS"};
        mwSize dims[2] = {1, 1};
        plhs[3] = mxCreateStructArray(2, dims, 5, fieldnames);
        
        const mwSize dims2[3] = {ny,nx,nz};
        mxArray *val = mxCreateNumericArray(3, dims2, mxDOUBLE_CLASS, mxREAL);
        double* res = mxGetPr(val);
        
        double mean = 0.0, std = 0.0;
        for (i=0; i<data.nValid; ++i) {
            res[data.idx[i]] = resValid[i];
            mean += resValid[i];
        }
        std = sqrt((RSS-mean*mean/data.nValid)/(data.nValid-1));
        mean /= data.nValid;
        
        for (i=0; i<N-data.nValid; ++i) {
            res[nanIdx[i]] = mxGetNaN();
        }
        
        // A-D test, case 2: mean known
        unsigned char hAD = adtest(resValid, data.nValid, 2, 0.0, std, 0.05);
        mxSetFieldByNumber(plhs[3], 0, 0, val);
        mxSetFieldByNumber(plhs[3], 0, 1, mxCreateLogicalScalar(hAD));
        mxSetFieldByNumber(plhs[3], 0, 2, mxCreateDoubleScalar(mean));
        mxSetFieldByNumber(plhs[3], 0, 3, mxCreateDoubleScalar(std));
        mxSetFieldByNumber(plhs[3], 0, 4, mxCreateDoubleScalar(RSS));
    }
    
    // Jacobian
    if (nlhs > 4) {
        // convert row-major double* data.J->data to column-major double*
        plhs[4] = mxCreateDoubleMatrix(N, np, mxREAL);
        double *J = mxGetPr(plhs[4]);
        for (int k=0; k<np; ++k) {
            for (i=0; i<data.nValid; ++i) {
                J[data.idx[i]+k*N] = gsl_matrix_get(data.J, i, k);
            }
            for (i=0; i<N-data.nValid; ++i) {
                J[nanIdx[i]+k*N] = mxGetNaN();
            }
        }
    }
    
    free(resValid);
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    free(data.x_init);
    free(nanIdx);
    free(data.idx);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    free(data.gz);
    free(mode);
}

