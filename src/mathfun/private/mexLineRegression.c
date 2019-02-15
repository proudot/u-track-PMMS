/* MEX file for line regression according to generalized LS */

/* include to communicate with MatLab */
#include "mex.h"
#include "matrix.h"

/* other header files required */
#include "mexDataConverters.h"
#include "num.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   
   LHS
   ---
   plhs[0]: 2xn matrix with point residuals
   plhs[1]: 2x1 vector with line parameters
   plhs[2]: 2x1 vector with line parameter precision
   plhs[3]: sigma value 
   plhs[4]: number of iterations

   RHS
   ---
   prhs[0] : 2xn matrix with double point coordinates
   prhs[1] : 2xn matrix with double point coordinate uncertainty
   prhs[2] : 2x1 vector with approximations for the line parameters
   prhs[3] : 2x1 vector with the uncertainty of the parameter values
   prhs[4] : scalar with sigma0
   */
{
  double *xData,*sXData,*xRes;
  double *sRhoAndAlpha,*rhoAndAlpha;
  double sigma0;
  int n,m,sn,sm,nIter;

  /* unpack coordinate data */
  if(!mat2dblMtx(&xData,&n,&m,prhs[0]))
     mexErrMsgTxt("Invalid coordinate matrix entered");
  if(m!=2)
    mexErrMsgTxt("Invalid coordinate matrix entered");
  if(!mat2dblMtx(&sXData,&sn,&sm,prhs[1]))
    mexErrMsgTxt("Invalid uncertainty matrix entered"); 
  if((sm!=2)||(sn!=n))
    mexErrMsgTxt("Invalid uncertainty matrix entered");
  if(!mat2dblVec(&rhoAndAlpha,&m,prhs[2]))
    mexErrMsgTxt("Invalid parameter vector entered");
  if(!mat2dblVec(&sRhoAndAlpha,&m,prhs[3]))
    mexErrMsgTxt("Invalid uncertainty vector for parameters entered");
  if(!mat2dblScalar(&sigma0,prhs[4]))
    mexErrMsgTxt("Invalid sigma0 entered");
  
  /* allocate residual vector */
  xRes = (double*)mxMalloc(2*n*sizeof(double));

  nIter = linefit(xData,xData+n,sXData,sXData+n,n,rhoAndAlpha,sRhoAndAlpha,
		  &sigma0,xRes,xRes+n);

  /* pack data to send back */
  plhs[0] = dblMtx2mat(xRes,n,2);
  plhs[1] = dblMtx2mat(rhoAndAlpha,1,2);
  plhs[2] = dblMtx2mat(sRhoAndAlpha,1,2);
  plhs[3] = dblScalar2mat(sigma0);
  plhs[4] = intScalar2mat(nIter);
  

  mxFree(xRes);
  mxFree(rhoAndAlpha);
  mxFree(sRhoAndAlpha);
  mxFree(xData);
  mxFree(sXData);
}


