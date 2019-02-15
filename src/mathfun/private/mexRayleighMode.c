/* MEX file for getting the Rayliegh mode of a map */

/* include to communicate with MatLab */
#include "mex.h"
#include "matrix.h"

/* other header files required */
#include <math.h>
#include "mexUtils.h"
#include "tools.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   
   LHS
   ---
   plhs[0]: scalar with noise estimate
   plhs[1]: scalar with threshold
   plhs[2]: nInGradBufx1 matrix with double valued gradient values

   RHS
   ---
   prhs[0] : nx x ny matrix with a Rayleigh distributed map
   prhs[1] : scalar with risk probability

   */
{
  int status;
  double *buf,*gradMagn = NULL;
  double c,mode,thresh;
  int     nx,ny,nVals;

  status = mat2dblMtx(&buf,&nx,&ny,prhs[0]);
  if(!status)
    return;

  status = mat2dblScalar(&c,prhs[1]);

  tbx_dGetRayleighMode(buf,nx*ny,&gradMagn,&nVals,&mode);
  thresh = sqrt(-2.0*log(c))*mode;
  
  plhs[2] = dblMtx2mat(gradMagn,nVals,1);

  if(gradMagn) 
    free(gradMagn);

  plhs[0]=dblScalar2mat(mode);
  plhs[1]=dblScalar2mat(thresh);

}
