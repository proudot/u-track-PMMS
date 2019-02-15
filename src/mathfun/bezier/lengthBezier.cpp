#include <mex.h>

#include <Wm5BezierCurve2.h>
#include <Wm5BezierCurve3.h>
#include <Wm5Vector2.h>
#include <Wm5Vector3.h>

using namespace Wm5;

// Get Wild Magic: from http://www.geometrictools.com/Downloads/Downloads.html

// Windows: Compile LibCore and LibMathematics in 64 bit RELEASE mode with Visual Studio!
// Windows: mex -I. -I..\..\..\extern\mex\include\wildmagic-5.7\include -L..\..\..\extern\mex\lib -lWm5Core90 -lWm5Mathematics90 lengthBezier.cpp

// Linux: 1. Remove lines in WildMagic5/makefile.wm5 to compile only LibCore and LibMathematics
// Linux: 2. Modify "CFLAGS := -c -D__LINUX__" in WildMagic5/LibCore/makeprj.wm5 and WildMagic5/LibMathematics/makeprj.wm5 to "CFLAGS := -c -D__LINUX__ -fPIC"
// Linux: 3. Compile in RELEASE mode: make CFG=Release -f makefile.wm5
// Linux: mex -I. -I../../../extern/mex/include/wildmagic-5.7/include -L../../../extern/mex/lib -lWm5Core -lWm5Mathematics lengthBezier.cpp

double lengthBezier(double *cP, double t0, double t1, int cPdim, int nCP);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Check number of input arguments
    if (nrhs != 1 && nrhs !=3) mexErrMsgTxt("1 or 3 input arguments required.");
    
    // Check number of output arguments
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    
    // Check the dimension of the control points
    const mwSize *inputSize = mxGetDimensions(prhs[0]);
    int nCP = inputSize[0]; // Number of control points
    int cPdim = inputSize[1]; // Control point dimension
    
    if (cPdim > 3 || cPdim < 1) mexErrMsgTxt("Only 1D, 2D, and 3D control points supported!");
	if (nCP < 2) mexErrMsgTxt("At least two control points needed!");

	// Variables
    double *cP; // Control points
    double *length; // Length of the Bezier curve
    double t0, t1; // Parametrization interval
    
    // Associate inputs
    cP = mxGetPr(prhs[0]);
    
    if (nrhs == 3) {
        t0 = *mxGetPr(prhs[1]);
        t1 = *mxGetPr(prhs[2]);
    } else {
		t0 = 0.0;
		t1 = 1.0;
	}

    // Associate outputs
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    length = mxGetPr(plhs[0]);

    *length = lengthBezier(cP, t0, t1, cPdim, nCP);

}

double lengthBezier(double *cP, double t0, double t1, int cPdim, int nCP) {
	if (cPdim == 1) {
		// Put the control points into a vector array
		Vector2d *cPVector2 = new Vector2d[nCP]; // Control points in a Vector2 array
		for (int i=0;i<nCP;i++) {
			cPVector2[i] = Vector2d(cP[0*nCP+i], 0);
		}

		// Create the Bezier curve, compute its length and set it as output
		BezierCurve2d bezierCurve2 = BezierCurve2d(nCP-1, cPVector2);
		return fabs(bezierCurve2.GetLength(t0, t1));

		// BezierCurve2d accepts responsibility for deleting the input arrays
	} else if (cPdim == 2) {
		// Put the control points into a vector array
		Vector2d *cPVector2 = new Vector2d[nCP]; // Control points in a Vector2 array
		for (int i=0;i<nCP;i++) {
			cPVector2[i] = Vector2d(cP[0*nCP+i], cP[1*nCP+i]);
		}

		// Create the Bezier curve, compute its length and set it as output
		BezierCurve2d bezierCurve2 = BezierCurve2d(nCP-1, cPVector2);
		return fabs(bezierCurve2.GetLength(t0, t1));

		// BezierCurve2d accepts responsibility for deleting the input arrays
	} else if (cPdim == 3) {
		// Put the control points into a vector array
		Vector3d *cPVector3 = new Vector3d[nCP]; // Control points in a Vector3 array
		for (int i=0;i<nCP;i++) {
			cPVector3[i] = Vector3d(cP[0*nCP+i], cP[1*nCP+i], cP[2*nCP+i]);
		}

		// Create the Bezier curve, compute its length and set it as output
		BezierCurve3d bezierCurve3 = BezierCurve3d(nCP-1, cPVector3);
		return fabs(bezierCurve3.GetLength(t0, t1));

		// BezierCurve3d accepts responsibility for deleting the input arrays
	}
}

/*int main() {
    double cP[] = {0,0,0,1,1,1,2,2,2,4,5,9,1,9,4,4,4,-1};
    double t0 = 0.1;
    double t1 = 0.9;
    int cPdim = 3;
    int nCP = 6;
    double len = lengthBezier(cP, t0, t1, cPdim, nCP);
}*/


