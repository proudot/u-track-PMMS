#include <mex.h>
#include <math.h>
#include <stdio.h>

#include "sisl.h"

// Windows: mex -I. -I..\..\..\extern\mex\include\sisl-4.5.0\include -L..\..\..\extern\mex\lib -lsisl distancePointBezier.cpp
// Linux: mex -I. -I../../../extern/mex/include/sisl-4.5.0/include -L../../../extern/mex/lib -lsisl distancePointBezier.cpp

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd);
void extend(double *cP, int cPdim, int nCP, double tStart, double tEnd);
void transposeArray(double *a, int sizeX, int sizeY);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// =====================
	// Check inputs
	// =====================

	// Check number of input arguments
	if (nrhs != 2 && nrhs != 4) mexErrMsgTxt("2 or 4 input arguments required.");

	// Check number of output arguments
	if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

	// Check the dimension of the control points
	const mwSize *inputSize = mxGetDimensions(prhs[0]);
	int nCP = inputSize[0]; // Number of control points
	int cPdim = inputSize[1];
	if (cPdim < 1 || nCP < 1) mexErrMsgTxt("At least 1 1D control point is needed!");

	// Check the dimension of the free point
	inputSize = mxGetDimensions(prhs[1]);
	int nPoints = inputSize[0];
	int pointsDimension = inputSize[1];
	if (nPoints != 1) mexErrMsgTxt("Distance computation of only one point supported!");
	if (pointsDimension != cPdim) mexErrMsgTxt("The dimension of the point and the control points should be the same!");

	// =====================
	// Inputs/Outputs
	// =====================

	// Variables
	double *cP; // Control points
	double *point; // Point
	double *distanceMin; // Minimal distance between the point and the Bezier curve
	double *tDistMin; // Value of the parameter corresponding to the minimal distance

	// Associate inputs
	double *cPtmp = mxGetPr(prhs[0]);
	point = mxGetPr(prhs[1]);

	// Transpose control points
	cP = new double[cPdim*nCP];
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = cPtmp[i];
	transposeArray(cP, nCP, cPdim);

	// Associate outputs
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	distanceMin = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	tDistMin = mxGetPr(plhs[1]);

	// Optional inputs
	if (nrhs == 4) {
		double tStart, tEnd;
		tStart = *mxGetPr(prhs[2]);
		tEnd = *mxGetPr(prhs[3]);

		if (tStart >= 0 && tEnd <=1) {
			truncate(cP, cPdim, nCP, tStart, tEnd);
		} else {
			extend(cP, cPdim, nCP, tStart, tEnd);
		}

	}

	// =====================
	// Compute the distance
	// =====================

	if (nCP == 1) { // Single control point
		double dist = 0;
		for (int i = 0; i < cPdim; i++) {
			dist = dist + pow(cP[i] - point[i],2);
		}
		*distanceMin = sqrt(dist);
		*tDistMin = 0;
	} else { // A curve

		// Convert the Bézier curve into a B-spline
		const int number = nCP;
		const int order = nCP;

		int nKnots = number+order;

		// Define the knots of the B-spline
		double *knots = new double[nKnots];
		for (int i=0;i<nKnots/2;i++) {
			knots[i] = 0;
		}
		for (int i=nKnots/2;i<nKnots;i++) {
			knots[i] = 1;
		}

		// Inputs
		SISLCurve *curve = newCurve(number, // number of control points
			order,							// order of spline curve (degree + 1)
			knots,							// pointer to knot vector (parametrization)
			cP,								// pointer to coefficient vector (control points)
			3,								// kind = Polynomial Bezier curve
			cPdim,							// dimension
			0);								// no copying of information, 'borrow' arrays
		double epsge = 1.0e-5;				// geometric tolerance
		double epsco = 0;					// Not used

		// Outputs
		int numintpt; 
		double *intpar;
		int numintcu;
		SISLIntcurve **intcurve;
		int jstat = 0; // Status flag

		// Compute the curve parameter of the closest point (See the SISL manual)
		s1953(curve,point,cPdim,epsco,epsge,&numintpt,&intpar,&numintcu,&intcurve,&jstat);
		if (jstat < 0) mexErrMsgTxt("SISL library error in function s1953!");

		if (numintpt > 0 && numintcu == 0) {
			*tDistMin = *intpar;
			delete intpar;
		} else if (numintpt == 0 && numintcu > 0) { // The libary retruns a intersection curve object
			if ((*intcurve)->ipar1 == 1 && (*intcurve)->ipar2 == 0) {
				int nPoints = (*intcurve)->ipoint;
				*tDistMin = 0;
				for (int i=0;i < nPoints; i++) {
					*tDistMin += (*intcurve)->epar1[i];
				}
				*tDistMin /= nPoints;
			} else {
				mexErrMsgTxt("These values are suspicious!");
			}
			freeIntcurve(*intcurve);
		} else {
			mexErrMsgTxt("These values are suspicious!");
		}

		// Inputs
		int leftknot = nKnots/2 - 1;

		// Evaluate the curve (See the SISL manual)
		double *curvePoint = new double[cPdim];
		s1227(curve, 0, *tDistMin, &leftknot, curvePoint, &jstat);
		if (jstat < 0) mexErrMsgTxt("SISL library error in function s1227!");

		// Compute the distance to the closest point
		double dist = 0;
		for (int i = 0; i<cPdim; i++) {
			dist = dist + pow(curvePoint[i]-point[i],2);
		}

		*distanceMin = sqrt(dist);

		// Clean up
		freeCurve(curve);
		
		delete [] knots;
		delete [] curvePoint;
		delete [] cP;
	}
}

void truncate(double *cP, int cPdim, int nCP, double tStart, double tEnd) {
	// Convert the Bézier curve into a B-spline
	const int number = nCP;
	const int order = nCP;
	const int nKnots = number+order;

	// Define the knots of the B-spline
	double *knots = new double[nKnots];
	for (int i=0;i<nKnots/2;i++) knots[i] = 0;
	for (int i=nKnots/2;i<nKnots;i++) knots[i] = 1;

	// Inputs
	SISLCurve *curve = newCurve(number, // number of control points
		order,							// order of spline curve (degree + 1)
		knots,							// pointer to knot vector (parametrization)
		cP,								// pointer to coefficient vector (control points)
		3,								// kind = Polynomial Bezier curve
		cPdim,							// dimension
		0);								// no copying of information, 'borrow' arrays

	// Outputs
	int jstat = 0;	// Status flag
	SISLCurve *newcurve = newCurve(number,	// number of control points
		order,								// order of spline curve (degree + 1)
		knots,								// pointer to knot vector (parametrization)
		cP,									// pointer to coefficient vector (control points)
		3,									// kind = Polynomial Bezier curve
		cPdim,								// dimension
		0);									// no copying of information, 'borrow' arrays

	// Truncate the curve (See SISL manual)
	s1712(curve, tStart, tEnd, &newcurve, &jstat);
	if (jstat < 0) mexErrMsgTxt("SISL library error in function s1712!");

	// Copy the control points
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = newcurve->ecoef[i];

	// Clean up
	freeCurve(curve);
	freeCurve(newcurve);
	delete [] knots;
}

void extend(double *cP, int cPdim, int nCP, double tStart, double tEnd) {
	// Convert the Bézier curve into a B-spline
	const int number = nCP;
	const int order = nCP;
	const int nKnots = number+order;

	// Define the knots of the B-spline
	double *knots = new double[nKnots];
	for (int i=0;i<nKnots/2;i++) knots[i] = 0;
	for (int i=nKnots/2;i<nKnots;i++) knots[i] = 1;

	// Inputs
	SISLCurve *curve = newCurve(number, // number of control points
		order,							// order of spline curve (degree + 1)
		knots,							// pointer to knot vector (parametrization)
		cP,								// pointer to coefficient vector (control points)
		3,								// kind = Polynomial Bezier curve
		cPdim,							// dimension
		0);								// no copying of information, 'borrow' arrays

	// Sample the curve
	int leftknot = nKnots/2 - 1;
	double *curvePoints = new double[cPdim*nCP];
	int jstat = 0;

	double *epar = new double[nCP];
	for (int i=0; i < nCP; i++) {
		epar[i] = (tStart <= tEnd) ? i*(tEnd-tStart)/(nCP-1)+tStart : i*(tStart-tEnd)/(nCP-1)+tEnd;

		// Evaluate the curve (See the SISL manual)
		s1227(curve, 0, epar[i], &leftknot, &curvePoints[i*cPdim], &jstat);
		if (jstat < 0) mexErrMsgTxt("SISL library error in function s1227!");
	}

	// Interpolate these points to find the new curve
	// Inputs
	int *ntype = new int[nCP];
	for (int i=0;i<nCP;i++) ntype[i] = 1;
	int icnsta = 0;
	int icnend = 0;
	int iopen = 1;
	int ik = order;
	double astpar = (tStart <= tEnd) ? tStart : tEnd;

	// Outputs
	double cendpar;
	double *gpar = new double[nCP];
	int jnbpar;

	s1357(curvePoints, nCP, cPdim, ntype, epar, icnsta, icnend, iopen, ik, astpar, &cendpar, &curve, &gpar, &jnbpar, &jstat);
	if (jstat < 0) mexErrMsgTxt("SISL library error in function s1357!");

	// Reverse the orientation if needed
	if (tEnd <= tStart){
		s1706(curve);
	}

	// Copy the control points
	for (int i = 0; i < cPdim*nCP; i++) cP[i] = curve->ecoef[i];

	delete [] curvePoints;
	delete [] ntype;
	delete [] epar;
	delete [] gpar;
	freeCurve(curve);
}

void transposeArray(double *a, int sizeX, int sizeY) {
	double *aTrans = new double[sizeX * sizeY];
	for (int x = 0; x < sizeX; x++) {
		for (int y = 0; y < sizeY; y++) {
			aTrans[x * sizeY + y] = a[y * sizeX + x];
		}
	}
	for (int i = 0; i < sizeX * sizeY; i++) a[i] = aTrans[i];
	delete [] aTrans;
}