#include <mex.h>
#include <math.h>
#include <stdio.h>

// Compilation: mex -I. distancePointBezierOld.cpp

#include "rpoly_ak1.cpp" // http://www.akiti.ca/rpoly_ak1_Intro.html

double distPointToPointOnLinBez(double *a, double *b, double t, double *point);
double distPointToPointOnQuadBez(double *a, double *b, double *c, double t, double *point);
double distPointToPointOnCubBez(double *a, double *b, double *c, double *d, double t, double *point);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Check number of input arguments
    if (nrhs != 2 && nrhs != 4) {
        mexErrMsgTxt("2 or 4 input arguments required.");
    }
    
    // Check number of output arguments
    if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    // Check the dimension of the control points
    const mwSize *inputSize = mxGetDimensions(prhs[0]);
    int nCP = inputSize[0]; // Number of control points
    int controlPointsDimension = inputSize[1];
    
    if (controlPointsDimension != 3 || nCP < 1 || nCP > 4) {
        mexErrMsgTxt("Only 2,3 or 4 3D control points supported!");
    }
    
    // Check the dimension of the free point
    inputSize = mxGetDimensions(prhs[1]);
    int nPoints = inputSize[0];
    int pointsDimension = inputSize[1];
    
    if (pointsDimension != 3 || nPoints != 1) {
        mexErrMsgTxt("Distance computation of only one 3D point supported!");
    }
    
    // Variables
    double *cP; // Control points
    double *point; // Point
    double *distanceMin; // Minimal distance between the point and the Bezier curve
    double *tDistMin; // Value of the parameter corresponding to the minimal distance
    double t0, t1; // Parametrization interval
    
    // Associate inputs
    cP = mxGetPr(prhs[0]);
    point = mxGetPr(prhs[1]);
    
    if (nrhs == 4) {
        t0 = *mxGetPr(prhs[2]);
        t1 = *mxGetPr(prhs[3]);
    } else {
        t0 = 0.0;
        t1 = 1.0;
    }
    
    // Associate outputs
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    distanceMin = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    tDistMin = mxGetPr(plhs[1]);
    
    // =====================
    // Distance between two points
    // =====================
    if (nCP == 1) {
        *distanceMin = sqrt(pow(cP[0]-point[0],2)+pow(cP[1]-point[1],2)+pow(cP[2]-point[2],2));
    }
    
    // =====================
    // Linear Bezier curve
    // =====================
    else if (nCP == 2) {
        // B[t_] = (1-t)*P0 + t*P1;
        // B[t_] = a*t + b;
        double a[3], b[3];
        
        for (int i=0;i<3;i++) {
            // a = -P0 + P1
            a[i] = -cP[i*nCP+0]+cP[i*nCP+1];
            // b = P0
            b[i] = cP[i*nCP+0];
        }
        
        // d distance^2 / dt ~
        //   a(b-Q) + a^2t |x
        // + a(b-Q) + a^2t |y
        // + a(b-Q) + a^2t |z
        // = poly = 0
        double term0[3], term1[3];
        
        for (int i=0;i<3;i++) {
            // term0 = a(b-Q)
            term0[i] = a[i]*(b[i]-point[i]);
            // term1 = a^2
            term1[i] = a[i]*a[i];
        }
            
        // Put together the terms with the same degree
        double polyCoeff[2];
        polyCoeff[1] = term0[0]+term0[1]+term0[3];
        polyCoeff[0] = term1[0]+term1[1]+term1[3];
        
        // Find the root
        double outZeroReal = -polyCoeff[1]/polyCoeff[0];
        
        // Evaluate the Bezier curve at the beginning
        double distance_t0 = distPointToPointOnLinBez(a, b, t0, point);
        
        // Evaluate the Bezier curve at the end
        double distance_t1 = distPointToPointOnLinBez(a, b, t1, point);
        
        // Determine the closest end point
        if ( distance_t0 > distance_t1 ) {
            *distanceMin = distance_t1;
            *tDistMin = t1;
        } else {
            *distanceMin = distance_t0;
            *tDistMin = t0;
        }
        
        // Check if the closest point is a point on the curve
        if (outZeroReal > t0 && outZeroReal < t1) {
            double t = outZeroReal;
            double distance_t = distPointToPointOnLinBez(a, b, t, point);
            
            // Determine the closest end point
            if ( *distanceMin > distance_t ) {
                *distanceMin = distance_t;
                *tDistMin = t;
            }
        }
        
        // Set the output value
        *distanceMin = sqrt(*distanceMin);
    }
    
    // =====================
    // Quadratic Bezier curve
    // =====================
    else if (nCP == 3) {
        // B[t_] = (1-t)^2*P0 + 2*(1-t)*t*P1 + t^2*P2;
        // B[t_] = a*t^2 + b*t + c;
        double a[3], b[3], c[3];
        
        for (int i=0;i<3;i++) {
            // a = P0-2P1+P2
            a[i] = cP[i*nCP+0]-2*cP[i*nCP+1]+cP[i*nCP+2];
            // b = -2P0+2P1
            b[i] = -2*cP[i*nCP+0]+2*cP[i*nCP+1];
            // c = P0
            c[i] = cP[i*nCP+0];
        }
        
        // d distance^2 / dt ~
        //   bc-bQ + (b^2+2ac-2aQ)t + 3abt^2 + 2a^2t^3 |x
        // + bc-bQ + (b^2+2ac-2aQ)t + 3abt^2 + 2a^2t^3 |y
        // + bc-bQ + (b^2+2ac-2aQ)t + 3abt^2 + 2a^2t^3 |z
        // = poly = 0
        double term0[3], term1[3], term2[3], term3[3];

        for (int i=0;i<3;i++) {
            // term0 = bc-bQ
            term0[i] = b[i]*c[i]-b[i]*point[i];
            // term1 = b^2+2ac-2aQ
            term1[i] = b[i]*b[i]+2*a[i]*c[i]-2*a[i]*point[i];
            // term2 = 3ab
            term2[i] = 3*a[i]*b[i];
            // term3 = 2a^2
            term3[i] = 2*a[i]*a[i];
        }
            
        // Put together the terms with the same degree
        double polyCoeff[4];
        polyCoeff[3] = term0[0]+term0[1]+term0[2];
        polyCoeff[2] = term1[0]+term1[1]+term1[2];
        polyCoeff[1] = term2[0]+term2[1]+term2[2];
        polyCoeff[0] = term3[0]+term3[1]+term3[2];
        
        // Polynomial root solver
        int degree = 3;
        double outZeroReal[3];
        double outZeroImag[3];
        int info[3];
        rpoly_ak1(polyCoeff, &degree, outZeroReal, outZeroImag);
        
        // Evaluate the Bezier curve at the beginning
        double distance_t0 = distPointToPointOnQuadBez(a, b, c, t0, point);
        
        // Evaluate the Bezier curve at the end
        double distance_t1 = distPointToPointOnQuadBez(a, b, c, t1, point);
        
        // Determine the closest end point
        if ( distance_t0 > distance_t1 ) {
            *distanceMin = distance_t1;
            *tDistMin = t1;
        } else {
            *distanceMin = distance_t0;
            *tDistMin = t0;
        }
        
        // Check if the closest point is a point on the curve
        for (int i=0;i<degree;i++) {
            if (outZeroReal[i] > t0 && outZeroReal[i] < t1 && outZeroImag[i] == 0) {
                double t = outZeroReal[i];
                double distance_t = distPointToPointOnQuadBez(a, b, c, t, point);
                
                // Determine the closest end point
                if ( *distanceMin > distance_t ) {
                    *distanceMin = distance_t;
                    *tDistMin = t;
                }
            }
        }
        
        // Set the output value
        *distanceMin = sqrt(*distanceMin);
    }
    
    // =====================
    // Cubic Bezier curve
    // =====================
    else if (nCP == 4) {
        // B[t_] = (1-t)^3*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3
        // B[t_] = a*t^3 + b*t^2 + c*t + d;
        double a[3], b[3], c[3], d[3];
        
        for (int i=0;i<3;i++) {
            // a = -P0 + 3P1 - 3P2 + P3
            a[i] = -cP[i*nCP+0]+3*cP[i*nCP+1]-3*cP[i*nCP+2]+cP[i*nCP+3];
            // b = 3P0 - 6P1 + 3P2
            b[i] = 3*cP[i*nCP+0]-6*cP[i*nCP+1]+3*cP[i*nCP+2];
            // c = -3P0 + 3P1
            c[i] = -3*cP[i*nCP+0]+3*cP[i*nCP+1];
            // d = P0
            d[i] = cP[i*nCP+0];
        }

        // d distance^2 / dt ~
        //   cd-cQ + (c^2+2bd-2bQ)t + (3bc+3ad-3aQ)t^2 + (2b^2+4ac)t^3 + 5abt^4 + 3a^2t^5 |x
        // + cd-cQ + (c^2+2bd-2bQ)t + (3bc+3ad-3aQ)t^2 + (2b^2+4ac)t^3 + 5abt^4 + 3a^2t^5 |y
        // + cd-cQ + (c^2+2bd-2bQ)t + (3bc+3ad-3aQ)t^2 + (2b^2+4ac)t^3 + 5abt^4 + 3a^2t^5 |z
        // = poly = 0
        double term0[3], term1[3], term2[3], term3[3], term4[3], term5[3];

        for (int i=0;i<3;i++) {
            // term0 = cd-cQ
            term0[i] = c[i]*d[i]-c[i]*point[i];
            // term1 = c^2+2bd-2bQ
            term1[i] = c[i]*c[i]+2*b[i]*d[i]-2*b[i]*point[i];
            // term2 = 3bc+3ad-3aQ
            term2[i] = 3*b[i]*c[i]+3*a[i]*d[i]-3*a[i]*point[i];
            // term3 = 2b^2+4ac
            term3[i] = 2*b[i]*b[i]+4*a[i]*c[i];
            // term4 = 5ab
            term4[i] = 5*a[i]*b[i];
            // term5 = 3a^2
            term5[i] = 3*a[i]*a[i];
        }

        // Put together the terms with the same degree
        double polyCoeff[6];
        polyCoeff[5] = term0[0]+term0[1]+term0[2];
        polyCoeff[4] = term1[0]+term1[1]+term1[2];
        polyCoeff[3] = term2[0]+term2[1]+term2[2];
        polyCoeff[2] = term3[0]+term3[1]+term3[2];
        polyCoeff[1] = term4[0]+term4[1]+term4[2];
        polyCoeff[0] = term5[0]+term5[1]+term5[2];
        
        // Polynomial root solver
        int degree = 5;
        double outZeroReal[5];
        double outZeroImag[5];
        int info[5];
        rpoly_ak1(polyCoeff, &degree, outZeroReal, outZeroImag);
        
        // Evaluate the Bezier curve at the beginning
        double distance_t0 = distPointToPointOnCubBez(a, b, c, d, t0, point);
        
        // Evaluate the Bezier curve at the end
        double distance_t1 = distPointToPointOnCubBez(a, b, c, d, t1, point);
        
        // Determine the closest end point
        if ( distance_t0 > distance_t1 ) {
            *distanceMin = distance_t1;
            *tDistMin = t1;
        } else {
            *distanceMin = distance_t0;
            *tDistMin = t0;
        }
        
        // Check if the closest point is a point on the curve
        for (int i=0;i<degree;i++) {
            if (outZeroReal[i] > t0 && outZeroReal[i] < t1 && outZeroImag[i] == 0) {
                double t = outZeroReal[i];
                double distance_t = distPointToPointOnCubBez(a, b, c, d, t, point);
                
                // Determine the closest end point
                if ( *distanceMin > distance_t ) {
                    *distanceMin = distance_t;
                    *tDistMin = t;
                }
            }
        }

        // Set the output value
        *distanceMin = sqrt(*distanceMin);
    }
}

// Computes the distance between a point and a point on a linear bezier curve
double distPointToPointOnLinBez(double *a, double *b, double t, double *point) {
    double point_t[3];
    
    // Evaluate the curve
    for (int i=0;i<3;i++)
        point_t[i] = a[i]*t + b[i];
    
    // Compute the distance
    return pow((point[0]-point_t[0]), 2) + pow((point[1]-point_t[1]), 2) + pow((point[2]-point_t[2]), 2);
}

// Computes the distance between a point and a point on a quadratic bezier curve
double distPointToPointOnQuadBez(double *a, double *b, double *c, double t, double *point) {
    double point_t[3];
    
    // Evaluate the curve
    for (int i=0;i<3;i++)
        point_t[i] = a[i]*pow(t, 2) + b[i]*t + c[i];
    
    // Compute the distance
    return pow((point[0]-point_t[0]), 2) + pow((point[1]-point_t[1]), 2) + pow((point[2]-point_t[2]), 2);

}

// Computes the distance between a point and a point on a cubic bezier curve
double distPointToPointOnCubBez(double *a, double *b, double *c, double *d, double t, double *point) {
    double point_t[3];
    
    // Evaluate the curve
    for (int i=0;i<3;i++)
        point_t[i] = a[i]*pow(t, 3) + b[i]*pow(t, 2) + c[i]*t + d[i];
    
    // Compute the distance
    return pow((point[0]-point_t[0]), 2) + pow((point[1]-point_t[1]), 2) + pow((point[2]-point_t[2]), 2);
}

