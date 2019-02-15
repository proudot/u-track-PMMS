/*
 * k-sample Anderson-Darling test
 *
 * Implements the test described in
 * [1] Scholz and Stephens, J. Am. Stat. Assoc. 82(399), 1987
 *
 * Created by Francois Aguet on 03/08/2012.
 */

// Compilation:
// OS X: use Xcode project file
// OS X/Linux: mex -I/usr/local/include -I../../mex/include adtestk.cpp
// Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\mex\include" -output adtestk adtestk.cpp


#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include "mex.h"

using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs<1 || nrhs>2) mexErrMsgTxt("Required input: cell array of samples; optional input: alpha value.");
    if (!mxIsCell(prhs[0])) mexErrMsgTxt("Input must be a cell array of samples");
    int ai = 2; // index for alpha = 0.05
    if (nrhs==2) {
        if (!mxIsDouble(prhs[1])) {
            mexErrMsgTxt("Alpha value must be a double: 0.25, 0.1, 0.05, 0.025, 0.01");
        } else {
            double alpha = *mxGetPr(prhs[1]);
            
            static const double alphaVec[5] = {0.25, 0.1, 0.05, 0.025, 0.01};
            ai = 5;
            for (int i=0;i<5;++i) {
                if ( fabs(alpha-alphaVec[i]) < 1e-10 ) {
                    ai = i;
                    break;
                }
            }
            if (ai==5) {
                mexErrMsgTxt("Admissible alpha values: 0.25, 0.1, 0.05, 0.025, 0.01");
            }
        }
    }
    
    int k = (int)mxGetNumberOfElements(prhs[0]);
    if (k<2) mexErrMsgTxt("Input must contain at least two sample arrays.");
    // get pointers to sample arrays
    mxArray* p;
    
    // size of sample vectors
    vector<int> n(k);
    int N = 0; // total number of samples
    for (int i=0;i<k;++i) {
        n[i] = (int)mxGetNumberOfElements(mxGetCell(prhs[0], i));
        N += n[i];
    }
    
    // sample vectors
    vector< vector<double> > samples;
    samples.resize(k);
    
    // vector of all samples combined
    vector<double> pooledSamples(N);
    
    N = 0; // reset counter
    for (int i=0;i<k;++i) {
        // pointer to cell data
        p = mxGetCell(prhs[0], i);
        double* dp = mxGetPr(p);
        
        // copy sample vectors
        samples[i].resize(n[i]);
        copy(dp, dp+n[i], samples[i].begin());
        sort(samples[i].begin(), samples[i].end()); // will need multiplity later
        
        // pooled vector
        copy(dp, dp+n[i], pooledSamples.begin()+N);
        
        N += n[i];
    }
    
    // unique samples in Z
    sort(pooledSamples.begin(), pooledSamples.end());
    vector<double>::iterator it;
    it = unique(pooledSamples.begin(), pooledSamples.end());
    pooledSamples.resize(it - pooledSamples.begin());
    int L = (int)pooledSamples.size();
    
    // Variance (sigma_N^2, eq. 4): depends only on sample sizes
    double H = 0.0;
    for (int i=0;i<k;++i) {
        H += 1.0/n[i];
    }
    double h = 0.0;
    for (int i=1;i<N;++i) {
        h += 1.0/i;
    }
    double g = 0.0;
    for (int i=1;i<=N-2;++i) {
        for (int j=i+1;j<N;++j) {
            g += 1.0/((N-i)*j);
        }
    }
    int k2 = k*k;
    double a = (4.0*g-6.0)*(k-1.0) + (10.0-6.0*g)*H;
    double b = (2.0*g-4.0)*k2 + 8.0*h*k + (2.0*g-14.0*h-4.0)*H - 8.0*h + 4.0*g - 6.0;
    double c = (6.0*h+2.0*g-2.0)*k2 + (4.0*h-4.0*g+6.0)*k + (2.0*h-6.0)*H + 4.0*h;
    double d = (2.0*h+6.0)*k2 - 4.0*h*k;
    double varA = (((a*N + b)*N + c)*N + d) / ( (N-1.0)*(N-2.0)*(N-3.0) );
    
    // loop through current sample; if value not equal pooled sample value, advance pooled sample index j
    vector< vector<int> > f(k);
    int j;
    for (int i=0;i<k;++i) {
        f[i].resize(L);
        j = 0;
        for (int ki=0;ki<n[i];++ki) {
            while (pooledSamples[j]!=samples[i][ki])
                j++;
            //can increment count
            f[i][j]++;
        }
    }
    // multiplicity of Z* in Z
    vector<int> l(L, 0);
    for (int j=0;j<L;++j) {
        for (int i=0;i<k;++i) {
            l[j] += f[i][j];
        }
    }
    
    vector<double> Ba(L);
    Ba[0] = l[0];
    for (int j=1;j<L;++j) {
        Ba[j] = Ba[j-1]+l[j];
    }
    for (int j=0;j<L;++j) {
        Ba[j] -= l[j]/2.0;
    }
    
    
    // Compute test statistic (eq. 7)
    double A2 = 0.0;
    double t1, t2;
    vector<double> Mai(L);
    for (int i=0;i<k;++i) {
        // Compute M_ai for ith sample
        Mai[0] = f[i][0];
        for (int j=1;j<L;++j) {
            Mai[j] = Mai[j-1]+f[i][j];
        }
        for (int j=0;j<L;++j) {
            Mai[j] -= f[i][j]/2.0;
        }
        
        t2 = 0.0;
        for (int j=0;j<L;++j) {
            t1 = N*Mai[j]-n[i]*Ba[j];
            t2 += l[j]*t1*t1/(Ba[j]*(N-Ba[j])-N*l[j]/4.0);
        }
        t2 /= n[i];
        A2 += t2;
    }
    A2 *= (N-1.0)/(N*N);
    
    // test statistic
    double T = (A2 - (k-1.0)) / sqrt(varA); // p. 921
    
    int m[] = {1, 2, 3, 4, 6, 8, 10};
    bool lookup = false;
    int mi;
    for (mi=0;mi<7;++mi) {
        if (m[mi]==k-1) {
            lookup = true;
            break;
        }
    }
    double cval;
    if (lookup) {
        // Table I
        static const double cTable[][5] = {{0.326, 1.225, 1.960, 2.719, 3.752},
                                           {0.449, 1.309, 1.945, 2.576, 3.414},
                                           {0.498, 1.324, 1.915, 2.493, 3.246},
                                           {0.525, 1.329, 1.894, 2.438, 3.139},
                                           {0.557, 1.332, 1.859, 2.365, 3.005},
                                           {0.576, 1.330, 1.839, 2.318, 2.920},
                                           {0.590, 1.329, 1.823, 2.284, 2.862}};
        cval = cTable[mi][ai];
    } else { // interpolate critical value
        static const double b[][3] = {{0.675,-0.245, -0.105},
                                      {1.281, 0.250, -0.305},
                                      {1.645, 0.678, -0.362},
                                      {1.960, 1.149, -0.391},
                                      {2.326, 1.822, -0.396}};
        cval = b[ai][0] + b[ai][1]/sqrt(k-1.0) + b[ai][2]/(k-1.0);
    }
    int hval = T>=cval;
    
    // Return results to Matlab: hval, T, A2, cval
    if (nlhs>0) {
        plhs[0] = mxCreateDoubleScalar(hval);
    }
    if (nlhs>1) {
        plhs[1] = mxCreateDoubleScalar(T);
    }
    if (nlhs>2) {
        plhs[2] = mxCreateDoubleScalar(A2);
    }
    if (nlhs>3) {
        plhs[3] = mxCreateDoubleScalar(cval);
    }
}

