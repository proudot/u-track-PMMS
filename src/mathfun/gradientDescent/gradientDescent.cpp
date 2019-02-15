/* [pts, values] = gradientDescent(F,X,Y);
 *
 * Sylvain Berlemont, 2010 (last modified Aug 3, 2011)
 *
 * Compilation:
 * Mac/Linux: mex -I.  -I../../mex/include/c++ gradientDescent.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"." -I"..\..\mex\include\c++" -output gradientDescent gradientDescent.cpp
 */

#include <mex.h>

#include <vector>

#include <image.hpp>
#include <mx_wrapper.hpp>


static void gradientDescent(const image<2,double> & f,
        const vector<2,double>& start,
        std::vector< vector<2, double> > & pts,
        std::vector<double> & values)
        
{
    static const double eps = 1e-5;
    
    vector<2,double> q, p(start);
    
    double dx, dy;
    
    int i, j;
    
    do
    {
        i = (int) floor(p[0]);
        j = (int) floor(p[1]);
        
        pts.push_back(p);
        values.push_back(f[p]);
        
        q = p;
        
        dx = .5 * (f(p[0] + 1, p[1]) - f(p[0] - 1, p[1]));
        dy = .5 * (f(p[0], p[1] + 1) - f(p[0], p[1] - 1));
        
        p[0] -= dx;
        p[1] -= dy;
    }
    while (f.contains(i, j) && f.contains(i + 1, j + 1) &&
            /*f[p] > eps && */dist(p, q) > eps);
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // Check number of input and output parameters
    
    if (nrhs != 3)
        mexErrMsgTxt("Three input argument are required.");
    
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    
    // Check the argument types
    if (!mxIsDouble(prhs[0]))
    {
        std::string err_msg("Undefined function or method 'gradientDescent' for input arguments of type ");
        err_msg += mxGetClassName(prhs[0]);
        mexErrMsgTxt(err_msg.c_str());
    }
    
    // Check the number of dimensions
    if (mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("1st argument is not a 2-dimensional matrix.");
    
    if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
        mexErrMsgTxt("2nd and 3rd arguments must have the same size.");
    
    // Read the image
    int f_size[2];
    sizeWrapper<2>::convert(mxGetDimensions(prhs[0]), f_size);
    image<2,double> f(f_size,2);
    f.fill(mxGetPr(prhs[0]));
    f.border_replicate(f.margin());
    
    // Compute gradient descent
    int n = mxGetNumberOfElements(prhs[1]);
    
    std::vector< std::vector< vector<2,double> > > pts_list(n);
    std::vector< std::vector<double> > values_list(n);
    
    double* px = mxGetPr(prhs[1]);
    double* py = mxGetPr(prhs[2]);
    
    for (int i = 0; i < n; ++i)
    {
        vector<2,double> pt;
        pt[0] = px[i] - 1;
        pt[1] = py[i] - 1;
        
        int x = (int) floor(pt[0]);
        int y = (int) floor(pt[1]);
        
        if (f.contains(x, y) && f.contains(x + 1, y + 1))
            gradientDescent(f,pt,pts_list[i],values_list[i]);
    }
    
    // Output
    if (nlhs > 0)
    {
        plhs[0] = mxCreateCellMatrix(n,1);
        
        for (int i = 0; i < n; ++i)
        {
            int m = pts_list[i].size();
            
            mxArray* tmp = mxCreateDoubleMatrix(m, 2, mxREAL);
            double* p = mxGetPr(tmp);
            
            for (int j = 0; j < m; ++j)
            {
                const vector<2,double> & pt = pts_list[i][j];
                
                p[j] = pt[0] + 1;
                p[j + m] = pt[1] + 1;
            }
            
            mxSetCell(plhs[0], i, tmp);
        }
    }
    
    if (nlhs > 1)
    {
        plhs[1] = mxCreateCellMatrix(n,1);
        
        for (int i = 0; i < n; ++i)
        {
            int m = values_list[i].size();
            
            mxArray* tmp = mxCreateDoubleMatrix(m, 1, mxREAL);
            double * p = mxGetPr(tmp);
            
            for (int j = 0; j < m; ++j)
                p[j] = values_list[i][j];
            
            mxSetCell(plhs[1], i, tmp);
        }
    }
}
