/* 
 * Compilation:
 * Mac/Linux: mex -I. -I../../mex/include/c++/ bresenhamMEX.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"."  -I"..\..\mex\include\c++"  -output bresenhamMEX bresenhamMEX.cpp
 */

#include <iterator>
#include <vector>

#include <mex.h>

#include <bresenham.hpp>

template <typename Container>
class inserter
{
  typedef typename Container::value_type value_type;

public:
  inserter(Container & c) : c_(c) {}

  void operator()(value_type v) { c_.push_back(v); }

private:
  Container & c_;
};

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /////////////////////////////////////////////////
  // Check number of input and output parameters //
  /////////////////////////////////////////////////

  if (nrhs < 2)
    mexErrMsgTxt("Two input arguments required.");

  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  ////////////////////////////
  // Check input parameters //
  ////////////////////////////

  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    mexErrMsgTxt("Invalid dimension for p0.");

  if (mxGetNumberOfElements(prhs[0]) != 2)
    mexErrMsgTxt("p0 should be a 2-dimensional point.");

  if (mxGetNumberOfDimensions(prhs[1]) != 2)
    mexErrMsgTxt("Invalid dimension for p1.");

  if (mxGetNumberOfElements(prhs[1]) != 2)
    mexErrMsgTxt("p1 should be a 2-dimensional point.");

  //////////////////////////
  // Get input parameters //
  //////////////////////////

  double* ptr = mxGetPr(prhs[0]);

  vector<2,double> p0, p1;

  p0[0] = *ptr++;
  p0[1] = *ptr;

  ptr = mxGetPr(prhs[1]);

  p1[0] = *ptr++;
  p1[1] = *ptr;
  
  ///////////////////////
  // Compute bresenham //
  ///////////////////////

  std::vector< vector<2,double> > res;
  inserter< std::vector<vector<2,double> > > f(res);
  bresenham(p0, p1, f);

  ///////////////////////////
  // Output list of points //
  ///////////////////////////

  if (nlhs > 0)
    {
      unsigned n = res.size();

      plhs[0] = mxCreateDoubleMatrix(n, 2, mxREAL);
      ptr = mxGetPr(plhs[0]);

      for (unsigned i = 0; i < n; ++i)
	{
	  const vector<2,double> & v = res[i];
	  *ptr = v[0];
	  *(ptr++ + n) = v[1];
	}
    }
}
