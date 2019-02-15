#include <mex.h>

#include <boost/graph/adjacency_list.hpp>

#include <image.hpp>
#include <mx_wrapper.hpp>
#include <fast_marching.hpp>

// Type definition for graph vertex
template <int n>
struct vertex_t
{
  // default maximum degree = +inf
  vertex_t() : max_degree(std::numeric_limits<int>::max()) {}
  // the point where U[p] = 0
  vector<n, int> p;
  // maximum degree allowed
  int max_degree;
};

// Type definition for graph edge
template <int n>
struct edge_t
{
  // the saddle point between front r1 and r2
  vector<n, double> sp;
  // the 2 closest grid points to the saddle point belonging
  // respectively to r1 and r2.
  vector<n, int> sp_r1, sp_r2;
  // the minimal meeting time
  double time;
};

template <int n>
static void dispatch(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Graph type definition
  typedef boost::adjacency_list<boost::listS,
    boost::vecS,
    boost::undirectedS,
    vertex_t<n>,
    edge_t<n> > graph_t;

  // Populate graph from initial points
  const mwSize* x_size = mxGetDimensions(prhs[0]);
  double * ptr = mxGetPr(prhs[0]);

  graph_t g(x_size[0]);

  if (nrhs == 3 && !mxIsEmpty(prhs[2]))
    {
      // Check max_degree argument
      if (mxGetM(prhs[2]) != x_size[0])
	mexErrMsgTxt("1st and 3rd argument must have the same number of \
columns.");

      if (mxGetN(prhs[2]) != 1)
	mexErrMsgTxt("3rd argument must be a column vector.");

      double * deg_ptr = mxGetPr(prhs[2]);

      for (unsigned k = 0; k < num_vertices(g); ++k)
	g[k].max_degree = (int) deg_ptr[k];
    }

  for (unsigned i = 0; i < n; ++i)
    {
      int offset = x_size[0] * i;

      for (unsigned k = 0; k < num_vertices(g); ++k)
	g[k].p[i] = (int) ptr[k + offset] - 1;
    }

  int f_size[n];

  // get speed function size
  if (nrhs < 2 || mxIsEmpty(prhs[1]))
    {
      for (unsigned i = 0; i < n; ++i)
	{
	  f_size[i] = 0;
	  
	  for (unsigned k = 0; k < num_vertices(g); ++k)
	    {
	      vector<n, int> & p = g[k].p;
	      
	      if (p[i] + 1 > f_size[i]) f_size[i] = p[i] + 1;
	    }
	}     
    }
  else
    {
      if (mxGetNumberOfDimensions(prhs[1]) != n)
	mexErrMsgTxt("1st and 2nd argument's dimensions mismatch.");
      
      sizeWrapper<n>::convert(mxGetDimensions(prhs[1]), f_size);
    }

  fast_marching<n> fm(f_size);

  // Compute front propagation

  image<n, double> f(f_size);
  
  if (nrhs < 2 || mxIsEmpty(prhs[1]))
    f.fill(1);
  else
    f.fill(mxGetPr(prhs[1]));
  
  fm.compute(f, g);

  // Output U
  if (nlhs > 0)
    image2mxArray(fm.u(), plhs[0]);

  // Output E
  if (nlhs > 1)
    {
      typename graph_t::edge_iterator ei, ei_end;

      int m = num_edges(g);

      plhs[1] = mxCreateDoubleMatrix(m, 2, mxREAL);

      ptr = mxGetPr(plhs[1]);

      for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{
	  *ptr = (double) source(*ei, g) + 1;
	  *(ptr + m) = (double) target(*ei, g) + 1;

	  ptr++;
	}
    }

  // Output S
  if (nlhs > 2)
    {
      typename graph_t::edge_iterator ei, ei_end;

      int m = num_edges(g);

      plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);

      ptr = mxGetPr(plhs[2]);

      for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{
	  vector<n, double> & sp = g[*ei].sp;

	  for (int i = 0; i < n; ++i)
	    *(ptr + i * m) = (double) sp[i] + 1.0;

	  ptr++;
	}
    }

  // Output R
  if (nlhs > 3)
    {
      // TODO: add 1 to the map.
      image2mxArray(fm.r(), plhs[3]);
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Check number of input and output parameters

  if (nrhs < 1)
    mexErrMsgTxt("At least one input argument is required.");
  
  if (nrhs > 3)
    mexErrMsgTxt("Too many input arguments.");

  if (nlhs > 4)
    mexErrMsgTxt("Too many output arguments.");

  // Get the dimension

  int dim = mxGetN(prhs[0]);

  // Dispatch according to dimension.

  switch (dim)
    {
    case 2: dispatch<2>(nlhs, plhs, nrhs, prhs); break;
    case 3: dispatch<3>(nlhs, plhs, nrhs, prhs); break;
    default: mexErrMsgTxt("Invalid points dimension (must be 2d or 3d).");
    }
}
