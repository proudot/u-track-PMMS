/* 
 * Compilation:
 * Mac/Linux: mex  -I../../../extern/mex/include/lemon-1.2.1 maxWeightedMatching.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\lemon-1.2.1" -output maxWeightedMatching maxWeightedMatching.cpp
 */

#include <iostream>

#include <mex.h>

#include <lemon/smart_graph.h>
#include <lemon/matching.h>

using namespace lemon;
using namespace std;

// plhs[0] = M

// prhs[0] = numVertices
// prhs[1] = edges
// prhs[2] = weights

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Read input arguments
	
  if (nrhs != 3) {
    mexErrMsgTxt("3 input arguments required.");
  }
	
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
	
  int num_nodes = (int) *mxGetPr(prhs[0]);
	
  int ndim = mxGetNumberOfDimensions(prhs[2]);
	
  if (ndim != 2) {
    mexErrMsgTxt("Invalid number of dimensions.");
  }
	
  const mwSize* size = mxGetDimensions(prhs[1]);
	
  int num_edges = size[0];

  if (size[1] != 2) {
    mexErrMsgTxt("Invalid number of columns for the 2nd argument.");
  }

  size = mxGetDimensions(prhs[2]);
	
  if (size[0] != num_edges) {
    mexErrMsgTxt("Invalid number of elements for the 3rd argument.");
  }
	
  SmartGraph g;

  g.reserveNode(num_nodes);
  g.reserveEdge(num_edges);

  for (int i = 0; i < num_nodes; ++i)
    g.addNode();
  
  double* p = mxGetPr(prhs[1]);
  double* q = mxGetPr(prhs[2]);

  typedef SmartGraph::EdgeMap<double> WeightMap;

  WeightMap weight_map(g);
  
  for (int i = 0; i < num_edges; ++i)
    {
      SmartGraph::Node u = SmartGraph::nodeFromId(p[i] - 1);
      SmartGraph::Node v = SmartGraph::nodeFromId(p[num_edges + i] - 1);

      SmartGraph::Edge e = g.addEdge(u,v);

      weight_map[e] = q[i];
    }

  MaxWeightedMatching<SmartGraph, WeightMap> matching(g, weight_map);

  matching.run();

  plhs[0] = mxCreateLogicalMatrix(num_edges, 1);

  mxLogical *r = mxGetLogicals(plhs[0]);

  for(SmartGraph::EdgeIt e(g); e != INVALID; ++e)
    r[SmartGraph::id(e)] = matching.matching(e);

  g.clear();
}
