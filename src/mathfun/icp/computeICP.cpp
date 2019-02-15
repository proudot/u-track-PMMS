 /* [T R] = computeICP(X1, X2, numIter, tol)
 *
 * (c) Sylvain Berlemont, 2011 (last modified Oct 7, 2011)
 *
 * Compilation:
 * Mac/Linux: mex -I. -I../../deprecated/kdtree -I../../mex/include/c++ computeICP.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP" -I"." -I"..\..\deprecated\kdtree" -I"..\..\mex\include\c++" -output computeICP computeICP.cpp
 */

# include <mex.h>

# include <vector>

# include <vector.hpp>
# include <matrix.hpp>
# include <quaternion.hpp>
# include <KDTree.hpp>
# include <jacobi.hpp>

// Compilation line:
// mex -I. -I../kdtree -I../../mex/include/c++ computeICP.cpp

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Check input/output parameter

  if (nrhs != 4)
    mexErrMsgTxt("Four input arguments required.");

  if (nlhs > 2)
    mexErrMsgTxt("Too many output arguments.");

  int dim = mxGetN(prhs[0]);

  if (dim != 3 || dim != mxGetN(prhs[1]))
    mexErrMsgTxt("Input arguments must be a row vectors of 3D points.");

  // Get the parameters

  int n_model = mxGetM(prhs[0]);
  int n_data = mxGetM(prhs[1]);

  std::vector< vector<3, double> > X, P;

  double *ptr = mxGetPr(prhs[0]);
  vector<3, double> v;

  for (int i = 0; i < n_model; ++i)
    {
      v[0] = ptr[i];
      v[1] = ptr[n_model + i];
      v[2] = ptr[2 * n_model + i];
      X.push_back(v);
    }

  ptr = mxGetPr(prhs[1]);
  for (int i = 0; i < n_data; ++i)
    {
      v[0] = ptr[i];
      v[1] = ptr[n_data + i];
      v[2] = ptr[2 * n_data + i];
      P.push_back(v);
    }

  int max_iter = (int) *mxGetPr(prhs[2]);
	
  double tol = *mxGetPr(prhs[3]);

  // Build the kd-tree
  KDTree<3, double> kdtree(X);

  // Resolve ICP

  static const double eps = std::numeric_limits<double>::epsilon();
  
  std::vector< vector<3, double> > Y(n_data), Pk(P);

  vector<3, double> a, mu_Y, mu_P, qt;
  quaternion qr(1, 0, 0, 0);
  matrix<3, 3, double> Spy;
  matrix<4, 4, double> Q;

  int iter = 0;
	
  // Compute mu_P
  for (int i = 0; i < n_data; ++i)
    mu_P += P[i];
  mu_P /= (double) n_data;

  // 1. Compute the the closest point Y = C(Pk, X)
  double dk_old = std::numeric_limits<double>::max();
  double dk = 0;
  for (int i = 0; i < n_data; ++i)
    {
      KDTree<3,double>::pair_type pair = kdtree.closest_point(Pk[i]);

      dk += pair.first;
      Y[i] = X[pair.second];
    }
  dk /= (double) n_data;
  
  while (iter++ < max_iter && fabs(dk) > eps && (dk_old - dk) / dk > tol)
    {
      // 2. Compute the registration (qk, dk) = Q(P0, Yk)

      // Compute mu_Y
      mu_Y.set_all(0);
      for (int i = 0; i < n_data; ++i)
	mu_Y += Y[i];
      mu_Y /= (double) n_data;

      // Compute Spy
      Spy.set_all(0);
      for (int i = 0; i < n_data; ++i)
	Spy += mtimes(P[i], X[i]);
      Spy /= (double) n_data;
      Spy -= mtimes(mu_P, mu_Y);

      // Compute a
      a[0] = Spy(1,2) - Spy(2,1);
      a[1] = Spy(2,0) - Spy(0,2);
      a[2] = Spy(0,1) - Spy(1,0);

      // Compute Q
      Q(0, 0) = Spy.trace();

      for (int i = 1; i < 4; ++i)
	{
	  Q(i, 0) = a[i - 1];
	  Q(0, i) = a[i - 1];
	}

      for (int i = 1; i < 4; ++i)
	Q(i, i) = 2 * Spy(i - 1, i - 1) - Q(0, 0);

      Q(1,2) = Spy(0,1) + Spy(1,0);
      Q(2,1) = Spy(0,1) + Spy(1,0);
      Q(1,3) = Spy(0,2) + Spy(2,0);
      Q(3,1) = Spy(0,2) + Spy(2,0);
      Q(2,3) = Spy(1,2) + Spy(2,1);
      Q(3,2) = Spy(1,2) + Spy(2,1);

      // Find the eigen vector associated with the largest eigen value
      // of Q.
      qr = jacobi(Q);

      qt = mu_Y - mu_P;

      // 3. Apply the registration Pk+1 = qk(P0)

      for (int i = 0; i < n_data; ++i)
	Pk[i] = qt + qr.rotate(P[i] - mu_P) + mu_P;
      
      // 1. Compute the closest point Y = C(Pk, X) and the new
      // distance between the 2 sets.

      dk_old = dk;

      dk = 0;
      for (int i = 0; i < n_data; ++i)
	{
	  KDTree<3,double>::pair_type pair = kdtree.closest_point(Pk[i]);

	  dk += pair.first;
	  Y[i] = X[pair.second];
	}
      dk /= (double) n_data;
    }

  // Return T and R

  if (nlhs > 0)
    {
      plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
      ptr = mxGetPr(plhs[0]);
      for (unsigned i = 0; i < 3; ++i)
	ptr[i] = qt[i];

      if (nlhs > 1)
	{
	  plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
	  ptr = mxGetPr(plhs[1]);
 	  matrix<3,3,double> mat = qr.to_matrix();
	  unsigned c = 0;
	  for (unsigned j = 0; j < 3; ++j)
	    for (unsigned i = 0; i < 3; ++i)
	      ptr[c++] = mat(i, j);
	}	
		
    }
}
