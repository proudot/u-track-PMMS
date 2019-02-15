# ifndef JACOBI_HPP
# define JACOBI_HPP

# include <vector.hpp>
# include <matrix.hpp>

// Note: this comes from eigen_sym.h from Numerical Recipes in
// C++. pp574-575.

static inline
void rot(matrix<4,4,double> & a, double s, double tau,
	 int i, int j, int k, int l)
{
  double g = a(i, j);
  double h = a(k, l);

  a(i, j) = g - s * (h + g * tau);
  a(k, l) = h + s * (g - h * tau);
}

vector<4, double> jacobi(matrix<4,4,double> & a)
{
  const static double eps = 1e-12;

  double d[4];
  matrix<4,4,double> v = matrix<4,4,double>::Id;
  double theta, tau, t, sm, s, h, g, c, b[4], z[4];

  for (int ip = 0; ip < 4; ip++)
    {
      b[ip] = d[ip] = a(ip, ip);
      z[ip] = 0.0;
    }

  for (int i = 1; i <= 50; ++i)
    {
      sm = 0.0;
      for (int ip = 0; ip < 3; ip++)
	for (int iq = ip + 1; iq < 4; iq++)
	  sm += fabs(a(ip, iq));

      if (sm < eps)
	{
	  double dd = d[0];
	  int iq = 0;

	  for (int ip = 1; ip < 4; ip++)
	    if (d[ip] > dd)
	      {
		iq = ip;
		dd = d[ip];
	      }

	  vector<4,double> res;

	  for (unsigned j = 0; j < 4; ++j)
	    res[j] = v(j, iq);

	  return res.normalize();
	}

      double tresh = i < 4 ? 0.0125 * sm : 0.0;

      for (int ip = 0; ip < 3; ip++)
	{
	  for (int iq = ip + 1; iq < 4; iq++)
	    {
	      g = 100.0 * fabs(a(ip, iq));

	      if (i > 4 && g <= eps * fabs(d[ip]) && g <= eps * fabs(d[iq]))
		a(ip, iq) = 0.0;
	      else if (fabs(a(ip, iq)) > tresh)
		{
		  h = d[iq] - d[ip];

		  if (g <= eps * fabs(h))
		    t = (a(ip, iq)) / h;
		  else
		    {
		      theta = 0.5 * h / (a(ip, iq));
		      t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
		      if (theta < 0.0)
			t = -t;
		    }

		  c = 1.0 / sqrt(1 + t * t);
		  s = t * c;
		  tau = s / (1.0 + c);
		  h = t * a(ip, iq);
		  z[ip] -= h;
		  z[iq] += h;
		  d[ip] -= h;
		  d[iq] += h;
		  a(ip, iq) = 0.0;

		  for (int j = 0; j < ip; j++)
		    rot(a, s, tau, j, ip, j, iq);
		  for (int j = ip + 1; j <= iq - 1; j++)
		    rot(a, s, tau, ip, j, j, iq);
		  for (int j = iq + 1; j < 4; j++)
		    rot(a, s, tau, ip, j, iq, j);
		  for (int j = 0; j < 4; j++)
		    rot(v, s, tau, j, ip, j, iq);
		}
	    }
	}

      for (int ip = 0; ip < 4; ip++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }
  throw("Too many iterations in routine jacobi.");
}


#endif

