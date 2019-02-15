#ifndef		FAST_MARCHING_HPP
# define	FAST_MARCHING_HPP

# include <map>
# include <vector>
# include <limits>

# include <image.hpp>
# include <window.hpp>

// Maybe put this elsewhere:
template <int B, int N>
struct Pow
{
  enum { ret = B * Pow<B, N - 1>::ret };
};

template <int B>
struct Pow<B, 0>
{
  enum { ret = 1 };
};

template <int n>
class fast_marching
{
private:
  static const unsigned char ACCEPTED	= 0;
  static const unsigned char FAR	= 1;
  static const unsigned char TRIAL	= 2;

  typedef std::multimap<double, vector<n, int> > trial_set_t;

  typedef window<n, (n << 1), int>		win1_t;
  typedef window<n, Pow<3,n>::ret - 1, int>	win2_t;

  static const win1_t win1;
  static const win2_t win2;

public:
  fast_marching(const int size[n]) : u_(size), r_(size), m_(size)
  {
  }

  const image<n, double> & u() const { return u_; }
  const image<n, int> & r() const { return r_; }

  template <typename G>
  void compute(const image<n, double> & f, G & g)
  {
    typename boost::graph_traits<G>::edge_descriptor e;

    //
    // INITIALIZATION
    //
    
    u_.fill(std::numeric_limits<double>::max());
    r_.fill(-1);
    m_.fill(FAR);
    
    u_.border_replicate(u_.margin());
    r_.border_replicate(u_.margin());
    m_.border_replicate(u_.margin());

    trial_.clear(); // safe. should be already empty.

    // Initialization of actions, ancestors and status maps for
    // initial Points.

    for (unsigned k = 0; k < num_vertices(g); ++k)
      {
	const vector<n,int>& p = g[k].p;

	u_[p] = 0.0;
	r_[p] = k;
	m_[p] = ACCEPTED;
      }

    // Insert neightbors of initial Points into the trial set.

    for (unsigned k = 0; k < num_vertices(g); ++k)
      {
	const vector<n,int>& p = g[k].p;
	
	for (unsigned i = 0; i < win1.size(); ++i)
	  {
	    const vector<n,int> q = p + win1.point(i);
	    
	    if (u_.contains(q))
	      {
		switch (m_[q])
		  {
		  case FAR:
		    {
		      upwind_update_(f, q);
		      trial_.insert(std::make_pair(u_[q], q));
		      m_[q] = TRIAL;
		      break;
		    }
		  case ACCEPTED:
		    {
		      bool is_here = false;

		      const int r1 = r_[p];
		      const int r2 = r_[q];

		      tie(e, is_here) = edge(r1, r2, g);

		      if (is_here == false &&
			  degree(r1, g) < g[r1].max_degree &&
			  degree(r2, g) < g[r2].max_degree)
			{
			  // Create a link between p and q
			  e = add_edge(r1, r2, g).first;

			  g[e].sp_r1 = p;
			  g[e].sp_r2 = q;
			}
		    }
		  }
	      }
	  }     
      }

    //
    // MAIN LOOP
    //

    while (!trial_.empty())
      {
	// Move the Point with the lowest time arrival from trial to
	// alive set.
	
	const vector<n,int> p = trial_.begin()->second;
	trial_.erase(trial_.begin());
	m_[p] = ACCEPTED;

	// Here, we check if there is a saddle point between front
	// r_[p] and another front r_[q]. A saddle point represent the
	// location where 2 fronts meet for the first time,
	// i.e. u_[sp] from p and q is minimal. sp is computed after
	// the propagation is finished. There can be multiple saddle
	// point whether p hit multiple fronts at the same time.

	int r1 = r_[p];

	for (unsigned i = 0; i < win2.size(); ++i)
	  {
	    const vector<n,int> q = p + win2.point(i);

	    int r2 = r_[q];
	    
	    if (u_.contains(q) &&
		m_[q] == ACCEPTED &&
		r2 != r1 &&
		edge(r1, r2, g).second == false &&
		degree(r1, g) < g[r1].max_degree &&
		degree(r2, g) < g[r2].max_degree)
	      {
		// There is a saddle point between r1 and r2
		// fronts. sp is between p and q.

		e = add_edge(r1, r2, g).first;
		g[e].sp_r1 = p;
		g[e].sp_r2 = q;
	      }
	  }
	
	// continue with the classic fast marching

	for (unsigned i = 0; i < win1.size(); ++i)
	  {
	    const vector<n,int> q = p + win1.point(i);	    

	    if (u_.contains(q) && m_[q] == FAR)
	      {
		upwind_update_(f, q);
		trial_.insert(std::make_pair(u_[q], q));
		m_[q] = TRIAL;
	      }
	  }
      }

    //
    // CALCULATE SADDLE POINTS
    //

    typename boost::graph_traits<G>::edge_iterator ei, ei_end;

    // TODO: this needs to be replaced by an estimation of the
    // intersection of the 2 fronts.
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
      {
	g[*ei].sp = (g[*ei].sp_r1 + g[*ei].sp_r2) / 2.0;
	g[*ei].time = (*const_cast<const image<n, double> *>(&u_))[g[*ei].sp];
      }
  }

private:
  void upwind_update_(const image<n, double>& f, const vector<n,int>& p);

private:
  image<n, double> u_;
  image<n, int> r_;
  image<n, unsigned char> m_;

  trial_set_t trial_;
};

template <>
const fast_marching<2>::win1_t fast_marching<2>::win1 = neighb_c4<int>();

template <>
const fast_marching<2>::win2_t fast_marching<2>::win2 = neighb_c8<int>();

template <>
const fast_marching<3>::win1_t fast_marching<3>::win1 = neighb_c6<int>();

template <>
const fast_marching<3>::win2_t fast_marching<3>::win2 = neighb_c26<int>();

template <>
void fast_marching<2>::upwind_update_(const image<2,double> & f,
				      const vector<2,int> & p)
{
  assert(fabs(f[p]) >= std::numeric_limits<double>::epsilon());

  double pot = 1.0 / f[p];
    
  vector<2,int> a(p);
  vector<2,int> b(p);

  if (u_(p[0] - 1, p[1]) < u_(p[0] + 1, p[1]))
    a[0]--;
  else
    a[0]++;

  if (u_(p[0], p[1] - 1) < u_(p[0], p[1] + 1))
    b[1]--;
  else
    b[1]++;

  double ua = u_[a];
  double ub = u_[b];

  if (ua > ub)
    {
      std::swap(ua, ub);
      std::swap(a, b);
    }
    
  assert(ua != std::numeric_limits<double>::max() ||
	 ub != std::numeric_limits<double>::max());

  r_[p] = r_[a];

  u_[p] = (pot > (ub - ua)) ? 0.5 *
    (ub + ua + sqrt(2 * pot * pot - (ub - ua) * (ub - ua))) :
    ua + pot;
}

template <>
void fast_marching<3>::upwind_update_(const image<3, double> & f,
				      const vector<3,int> & p)
{
  assert(fabs(f[p]) >= std::numeric_limits<double>::epsilon());

  double pot = 1.0 / f[p];

  vector<3,int> a(p);
  vector<3,int> b(p);
  vector<3,int> c(p);

  if (u_(p[0] - 1, p[1], p[2]) < u_(p[0] + 1, p[1], p[2]))
    a[0]--;
  else
    a[0]++;

  if (u_(p[0], p[1] - 1, p[2]) < u_(p[0], p[1] + 1, p[2]))
    b[1]--;
  else
    b[1]++;

  if (u_(p[0], p[1], p[2] - 1) < u_(p[0], p[1], p[2] + 1))
    c[2]--;
  else
    c[2]++;

  double ua = u_[a];
  double ub = u_[b];
  double uc = u_[c];
  
  if (ua > ub)
    {
      std::swap(ua, ub);
      std::swap(a, b);
    }

  if (ua > uc)
    {
      std::swap(ua, uc);
      std::swap(a, c);
    }

  if (ub > uc)
    {
      std::swap(ub, uc);
      std::swap(b, c);
    }

  assert(ua != std::numeric_limits<double>::max() ||
	 ub != std::numeric_limits<double>::max() ||
	 uc != std::numeric_limits<double>::max());

  r_[p] = r_[a];
  
  double beta = -2 * (ua + ub + uc);
  double delta = beta * beta - 12 * (ua * ua + ub * ub + uc * uc - pot * pot);

  if (delta >= 0)
    {
      double u = (-beta + sqrt(delta)) / 6.0;

      if (u > uc)
	u_[p] = u;
      else
	u_[p] = (pot > (ub - ua)) ? 0.5 *
	  (ub + ua + sqrt(2 * pot * pot - (ub - ua) * (ub - ua))) :
	  ua + pot;
    }
  else
    u_[p] = (pot > (ub - ua)) ? 0.5 *
      (ub + ua + sqrt(2 * pot * pot - (ub - ua) * (ub - ua))) :
      ua + pot;
}

// template <typename InsertIterator>
// void back_propagate(const image<double> & u,
// 		    const image<int> & r,
// 		    const vector<2,int>& p1,
// 		    const vector<2,int>& p2,
// 		    InsertIterator it)
// {
//   static const window<8,int> w8 = neighb_c8<int>();

//   vector<2,int> q, q_min, p(p1);

//   int rpk = r[p2];

//   while (!(p == p2))
//     {
//       it = p;

//       q_min = p;

//       for (unsigned i = 0; i < w8.size(); ++i)
// 	{
// 	  q = p + w8.point(i);

// 	  if (u[q] <= u[q_min] && r[q] == rpk)
// 	    q_min = q;
// 	}

//       assert(!(q_min == p));
	
//       p = q_min;
//     }

//   it = p;
// }

// template <typename InsertIterator>
// void back_propagate_sp(const image<double> & u,
// 		       const image<int> & r,
// 		       const vector<2,int>& p1,
// 		       const vector<2,int>& p2,
// 		       InsertIterator it)
// {
//   static const double eps = 1e-5;
	
//   vector<2,double> q, p(p1);

//   double dx, dy;

//   do
//     {
//       it = p;
//       q = p;
//       dx = .5 * (u(p[0] + 1, p[1]) - u(p[0] - 1, p[1]));
//       dy = .5 * (u(p[0], p[1] + 1) - u(p[0], p[1] - 1));
//       p[0] -= dx;
//       p[1] -= dy;
//     }
//   while (dist(p, p2) > 1 && dist(p, q) > eps);

//   it = p2;
// }

#endif	/* !FAST_MARCHING_HPP */
