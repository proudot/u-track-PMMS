#ifndef BRESENHAM_HPP
# define BRESENHAM_HPP

# include <vector.hpp>

template <typename T, typename F>
void bresenham(const vector<2,T>& p1,
	       const vector<2,T>& p2,
	       F& func)
{
  int l, q, dqr, dqru;
      
  vector<2,T> p, incr1, incr2;
      
  int dy = p2[1] - p1[1];
  int dx = p2[0] - p1[0];
      
  incr1[1] = dy < 0 ? (dy = -dy, -1) : 1;
  incr1[0] = dx < 0 ? (dx = -dx, -1) : 1;
      
  if (dy >= dx)
    {
      dqr = dx << 1;
      dqru = dqr - (dy << 1);
      q = dqr - dy;
      l = dy;
      incr2[1] = incr1[1];
    }
  else
    {
      dqr = dy << 1;
      dqru = dqr - (dx << 1);
      q = dqr - dx;
      l = dx;
      incr2[0] = incr1[0];
    }
      
  p = p1;
      
  for (int d = l; d >= 0; --d)
    {
      func(p);
	  
      if (q > 0)
	{
	  p += incr1;
	  q += dqru;
	}
      else
	{
	  p += incr2;
	  q += dqr;
	}
    }
}

#endif /* BRESENHAM_HPP */
