#ifndef CHISQ_HPP
#define CHISQ_HPP

namespace conan {

  namespace detail {

#ifdef GOODNESS_OF_FIT
    double poz(
        double z
        )
    {
      double const Z_MAX = 6.0; /* Maximum meaningful z value */
      double y, x, w;

      if (z == 0.0)
        x = 0.0;
      else
      {
        y = 0.5 * abs(z);
        if (y >= (Z_MAX * 0.5))
          x = 1.0;
        else if (y < 1.0)
        {
          w = y * y;
          x = ((((((((0.000124818987 * w
                   - 0.001075204047) * w + 0.005198775019) * w
                   - 0.019198292004) * w + 0.059054035642) * w
                   - 0.151968751364) * w + 0.319152932694) * w
                   - 0.531923007300) * w + 0.797884560593) * y * 2.0;
        }
        else
        {
          y -= 2.0;
          x = (((((((((((((-0.000045255659 * y
                         + 0.000152529290) * y - 0.000019538132) * y
                         - 0.000676904986) * y + 0.001390604284) * y
                         - 0.000794620820) * y - 0.002034254874) * y
                         + 0.006549791214) * y - 0.010557625006) * y
                         + 0.011630447319) * y - 0.009279453341) * y
                         + 0.005353579108) * y - 0.002141268741) * y
                         + 0.000535310849) * y + 0.999936657524;
        }
      }
      return z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5);
    }

    double ex(
        double x
        )
    {
      double const BIGX = 20.0; /* max value to represent exp(x) */
      return (x < -BIGX) ? 0.0 : exp(x);
    }

    double pochisq(
        double x,
        int df)
    {
      double const LOG_SQRT_PI = 0.5723649429247000870717135, /* log(sqrt(pi)) */
                   I_SQRT_PI = 0.5641895835477562869480795,   /* 1 / sqrt(pi) */
                   BIGX = 20.0; /* max value to represent exp(x) */
      double a, y, s, e, c, z;
      bool even;

      if (x <= 0.0 || df < 1)
        return 1.0;

      a = 0.5 * x;
      even = df % 2; // = !(df & 1);
      if (df > 1)
        y = ex(-a);

      s = (even ? y : (2.0 * poz(-sqrt(x))));
      if (df > 2)
      {
        x = 0.5 * (df - 1.0);
        z = (even ? 1.0 : 0.5);
        if (a > BIGX)
        {
          e = (even ? 0.0 : LOG_SQRT_PI);
          c = log(a);
          while (z <= x)
          {
            e += log(z);
            s += ex(c * z - a - e);
            z += 1.0;
          }
          return s;
        }
        else
        {
          e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
          c = 0.0;
          while (z <= x)
          {
            e = e * (a / z);
            c += e;
            z += 1.0;
          }
          return c * y + s;
        }
      }
      else
      {
        return s;
      }
    }

    double critchi(
        double alpha,
        int df
        )
    {
      // constants
      double const CHI_EPSILON = 0.000001,
                   CHI_MAX = 99999.0;
      // variables
      double minchisq = 0.0,
             maxchisq = CHI_MAX,
             chisqval = 0.0;

      if (alpha <= 0.0)
        return maxchisq;
      else if (alpha >= 1.0)
        return 0.0;

      chisqval = df / sqrt(alpha);
      while ((maxchisq - minchisq) > CHI_EPSILON)
      {
        if (pochisq(chisqval, df) < alpha)
          maxchisq = chisqval;
        else
          minchisq = chisqval;
        chisqval = (maxchisq + minchisq) / 2;
      }
      return chisqval;
    }
#endif // GOODNESS_OF_FIT

  } // detail

} // conan

#endif // CHISQ_HPP
