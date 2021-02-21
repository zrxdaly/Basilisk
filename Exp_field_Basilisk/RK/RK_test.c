#include "runge-kutta.h"

static void du (scalar * ul, double t, scalar * kl)
{
  scalar u = ul[0], k = kl[0];
  foreach()
    k[] = t*u[];
}

int main()
{
  init_grid (1);

  for (int order = 1; order <= 4; order *= 2)
    for (double dt = 1e-2; dt <= 8e-2; dt *= 2) {
      scalar u[];
      foreach()
	u[] = 1.; // the initial condition
      double emax = 0.;
      for (t = 0; t <= 2.; t += dt) {
	foreach() {
	  double e = fabs (u[] - exp(t*t/2.));
	  if (e > emax)
	    emax = e;
	  printf ("%g %g %g\n", t, u[], u[] - exp(t*t/2.));
	}
	runge_kutta ({u}, t, dt, du, order);
      }
      printf("\n");
      fprintf (stderr, "%g %g %d\n\n", dt, emax, order);
    }
    printf("\n");
}
