#include "navier-stokes/centered.h"

// face vector muv[];
#define MAXLEVEL 8

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << MAXLEVEL);
  //  TOLERANCE = 1e-12;
  run();
}


event init (t = 0)
{
  foreach()
    foreach_dimension()
      u.x[] = 1;
  boundary ((scalar *){u});
}

event movie (t += 0.2; t <= 60) {
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, linear = true, file = "vort.mp4");
}

#if TREE
event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif