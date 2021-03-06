#include "grid/multigrid.h"
#if MAC
#  include "navier-stokes/mac.h"
#else
#  include "navier-stokes/centered.h"
#endif
#include "view.h"
double END = 120;

int main()
{ 
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  // number of grid points
  init_grid (64);
  // viscosity
#if MAC
  nu = 1e-3;
#else
  const face vector muc[] = {1e-3,1e-3};
  mu = muc;
#endif
  // maximum timestep
  DT = 0.1;
  // CFL number
  CFL = 0.8;
  run();
}

u.t[top] = dirichlet(1);

u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

#if !MAC
uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;
#endif

static double energy()
{
  double se = 0.;
  if (u.x.face)
    foreach(reduction(+:se))
      se += (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.*sq(Delta);
  else // centered
    foreach(reduction(+:se))
      se += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
  return se;
}

event init (i = 0) {
#if !MAC
  restore (file = "lid-restore.dump");
#endif
}

scalar un[];

event logfile (t += 0.1; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-5)
    return 1; /* stop */
  fprintf (stderr, "%f %.9f %g\n", t, energy(), du);
}

event outputfile (i += 100) output_matrix (u.x, stdout, N, linear = true);

#if !MAC
event snapshot (i = 1700)
  dump (file = "dump");
#endif

event profiles (t = end)
{
  FILE * fp = fopen("xprof", "w");
  for (double y = -0.5; y <= 0.5; y += 0.01)
    fprintf (fp, "%g %g\n", y, interpolate (u.x, 0, y));
  fclose (fp);
  
  fp = fopen("yprof", "w");
  for (double x = -0.5; x <= 0.5; x += 0.01)
    fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
  fclose (fp);
}

event mov(t += 1){
    squares("u.x", map = cool_warm);
    save("u.mp4");
}

event end(t=END) {
}