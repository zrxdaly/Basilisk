/**
# Vortex shedding behind a sphere at Reynolds = 300

![Animation of the $\lambda_2$ vortices coloured with the vorticity
 component aligned with the flow.](sphere/movie.mp4)(loop)

We solve the Navier--Stokes equations on an adaptive octree and use
embedded boundaries to define the sphere. */
#include <math.h>
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "view.h"

#include "tracer-particles.h"
#include "scatter2.h"

/**
We will use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $256^3$. */

int maxlevel = 8;
Particles P;
double R1 = 3, H = 15, Hg = 1;
/**
We need a new field to define the viscosity. */

face vector muv[];

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

int main()
{
  init_grid (64);
  size (32.);
  origin (0, 0, 0);
  mu = muv;
  run();
}


/**
The viscosity is just $1/Re$, because we chose a sphere of diameter
unity and an unit inflow velocity. */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/300.;
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

event init (t = 0) {

  /**
  We initially refine only in a sphere, slightly larger than the solid
  sphere. */

  refine (sq(x-16) + sq(y-16) + sq(z-16) < sq(2) && level < maxlevel);

  /**
  We define the unit sphere. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x-16) + sq(y-16) + sq(z-16) - sq(2);
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  We set the initially horizontal velocity to unity everywhere
  (outside the sphere). */
  
  foreach(){
    u.x[] = cs[] ? 1. : 0.;}

  P = init_tp_square (n = 10, ym = H + Hg, l = 2*R1, zp = L0/2);
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use Basilisk view to create the animated isosurface of $\lambda_2$
for $30 <= t <= 60$. */
double th = 0;
event movies (t = 0; t += 0.25; t <= 30)
{

  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */
  
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);

  view (fov = 30, theta = th, phi = 0.2,
	tx = -0.307321, ty = -0.22653, bg = {0.3, 0.3, 0.9},
	width = 802, height = 634);
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  scatter (P);
  box();
  save ("movie.mp4");
  th += 0.006;
}

event movies2 (t = 0; t += 0.25; t <= 30)
{ 
  view (fov = 40, theta = 0, phi = 0,
	tx = 0, ty = 0, bg = {0.3, 0.3, 0.9},
	width = 600, height = 600);
  squares ("u.x", n = {1,0,0}, alpha = 16, min = -1.1, max = 1.1, map = cool_warm);
  draw_vof("cs", fc = {98./256,78./256,44./256});
  // scatter (P);
  save ("movie_u.mp4");
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt (i++) {
  astats s = adapt_wavelet ({cs,u}, (double[]){1e-2,0.02,0.02,0.02},
			    maxlevel, 4);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
