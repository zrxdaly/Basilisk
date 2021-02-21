/**
# Vortex shedding behind a sphere at Reynolds = 300

![Animation of the $\lambda_2$ vortices coloured with the vorticity
 component aligned with the flow.](sphere/movie.mp4)(loop)

We solve the Navier--Stokes equations on an adaptive octree and use
embedded boundaries to define the sphere. */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "PointTriangle.h"
#include "view.h"

/**
We will use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $256^3$. */

double H = 1, W = 4, D = 0.3; //Height (y), width (z) and `depth`(x)
double xp = 10, zp = 0;

int maxlevel = 8;

/**
We need a new field to define the viscosity. */

face vector muv[];

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

int main()
{
  init_grid (128);
  size (20.);
//   X0 = -L0/2;
  X0 = Z0 = 0.;
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

// /**
// The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

// u.n[bottom] = dirichlet(0.);
// u.t[bottom] = dirichlet(0.);
// u.n[top] = dirichlet(0.);
// u.t[top] = dirichlet(1.);

// u.n[left]  = dirichlet(1.);
// p[left]    = dirichlet(0.);
// pf[left]   = dirichlet(0.);

// u.n[right] = neumann(0.);
// p[right]   = neumann(0.);
// pf[right]  = neumann(0.);

event init (t = 0) {
  vertex scalar phi[];
  D /= 2.;
  H -= D/2;
  W -= W/2;
  foreach_vertex() {
    coord cc = {x, y, z};
    if ((z - zp) > W/2) { // distance to a side edge
      coord p1 = {xp, Y0, zp + W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
      coord tmp1[1]; double tmp2[1];
      phi[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else if ((z - zp) < -W/2) {
      coord p1 = {xp, Y0, zp - W/2.}, p2 = {xp, Y0 + H, zp - W/2.};
      coord tmp1[1]; double tmp2[1];
      phi[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else if (y - Y0 > H) {
      coord p1 = {xp, Y0 + H, zp - W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
      coord tmp1[1]; double tmp2[1];
      phi[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else
      phi[] = fabs(cc.x - xp);
  }
  fractions (phi, cs, fs, D);
  boundary ({cs, fs});
  foreach()
    u.x[] = cs[] ? 1. : 0.;
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use Basilisk view to create the animated isosurface of $\lambda_2$
for $30 <= t <= 60$. */

event movies (t = 1; t += 0.25; t <= 30)
{
  #if (dimension == 2)
    scalar omega[];
    vorticity (u, omega);
    boundary ({omega});
    draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.1,0.2});
    draw_vof ("cs", "fs");
    squares ("omega", linear = true, map = cool_warm);
    mirror ({0,-1})
    cells();
  #elif (dimension == 3)
  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */

  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);

  view (fov = 10, theta = -0.8, phi = 0.4, 
	tx = 0, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);

  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  #endif
  save ("movie.mp4");

}

event dumper (t += 10) {
  
  p.nodump = false;
  char str[99];
  sprintf (str, "dump3D%g", t);
  dump(str);
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt (i++) {
  astats s = adapt_wavelet ({cs,u}, (double[]){1e-2,0.02,0.02,0.02},
			    maxlevel, 4);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
