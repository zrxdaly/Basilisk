#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion 

#include "PointTriangle.h"
#include "view.h"

#include "fan.h"			// Include a fan
#include "lambda2.h"

#if dimension == 3
	#include "SGS.h"
#endif 

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $256^3$. */

int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields

double TEND = 600;
double H = 3, W = 50, D = 1; //Height (y), width (z) and `depth`(x)
double xp = 75, zp = 150;

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .3 	// Kelvin per meter
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.1     // roughness heat length

#define WIND(s) (max(0.1*log(s/roughY0u),0.))   // log Wind profile 
#define QFLX 0. 	// 0 (0.001 = 20wm2)
#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define GFLX (-Lambda*(BSURF - bd))
double Lambda = 0.005, bd = 0.;   // Grass coupling
#define STRAT(s) gCONST/TREF*s*INVERSION

scalar b[];
scalar * tracers = {b};

// double crho = 1.;	
face vector av[]; 
// double crho = 1.;

/**
We need a new field to define the viscosity. */

// face vector muc[];

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

int main()
{
    minlevel = 4;
    maxlevel = 8;
    L0 = 300.;
    X0 = Y0 = Z0 = 0.;
    // mu = muc;
    a = av; 
    fan.prolongation = fraction_refine;		// Fan is a volume fraction
    rot.phit = 2*M_PI/100;
    rot.ST_phit = 0.; // full rotation
    init_grid(1<<6);

    run();
}

/**
The viscosity is just $1/Re$, because we chose a sphere of diameter
unity and an unit inflow velocity. */

// event properties (i++)
// {
//   foreach_face()
//     muc.x[] = fm.x[]/300.;
// }

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(WIND(y));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
#if dimension == 3
  Evis[bottom] = dirichlet(0.); // Flux is explicitly calculated
  Evis[top] = dirichlet(0.); 
  u.r[embed] = dirichlet(0.);
#endif 

b[bottom] = BSURF;
b[top] = dirichlet(STRAT(y));

// b[embed] = dirichlet(0.);

event init (t = 0) {
  rot.fan = true;		// Yes we want a fan
  rot.rotate = true;		// If we want it to rotate 
  rot.start = 0;
  rot.stop = 1020;

  periodic (left);
  #if dimension == 3
    periodic (front);
  #endif 

  if(rot.fan) {
    init_rotor();
    rotor_coord();
  }

  vertex scalar phw[];
  D /= 2.;
  H -= D/2;
  W -= W/2;
  foreach_vertex() {
    coord cc = {x, y, z};
    if ((z - zp) > W/2) { // distance to a side edge
      coord p1 = {xp, Y0, zp + W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
      coord tmp1[1]; double tmp2[1];
      phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else if ((z - zp) < -W/2) {
      coord p1 = {xp, Y0, zp - W/2.}, p2 = {xp, Y0 + H, zp - W/2.};
      coord tmp1[1]; double tmp2[1];
      phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else if (y - Y0 > H) {
      coord p1 = {xp, Y0 + H, zp - W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
      coord tmp1[1]; double tmp2[1];
      phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
    }
    else
      phw[] = fabs(cc.x - xp);
  }
  fractions (phw, cs, fs, D);
  boundary ({cs, fs});
  foreach(){
    b[] = STRAT(y);
    u.x[] = cs[]? WIND(y) : 0.;
    }
}

event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

mgstats mgb;
/* Diffusion */
event tracer_diffusion(i++){
    scalar r[];
    foreach() {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX)/sq(Delta); // div needed as normalization 
    }
    
    double flx = 0, bt = 0;
    double fctr = CP*TREF/gCONST;
    foreach_boundary(bottom reduction(+:flx) reduction(+:bt)) {
        flx = flx + (QFLX + GFLX) * sq(Delta);
         bt = bt + BSURF * sq(Delta);
    }
    bt = bt/sq(L0);
    flx = flx/sq(L0);
    fprintf(stderr, "soil=%g %g %g %d\n", t, fctr*flx, fctr*bt/CP, i);  
    
    mgb = diffusion(b, dt, mu, r = r);
}
/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);



event movies (t = 1; t += 0.25; t <= TEND)
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

  view (fov = 30, theta = 0.5, phi = 0.4, 
	tx = -0.1, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  box();
  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  translate (y = -5){
    squares ("u.x", n = {0,1,0}, alpha = 5.,
	     min = -1.2, max = 2, map = cool_warm);}
  translate (z = -L0/2.){
  squares ("b", n = {0,0,1}, alpha = L0/2.,
      min = -0.1, max = 0.6, map = cool_warm);}
  translate (x = -L0/2.){
  squares ("b", n = {1,0,0}, alpha = L0/2.,
      min = -0.1, max = 0.6, map = cool_warm);}
  #endif
  save ("movie.mp4");

}

event adapt (i++) {
  astats s = adapt_wavelet ({cs,fan, u, b}, (double[]){1e-2,1e-2,0.05,0.05,0.05, .2*9.81/273},
			    maxlevel, minlevel);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}

event end(t=TEND) {
}