/** Include required libraries */ 
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h" 		// For 3D
#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion 

#include "embed.h"
#include "PointTriangle.h"

/** Global variables */
int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields
double TEND = 10;
double H = 4, W = 100, D = 1; //Height (y), width (z) and `depth`(x)
double xp = 200, zp = 400;

#include "physics.h"			// Physics of the simulation 
#include "fan.h"			// Include a fan
#include "lambda2.h"

/** Initialisation */
int main() {	
    minlevel = 4;
    maxlevel = 9;

    L0 = 800.;
    X0 = Y0 = Z0 = 0.;


    init_grid(1<<8);
    a = av; 

    foreach_dimension() {
        u.x.refine = refine_linear;  		// Momentum conserved 
    }

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

    rot.phit = 2*M_PI/160;
    rot.ST_phit = 0.; // full rotation

  meps = 10.;					// Maximum adaptivity criterion
	DT = 0.5;					// For poisson solver 
  TOLERANCE=10E-4;				// For poisson solver 
	CFL = 0.8;					// CFL condition

    run();						// Start simulation 
}

/** Initialisation */
event init(t=0) {
    rot.fan = true;		// Yes we want a fan
    rot.rotate = true;		// If we want it to rotate 
    rot.start = 0;
    rot.stop = 1020;

    eps = .5;
    
    init_physics();

    u.n[embed] = dirichlet(0.);
    u.t[embed] = dirichlet(0.);
    #if dimension == 3
        u.r[embed] = dirichlet(0.);
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

    while(adapt_wavelet((scalar *){cs, u, b},(double []){0.05,eps,eps,eps,0.35*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
	    b[] = STRAT(y);
        u.x[] = cs[]? WIND(y): 0.;
	}
	rotor_coord();
    }
}

/** Return to standard tolerances and DTs for poisson solver */ 
// event init_change(i=10) {
//     TOLERANCE=10E-3;
//     DT = .5;
// }

/** Adaptivity */
event adapt(i++) {
    adapt_wavelet((scalar *){cs, u,b},(double []){0.05, 0.01,eps,eps,eps,.35*9.81/273},maxlevel,minlevel);
}

/** Progress event */
event progress(t+=5) {
    fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}
event mov (t += 0.5) {
#if (dimension == 2)
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  // draw_vof ("cs", "fs", filled = -1, fc = {0,0,0});
  // draw_vof ("cs", "fs");
  squares ("omega", linear = true, map = cool_warm);
  mirror ({0,-1})
  cells();
#elif (dimension == 3)
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
  isosurface ("l2", -0.01);
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
  save ("lam_3D.mp4");
}


// char dir_slices[61];

// event slices(t += 1) {
//     char nameSlice[91];
//     coord slice = {1., 0., 1.};
//     int res = L0/2;

//     for(double yTemp = 0.5; yTemp<=5; yTemp+=0.5) {
// 	    slice.y = yTemp/L0;

//     	snprintf(nameSlice, 90, "%st=%05gy=%03g", dir_slices, t, yTemp*10);
//     	FILE * fpsli = fopen(nameSlice, "w");
//     	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
//     	fclose(fpsli);
//     }
// }


/** End the simulation */
event end(t=TEND) {
}
