/** Include required libraries */
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h"                // For 3D
#include "view.h"                       // For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"                     // Tracers
#include "diffusion.h"                  // Diffusion 

/** Global variables */
int minlevel, maxlevel;                 // Grid depths
double meps, eps;                       // Maximum error and error in u fields
double TEND = 5;

char sim_ID[] = "krab";                 // Simulation identifier
// char sim_var[] = "Trot";             // Notes if a variable is varied over runs

#include "physics.h"                    // Physics of the simulation 
#include "fan.h"                        // Include a fan
#include "output_vlices.h"
#include "lambda2.h"
// #include "diagnostics.h"             // Perform diagnostics

/** Initialisation */
int main() {
    minlevel = 4;
    maxlevel = 8;

    L0 = 800.;
    X0 = Y0 = Z0 = 0.;

    // Possibility to run for variable changes
    // for(double tempVar=90; tempVar<110; tempVar+=3) {
    // //for(rot.theta=90*M_PI/180; rot.theta<106*M_PI/180.; rot.theta+=3*M_PI/180) {

    //     rot.theta = tempVar*M_PI/180;

    init_grid(1<<6);
    a = av;

    foreach_dimension() {
        u.x.refine = refine_linear;             // Momentum conserved 
    }

        fan.prolongation = fraction_refine;             // Fan is a volume fraction
        p.refine = p.prolongation = refine_linear;
        b.gradient = minmod2;                           // Flux limiter 

    rot.phit = 2*M_PI/240;
    // rot.ST_phit = M_PI;
    rot.ST_phit = 0.; // full rotation

        meps = 10.;                                     // Maximum adaptivity criterion
        DT = 10E-5;                                     // For poisson solver 
    TOLERANCE=10E-6;                            // For poisson solver 
        CFL = 0.8;                                      // CFL condition

//      sim_dir_create();                               // Create relevant dir's
        // out.sim_i++;                                 // Simulation iteration

    run();                                              // Start simulation 

}


/** Initialisation */
event init(t=0) {
    rot.fan = true;             // Yes we want a fan
    rot.rotate = false;          // If we want it to rotate 
    rot.start = 0;
    rot.stop = 1500;

    // rot.phi = 0;             // Reset for different runs
    eps = .5;

    init_physics();

    if(rot.fan) {
        init_rotor();
            rotor_coord();
    }

    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,0.35*9.81/273},maxlevel,minlevel).nf) {
        foreach() {
            b[] = STRAT(y);
        u.x[] = WIND(y);
        }
        rotor_coord();
    }
}

/** Return to standard tolerances and DTs for poisson solver */
event init_change(i=10) {
    TOLERANCE=10E-3;
    DT = .5;
}

/** Adaptivity */
event adapt(i++) {
    adapt_wavelet((scalar *){fan,u,b},(double []){0.1,eps,eps,eps,.35*9.81/273},maxlevel,minlevel);
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
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 30, theta = 0.5, phi = 0.4,
        tx = -0.1, ty = 0.1, bg = {65./256,157./256,217./256},
        width = 1080, height = 1080);
  box();
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


char dir_slices[61];

event slices(t += 1) {
    char B_Slice[91];
    char U_Slice[91];
    coord slice = {1., 0., 1.};
    int res = L0/2;

    for(double yTemp = 1.; yTemp<=10.; yTemp+=1.) {
            slice.y = yTemp/L0;

        snprintf(B_Slice, 90, "%st=%05gy=%03g", "./resultslice/buo/", t, yTemp);
        FILE * fpsli = fopen(B_Slice, "w");
        output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
        fclose(fpsli);
        snprintf(U_Slice, 90, "%st=%05gy=%03g", "./resultslice/vel/", t, yTemp);
        FILE * fpsli_u = fopen(U_Slice, "w");
        output_slice(list = (scalar *){u.x}, fp = fpsli_u, n = res, linear = true, plane=slice);
        fclose(fpsli_u);
    }
}

/** End the simulation */
event end(t=TEND) {
}