#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



#include "grid/octree.h" 		// For 3D
//#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion
#include "profile5b.h"


int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum and normal refinement criteria
double TEND = 250.;

// char sim_ID[] = "krab";		        // Simulation identifier
// char sim_var[] = "angle";  		// Notes if a variable is varied over runs

#include "physics.h"			// Physics of the simulation 
//#include "diagnostics.h"

int main() {	
    minlevel = 3;
    maxlevel = 5;

    L0 = 400.;
    X0 = Y0 = Z0 = 0.;

    // Possibility to run for variable changes
    for(double tempVar=265; tempVar<267; tempVar+=100) {
	
	//rot.theta = tempVar*M_PI/180;
	//rot.phit = 2*M_PI/tempVar;                      // Hence a 300 second rotation time

        init_grid(1<<4);
	a = av; 

        foreach_dimension() {
	    u.x.refine = refine_linear;  		// Momentum conserved 
	}

	//fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

  	meps = 10.;					// Maximum adaptivity criterion
	DT = 10E-5;					// For poisson solver 
        TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.8;					// CFL condition

	//sim_dir_create();				// Create relevant dir's
	//out.sim_i++;					// Simulation iteration
 
    	run();						// Start simulation 

    }
}


event init(t=0) {
    //rot.fan = true;		// Yes we want a fan
    //rot.rotate = true;		// If we want it to rotate 
    //rot.start = 30.;            // start time of rotor forcing
    //rot.stop = 1260.;           // stop time of rotor forcing

    //rot.phi = 0;	        // Reset for repeated different runs
    //eps = .5;
    
    init_physics();             // See imported phyisics module for this function

    //if(rot.fan) {               // See imported fan module for these functions
    //    init_rotor();
    //	rotor_coord();
    //}
    // Adapt grid to initial buoyancy and velocity fields.    
    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,0.25*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
	    b[] = STRAT(y);  // physics.h
            u.x[] = WIND(y);
	}
//	rotor_coord();
    }
}


event init_change(i=10) {
    TOLERANCE=10E-3;
    DT = .5;
}

event adapt(i++) {
    adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,.25*9.81/273},maxlevel,minlevel);
}

event progress(t+=5) {
    fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}

//event dumpfields(t=120; t+=120) {
//    char nameDump[90];
//    snprintf(nameDump, 90, "./%s/fielddump", out.dir);
//    dump(file = nameDump, list = all);
//}

event end(t=TEND) {
}

event profiler (t =120; t+=120) {
  char fname[180];
  sprintf (fname, "proft=%dT", (int)(t/300));
  profile ((scalar*){u, b}, fname);
}

