// This code test the goal function of different setting of fan
// #include "grid/octree.h"
// #include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "output_vlices.h"
#include "fan.h"

#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define WIND(s) 0.1*log((s+0.1)/0.1)
#define G 9.81
#define Tref 273.
#define INV 0.3
#define P_B(s) G/Tref*s*INV

double TEND = 1560.;

scalar b[], * tracers = {b};

face vector av[];
// const face vector muc[] = {5e-3, 5e-3};
face vector muv[];

int maxlevel = 9;

int main() {
  #if dimension == 2
    periodic (left);
  #elif dimension == 3
    periodic(front);
  #endif
  L0 = 800;
  X0 = -L0/2;

  rot.phit = 2*M_PI/140;
  // rot.ST_phit = M_PI;
  rot.ST_phit = -90.; // full rotation

  mu = muv;
  a = av;
  N = 128;
  run();
}

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/300.;
  boundary ((scalar*){muv});
}


u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(WIND(y));

b[bottom] = BSURF;
b[top] = dirichlet (P_B(y));

event init (t = 0) {
  rot.fan = true;		// Yes we want a fan
  rot.rotate = true;		// If we want it to rotate
  rot.start = 60.;            // start time of rotor forcing
  rot.stop = 1500.;           // stop time of rotor forcing
  
  if(rot.fan) {
    init_rotor();
    rotor_coord();
  }

  foreach() {
    b[] = P_B(y);
    u.x[] = WIND(y);
  }
  boundary (all);
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++) 
  diffusion (b, dt, mu);

event inflow(i++){
    double sides = 50;
    double relaxtime = dt/50;
    foreach(){
	if((x < sides || x > L0-sides) ||
	  //  (z < sides || z > L0-sides) ||
	   (y > L0-2*sides )) {
	    double a = (x < sides) ? x : fabs(x-L0);
	    a = 1.; 
	    u.x[] = u.x[] + a*(WIND(y)-u.x[])*relaxtime;
 	    b[] = b[] + a*(P_B(y) - b[])*relaxtime;
	    u.y[] = u.y[] - a*u.y[]*relaxtime;
		// u.z[] = u.z[] - a*u.z[]*relaxtime;
	}
    }
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
  view (fov = 15, theta = -0.8, phi = 0.4, 
	tx = -0.25, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  isosurface ("l2", -0.01);
  cells (alpha = -L0/2);
  // draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
#endif
  save ("vor_2D.mp4");
}

event mov (t += 0.5) {
  squares ("b", linear = true, min = -0.1, max = 0.5);
  mirror ({0,-1})
  cells();
  save ("mov.mp4");
}

char dir_slices[200];

event slices (t += 1) {
  char nameslice[91];
  int res = L0/2;
  coord slice = {1., 1.};
  snprintf(nameslice, 90, "%st=%05g", dir_slices, t);
  FILE * fpsli = fopen(nameslice, "w");
  fprintf (fpsli, "t=%05g\n", t);
  
  for (double yT = 0.5; yT <= 5.; yT+=0.5){
    slice.y = yT/L0;
    output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane =slice);
  } 
  fclose(fpsli);
}


event adapt (i++) 
  adapt_wavelet ((scalar *){fan, b, u}, (double[]){0.05, 0.05, 0.05, 0.05}, maxlevel, 5);

event end(t=TEND) {
}
