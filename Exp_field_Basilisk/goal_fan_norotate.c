// This code test the goal function of different setting of fan
// #include "grid/octree.h"
// #include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "output_vlices.h"
// #include "lambda2.h"
#include "fan.h"

#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define WIND(s) 0.1*log((s+0.1)/0.1)
#define G 9.81
#define Tref 273.
#define INV 0.3
#define P_B(s) G/Tref*s*INV
double TEND = 240.;

scalar b[], * tracers = {b};

face vector av[];
// const face vector muc[] = {5e-3, 5e-3};
face vector muv[];

int maxlevel = 12;

int main() {
  #if dimension == 2
    periodic (left);
  #elif dimension == 3
    periodic(front);
  #endif
  L0 = 1000;
  X0 = -L0/2;

  rot.phit = 2*M_PI/60;
  // rot.ST_phit = M_PI;
  rot.ST_phit = 0.; // full rotation

  mu = muv;
  a = av;
  N = 128;
  run();
}

// later can be modify to SGS code
event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/300.;
  boundary ((scalar*){muv});
}


u.n[left]  = dirichlet(WIND(y));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);
// u.r[embed] = dirichlet(0.);

// u.n[bottom] = dirichlet(0.);
// u.t[bottom] = dirichlet(0.);
// u.n[top] = dirichlet(0.);
// u.t[top] = dirichlet(WIND(y));

// u.n[embed] = dirichlet (0.);
// u.t[embed] = dirichlet (0.);

// b[bottom] = neumann(0.1*exp(-sq(x))/(mu.y[])); // A local heat source
b[bottom] = BSURF;
b[top] = dirichlet (P_B(y));

// #if dimension == 3
//   u.r[top] = dirichlet(WIND(y));
//   u.r[bottom] = dirichlet(0.); 
//   u.t[top] = neumann(0.); 

//   u.r[embed] = dirichlet (0.);
// #endif  


event init (t = 0) {
  rot.fan = true;		// Yes we want a fan
  rot.rotate = false;		// If we want it to rotate
  rot.start = 60.;            // start time of rotor forcing
  rot.stop = 240.;           // stop time of rotor forcing
  
  if(rot.fan) {
    init_rotor();
    rotor_coord();
  }

  foreach() {
    b[] = P_B(y);
    // u.x[] = cs[] ? WIND(y) : 0.;
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

// event damp (i++) {
//   double relaxtime = dt/50;
//   foreach() {
//     if (y > L0*3/4) 
//       u.x[] = u.x[] + (WIND(y)-u.x[])*relaxtime;
//       b[] = b[] + (P_B(y) - b[])*relaxtime;
//       u.y[] = u.y[] - u.y[]*relaxtime;
//       u.z[] = u.z[] - u.z[]*relaxtime;

//   }
//   boundary ((scalar*){u});
// }

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
  save ("vor_Fix.mp4");
}

event mov (t += 0.5) {
  // draw_vof ("cs", "fs", filled = -1, fc = {0,0,0});
  // draw_vof ("cs", "fs");
  squares ("b", linear = true, min = -0.1, max = 0.5);
  mirror ({0,-1})
  cells();
  save ("buo_Fix.mp4");
}

char dir_slices[60];

event slices (t += 0.5) {
  char nameslice[91];
  int res = L0/2;
  coord slice = {1., 1.};
  snprintf(nameslice, 90, "%st=%05g", dir_slices, t*10.);
  FILE * fpsli = fopen(nameslice, "w");
  fprintf (fpsli, "t=%05g\n", t);
  
  for (double yT = 0.5; yT <= 5.; yT+=0.5){
    slice.y = yT/L0;
    output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane =slice);
  } 
  fclose(fpsli);
}


event adapt (i++) 
  adapt_wavelet ((scalar *){fan, b, u}, (double[]){0.05, 0.05, 0.05, 0.05}, maxlevel, 4);

event end(t=TEND) {
}
