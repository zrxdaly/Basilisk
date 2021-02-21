#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "PointTriangle.h" 
#include "lambda2.h"

#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define WIND(s) 0.1*log((s+0.1)/0.1)
#define G 9.81
#define Tref 273.
#define INV 0.3
#define P_B(s) G/Tref*s*INV

// #define WALL (fabs(x + 10) - 0.5 - (y<6))
// #define WALL ((fabs(x + 10)) - 0.5 - (y>10))

// #if dimension == 3
// #define WALL (y < 10 && fabs(x + 10) < 0.5 && fabs(z - 25) < 10)
// #endif  
// scalar wall[];

double xw = -10, W = 0.5, H = 10; 

double TEND = 50.;

scalar b[], * tracers = {b};

face vector av[];

int maxlevel = 6;

face vector muv[];

int main() {
  // periodic (left);
  // #if dimension == 3
  periodic(front);
  // #endif
  L0 = 50;
  X0 = -L0/2;
  // mu = muc;
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

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

// b[bottom] = neumann(0.1*exp(-sq(x))/(5e-3)); // A local heat source
// // b[bottom] = BSURF;
// #if dimension == 3
b[bottom] = neumann(0.1*exp(-sq(x)-sq(z-L0/2.))/(5e-3)); // A local heat source
// #endif  

b[top] = dirichlet (P_B(y));

#if dimension == 3
  u.r[top] = dirichlet(WIND(y));
  u.r[bottom] = dirichlet(0.); 
  u.t[top] = neumann(0.); 

  u.r[embed] = dirichlet (0.);        
#endif  

// event init (t = 0) {
//   // refine (sq(x+10)  < sq(1.5) && y  < 6 && level < maxlevel);
//   // refine (WALL<0.1 && WALL>-0.1 && level < maxlevel);
//   // #if dimension == 3
//     // refine (sq(x+10)  < sq(1) && y  < 11 && sq(z-25)  < sq(11) && level < maxlevel);
//   // #endif  
//   vertex scalar phi[];
//   foreach_vertex()
//     phi[] = WALL;
//   boundary ({phi});
//   fractions (phi, cs, fs);

//   foreach() {
//     b[] = P_B(y);
//     // u.x[] = cs[] ? WIND(y) : 0.;
//     u.x[] = WIND(y);
//   }
// }

event init (t = 0) {
  vertex scalar phih[];
  coord wall1 = {xw, L0/2-5., Y0}, wall2 = {xw, L0/2+5., Y0 + H - W}, wall3 = {xw, L0/2+5., Y0};
  double temp2[1];
  double temp1[1];
  foreach_vertex() {
    coord cc = {x, z, y};
    phih[] = sqrt(PointTriangleDistance (&cc, &wall1, &wall2, &wall3,
				       temp1, temp2));
  }
  fractions (phih, cs, fs, W);
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
//       u.y[] = u.y[] - u.y[]*relaxtime;
//       // u.z[] = u.z[] - u.z[]*relaxtime;
//       b[] = b[] + (P_B(y) - b[])*relaxtime;

//   // if (WALL)
// 	// foreach_dimension()
// 	//   u.x[] *= exp(-dt*10);
//   // }
//   boundary ((scalar*){u});
// }




event mov (t += 0.5) {
#if (dimension == 2)
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  // scalar Wall[];
  // foreach()
  //   Wall[] = WALL? 0 : nodata;
  // squares ("Wall", map = gray, min = 0, max = 1);
  squares ("omega", linear = true, map = cool_warm);
  mirror ({0,-1})
  cells();
#elif (dimension == 3)
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 25, quat = {0.431384,-0.216693,-0.317091,0.816338},
	tx = -5, ty = 0, bg = {0.3,0.4,0.6}, width = 500,
	height = 500, samples = 3);
  isosurface ("l2", -0.01);
  cells (alpha = -L0/2);
  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
#endif
  // char str[99];
  // sprintf (str, "Re = %g, C = %g, ML= %d", Re, c, maxlevel);
  // draw_string (str, 1, lw = 3, lc = {1, 0, 1});
  save ("mov20003D.mp4");
}



event mov (t += 0.5) {
  draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  squares ("b", linear = true, min = -0.1, max = 0.5);
  // mirror ({0,-1});
  // cells();
  save ("mov.mp4");
}

event adapt (i++) 
  adapt_wavelet ({cs, b, u}, (double[]){0.01, 0.1, 0.05, 0.05, 0.05}, maxlevel, 3);


event end(t=TEND) {
}
