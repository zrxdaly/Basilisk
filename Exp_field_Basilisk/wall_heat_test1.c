// #include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "PointTriangle.h" 
// #include "lambda2.h"

#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define WIND(s) 0.1*log((s+0.1)/0.1)
#define G 9.81
#define Tref 273.
#define INV 0.3
#define P_B(s) G/Tref*s*INV
#define max3(a,b,c) max(max(a,b), c)

double Point2Rec (const coord *P,
                  double * min_x,
                  double * max_x,
                  double * min_y,
                  double * max_y,
                  double * min_z,
                  double * max_z)
{
  double d2;
  if ((*P).x > *min_x && (*P).x < *max_x && (*P).y > *min_y && (*P).y < *max_y && (*P).z > *min_z && (*P).z < *max_z){
    d2 = max3(max(*min_x-(*P).x, (*P).x-*max_x), max(*min_y-(*P).y, (*P).y-*max_y), max(*min_z-(*P).z, (*P).z-*max_z));
  }
  else{
    double dx = max3(*min_x-(*P).x, 0, (*P).x-*max_x);
    double dy = max3(*min_y-(*P).y, 0, (*P).y-*max_y);
    double dz = max3(*min_z-(*P).z, 0, (*P).z-*max_z);
    d2 = sqrt(sq(dx)+sq(dy)+sq(dz));
  }
  return d2;
}


double xw = -10, W = 1, H = 10; 

double TEND = 100.;

scalar b[], * tracers = {b};

face vector av[];

int maxlevel = 8;

face vector muv[];

int main() {
  #if dimension == 2
    periodic (left);
  #elif dimension == 3
    periodic(front);
  #endif
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

// u.n[bottom] = dirichlet(0.);
// u.t[bottom] = dirichlet(0.);
// u.n[top] = dirichlet(0.);
// u.t[top] = dirichlet(WIND(y));

// u.n[left]  = dirichlet(WIND(y));
// p[left]    = dirichlet(0.);
// pf[left]   = dirichlet(0.);

// u.n[right] = neumann(0.);
// p[right]   = neumann(0.);
// pf[right]  = neumann(0.);


// u.n[embed] = dirichlet (0.);
// u.t[embed] = dirichlet (0.);
// #if dimension == 3
//   u.r[top] = dirichlet(WIND(y));
//   u.r[bottom] = dirichlet(0.); 
//   // u.t[top] = neumann(0.); 
//   // u.r[left] = dirichlet (0.);   
//   u.r[embed] = dirichlet (0.);        
// #endif  

u.n[left]  = dirichlet (WIND(y));
uf.n[left] = WIND(y);
u.t[left]  = dirichlet (0);
p[left]    = dirichlet (0.);
pf[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
#if (dimension == 3)
u.r[embed] = dirichlet (0.);
u.r[left]  = dirichlet (0);
#endif

// // b[bottom] = BSURF;
#if dimension ==2
  b[bottom] = neumann(0.1*exp(-sq(x))/(5e-3)); // A local heat source
#elif dimension == 3
  b[bottom] = neumann(0.1*exp(-sq(x)-sq(z-L0/2.))/(5e-3)); // A local heat source
#endif  
b[top] = dirichlet (P_B(y));


event init (t = 0) {
  refine (sq(x+9.5)  < sq(1) && y  < 11 && level < maxlevel);
  // refine (sq(x+9.5)  < sq(1) && y  < 11 && sq(z-25)  < sq(11) && level < maxlevel);
  vertex scalar phi[];
  double min_x = xw, max_x = xw + W, min_y = Y0, max_y = Y0 + H, min_z = 0, max_z = 0; //L0/2-5
  foreach_vertex() {
    coord cc = {x, y};
    phi[] = Point2Rec(&cc, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);
  }
  boundary({phi});
  fractions (phi, cs, fs);
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



event mov (t = 1; t += 0.25; t <= 100) {
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
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 15, theta = -0.8, phi = 0.4, 
	tx = -0.25, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  isosurface ("l2", -0.01);
  cells (alpha = -L0/2);
  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
#endif
  save ("mov20003D.mp4");
}


event mov (t += 0.5) {
  draw_vof ("cs", "fs", filled = -1, fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  squares ("b", linear = true, min = -0.1, max = 0.5);
  mirror ({0,-1})
  cells();
  save ("mov.mp4");
}

event adapt (i++) 
  adapt_wavelet ({cs, b, u}, (double[]){0.01, 0.1, 0.05, 0.05}, maxlevel, 4);


event end(t=TEND) {
}
