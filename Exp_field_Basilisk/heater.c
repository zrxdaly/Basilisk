#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "lambda2.h"
#include "output_vlices.h"

#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define WIND(s) 0.1*log((s+0.1)/0.1)
#define G 9.81
#define Tref 273.
#define INV 0.3
#define P_B(s) G/Tref*s*INV

double TEND = 600.;
double inter = 5.;
double surface_T = -0.1;

scalar b[], * tracers = {b};

face vector av[];

int maxlevel = 8;

face vector muv[];

int main() {
  periodic (left);
  #if dimension == 3
    periodic(front);
  #endif
  L0 = 125;
  X0 = Y0 = Z0 = 0.;
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

// b[bottom] = neumann(10 * sin(x/5 * 2 * M_PI + M_PI/2)+10); // A local heat source
// b[bottom] = neumann(5 * (sin(x/inter*M_PI*2 + M_PI/2) + abs(sin(z/inter*M_PI*2 + M_PI/2)) + sin (z/inter*M_PI*2 + M_PI/2) + abs(sin(x/inter*M_PI*2 + M_PI/2))) - surface_T); // A local heat source

// b[bottom] = neumann(x>75? (10 * (sin (z/inter*M_PI*2 + M_PI/2) + abs(sin(x/inter*M_PI*2 + M_PI/2))) - surface_T):1); // A local heat source
// b[bottom] = neumann((5 * (sin(x/inter*M_PI*2 + M_PI/2) + abs(sin(z/inter*M_PI*2 + M_PI/2)) + sin(z/inter*M_PI*2 + M_PI/2) + abs(sin(x/inter*M_PI*2 + M_PI/2))))<10.1? (surface_T):(5 * (sin(x/inter*M_PI*2 + M_PI/2) + abs(sin(z/inter*M_PI*2 + M_PI/2)) + sin (z/inter*M_PI*2 + M_PI/2) + abs(sin(x/inter*M_PI*2 + M_PI/2))))); // A local heat source
b[bottom] = neumann((15 * (sin(x/inter*M_PI*2 + M_PI/2) + sin(z/inter*M_PI*2 + M_PI/2)) - 10)<0? surface_T:(15 * (sin(x/inter*M_PI*2 + M_PI/2) + sin(z/inter*M_PI*2 + M_PI/2)) - 10));
// 5 * (np.sin(X/5 * 2 * np.pi) - np.sin(Y/5 * 2 * np.pi)) + 10

// // b[bottom] = BSURF;


b[top] = dirichlet (P_B(y));

#if dimension == 3
  u.r[top] = dirichlet(WIND(y));
  u.r[bottom] = dirichlet(0.); 
  u.t[top] = neumann(0.); 
      
#endif  

event init (t = 0) {

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

mgstats mgb;
event tracer_diffusion (i++) 
  mgb = diffusion (b, dt, mu);


event inflow(i++){
    double sides = 12.5;
    double relaxtime = dt/30;
    foreach(){
	if((x < sides || x > L0-sides) ||
	   (z < sides || z > L0-sides) ||
	   (y > L0-2*sides )) {
	    double a = (x < sides) ? x : fabs(x-L0);
	    a = 1.; 
	    u.x[] = u.x[] + a*(WIND(y)-u.x[])*relaxtime;
 	    b[] = b[] + a*(P_B(y) - b[])*relaxtime;
	    u.y[] = u.y[] - a*u.y[]*relaxtime;
		  u.z[] = u.z[] - a*u.z[]*relaxtime;
	}
    }
}

event mov_b1 (t += 0.5) {
  view (fov = 23, theta = 0, phi = 0,
        tx = -0.50, ty = 0);
  squares ("b", n = {0,0,1}, alpha = inter *13., linear = true, min = -0.1, max = 0.5);
  mirror ({0,-1})
  cells(n = {0,0,1}, alpha = L0/2.);
  save ("mov_b1.mp4");
}

event mov_b2 (t += 0.5) {
  view (fov = 23, theta = -0.2, phi = 0.4,
        tx = -0.6, ty = 0.1);
  squares ("b", n = {0,1,0}, alpha = 3., linear = true, min = -0.1, max = 0.5);
  translate (z = 0.){
    cells(n = {0,0,1}, alpha = inter *13.);}
  save ("mov_b2.mp4");
}

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
    coord slice = {1., 0., 1.};
    int res = L0;

    for(double yTemp = 1.; yTemp<=5.; yTemp+=1.) {
            slice.y = yTemp/L0;

        snprintf(B_Slice, 90, "%st=%05gy=%03g", "./buo/", t, yTemp);
        FILE * fpsli = fopen(B_Slice, "w");
        output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
        fclose(fpsli);
    }
}

event adapt (i++) 
  adapt_wavelet ({b, u}, (double[]){0.1, 0.05, 0.05, 0.05}, maxlevel, 5);


event end(t=TEND) {
}
