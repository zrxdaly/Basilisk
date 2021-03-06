#define RAD (sqrt(sq(x) + sq(y)))
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "profile6.h"
#include "fan.h"

scalar b[], * tracers = {b};

double R = 0.005, B = -2e-3;
double  nuv = 1.5e-5, kappav = 1.5e-5;
b[embed] = neumann (B/kappav*((y/RAD) + 1)/2.); // azimuthal Flux condition
int maxlevel = 9;
double be = 5e-3, ue = 1e-3;
face vector nu[], av[];
double U;

int main() {
  periodic (left);
  L0 = 20.;
  X0 = Y0 = -L0/2;
  mu = nu;
  a = av;

  rot.phit = 2*M_PI/100;
  rot.ST_phit = 0.; // left 180 rotation

  DT = 0.005;
  run();
  U = 0.1;
  run();
}

event properties (i++) {
  foreach_face()
    nu.x[] = fs.x[]*nuv;
  boundary((scalar*){nu});
}

event init (t = 0) {
  rot.fan = true;		// Yes we want a fan
  rot.rotate = false;		// If we want it to rotate
  rot.start = 0.;            // start time of rotor forcing
  rot.stop = 1020.;   
  if(rot.fan) {
    init_rotor();
    rotor_coord();
  }
  refine (RAD < 2*R && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = RAD - R;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach(){
    u.x[] = U*(cs[] > 0);
    b[] = 0.24;
  }
    
  boundary (all);
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] = (b[0,-1] + b[])/2.;
}

event heating_forcing (i++) {
  double tau = 1.;
  foreach() {
    if (y < Y0 + 3*R || x < X0 + 3*R || x > X0 + L0 - 3*R) {
      b[] -= b[]*dt/tau;
      u.x[] -= (u.x[] - U)*dt/tau;
    }
  }
  boundary ({b, u.x});
}


event adapt (i++) {
  adapt_wavelet ({fan, cs, b, u}, (double[]){0.05, 1e-9, be, ue, ue, ue}, maxlevel);
}

event tracer_diffusion (i++)
  diffusion (b, dt, nu);

event mov (t += 0.05) {
  double bs;
  scalar bd[];
  foreach()
    bd[] = x < 0 ? b[] : nodata;
  interface_average ({b}, &bs, cs, fs);
  squares ("bd", min = bs - 0.01, max = -bs + 0.01, map = cool_warm );
  translate (z = 1e-5)
    draw_vof ("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
  translate (z = -1e-3)
    cells();
  char str[99];
  sprintf (str, "U = %g m/s", U);
  draw_string (str, 15);
  save ("b.mp4");
}


event logger (t += 0.1) {
  double bs;
  interface_average ({b}, &bs, cs, fs);
  printf ("%g %g\n", t, 27.3*bs);
}

event stop (t = 5) {
  char fname[99];
  sprintf (fname, "prof%g", U);
  profile ({b}, RAD, fname);
}