/**
# Building a Wall

Building a wall is an imporant theme. Here we build a wall without
sharp edges.
*/
#include "grid/multigrid3D.h"
#include "embed.h"
#include "run.h"
#include "view.h"
#include "fan.h"

double TEND = 5;

int main() {
  L0 = 20;
  X0 = Y0 = Z0 = 0.;
  N = 128;
  // fan.prolongation = fraction_refine;
  rot.phit = 2*M_PI/300;
  rot.ST_phit = 0.;
  run(); 
}


event init (t = 0) {
    rot.fan = true;		// Yes we want a fan
    rot.rotate = false;		// If we want it to rotate 
    rot.start = 0;
    rot.stop = 1020;

    if(rot.fan) {
      init_rotor();
	    rotor_coord();
    }
}

event progress(t+=1) {
    fprintf(stderr, "i=%d t=%g \n", i, t);
}

double th = 0;
event mov (t += 1) {
  view (bg = {0.3, 0.3, 0.9}, theta = th, phi = 0.2);
    // double thst = 0.5, then = 0;
    // double phist = 0.2, phien = 0.1;
    // double fovst = 15, foven = 10;
    // double txst = -0.1, txen = -0.25;
    // double tyst = 0, tyen = 0;
    // double fov, thetA, phI, tx, ty;
  // if (t<300){
  //   thetA = thst  + t/300*(then  - thst);
  //   phI   = phist + t/300*(phien - phist);
  //   fov   = fovst + t/300*(foven - fovst);
  //   tx    = txst + t/200*(txen - txst);
  //   ty    = tyst + t/300*(tyen - tyst);
  // view (fov = fov, theta = thetA, phi = phI,
  //       tx = tx, ty = ty, bg = {65./256,157./256,217./256},
  //       width = 1080, height = 1080);
  // }
  // if (t>=300 && t <600){
  //   thetA = thst  + t/300*(then  - thst);
  //   phI   = phien;
  //   fov   = foven + (t-300)/300*(-2);
  //   tx    = txst + t/200*(txen - txst);
  //   ty    = tyen;
  // view (fov = fov, theta = thetA, phi = phI,
  //       tx = tx, ty = ty, bg = {65./256,157./256,217./256},
  //       width = 1080, height = 1080);
  // }
  // if (t>=600){
  //   thetA = thst  + t/300*(then  - thst);
  //   phI   = phien;
  //   fov   = foven + (t-300)/300*(-2);
  //   tx    = txst + 600/200*(txen - txst);
  //   ty    = tyen;
  // view (fov = fov, theta = thetA, phi = phI,
  //       tx = tx, ty = ty, bg = {65./256,157./256,217./256},
  //       width = 1080, height = 1080);
  // }
  draw_vof("fan", edges = true);
  squares ("x + z", n = {0,1,0});
  save ("lam.mp4");
  th += 0.006;
}

event end(t=TEND) {
}