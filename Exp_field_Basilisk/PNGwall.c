/**
# Building a Wall

Building a wall is an imporant theme. Here we build a wall without
sharp edges.
*/
#include "grid/multigrid3D.h"
#include "embed.h"
#include "run.h"
#include "PointTriangle.h"
#include "view.h"

double TEND = 1140;

scalar wall[];
double H = 5, W = 150, D = 6; //Height (y), width (z) and `depth`(x)
double xp = 200, zp = 400;

#define max3(a,b,c) max(max(a,b), c)

// double Point2Rec (const coord *P,
//                   double * min_x,
//                   double * max_x,
//                   double * max_y,
//                   double * min_z,
//                   double * max_z)
// {
//   double d2;
//   if ((*P).x > *min_x && (*P).x < *max_x && (*P).y < *max_y && (*P).z > *min_z && (*P).z < *max_z){
//     d2 = min(min((*P).x-*min_x, *max_x-(*P).x), min((*P).z-*min_z, *max_z-(*P).z));
//   }
//   else if ((*P).x > *min_x && (*P).x < *max_x &&  (*P).y > *max_y && (*P).z > *min_z && (*P).z < *max_z){
//     double dy0 = (*P).y-*max_y;
//     double dx0 = min((*P).x-*min_x, *max_x-(*P).x);
//     double dz0 = min((*P).z-*min_z, *max_z-(*P).z);
//     d2 = sqrt(sq(dy0)+sq(min(dx0,dz0)));
//   }
//   else{
//     double dx = max3(*min_x-(*P).x, 0, (*P).x-*max_x);
//     double dy = max3(-(*P).y, 0, (*P).y-*max_y);
//     double dz = max3(*min_z-(*P).z, 0, (*P).z-*max_z);
//     d2 = sqrt(sq(dx)+sq(dy)+sq(dz));
//   }
//   return d2;
// }


int main() {
  L0 = 800;
  X0 = Y0 = Z0 = 0.;
  N = 64;
  run(); 
}


event init (t = 0) {
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
}
event mov (t += 4) {
    double thst = 0.5, then = 0;
    double phist = 0.2, phien = 0.1;
    double fovst = 15, foven = 10;
    double txst = -0.1, txen = -0.25;
    double tyst = 0, tyen = 0;
    double fov, thetA, phI, tx, ty;
  if (t<300){
    thetA = thst  + t/300*(then  - thst);
    phI   = phist + t/300*(phien - phist);
    fov   = fovst + t/300*(foven - fovst);
    tx    = txst + t/200*(txen - txst);
    ty    = tyst + t/300*(tyen - tyst);
  view (fov = fov, theta = thetA, phi = phI,
        tx = tx, ty = ty, bg = {65./256,157./256,217./256},
        width = 1080, height = 1080);
  }
  if (t>=300 && t <600){
    thetA = thst  + t/300*(then  - thst);
    phI   = phien;
    fov   = foven + (t-300)/300*(-2);
    tx    = txst + t/200*(txen - txst);
    ty    = tyen;
  view (fov = fov, theta = thetA, phi = phI,
        tx = tx, ty = ty, bg = {65./256,157./256,217./256},
        width = 1080, height = 1080);
  }
  if (t>=600){
    thetA = thst  + t/300*(then  - thst);
    phI   = phien;
    fov   = foven + (t-300)/300*(-2);
    tx    = txst + 600/200*(txen - txst);
    ty    = tyen;
  view (fov = fov, theta = thetA, phi = phI,
        tx = tx, ty = ty, bg = {65./256,157./256,217./256},
        width = 1080, height = 1080);
  }
  draw_vof("cs", "fs", edges = true);
  squares ("x + z", n = {0,1,0});
  save ("lam.mp4");
}

event end(t=TEND) {
}