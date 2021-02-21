/**
# Building a Wall

Building a wall is an imporant theme. Here we build a wall without
sharp edges.
*/
#include "grid/multigrid3D.h"
#include "embed.h"
#include "run.h"
// #include "PointTriangle.h"
#include "view.h"

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

int main() {
  L0 = 50;
  X0 = -L0/2;
  N = 256;
  run(); //set cs and fs attributes
}


event init (t = 0) {
  vertex scalar phi[];
  double min_x = xw, max_x = xw + W, min_y = Y0, max_y = Y0 + H, min_z = L0/2-5, max_z = L0/2+5; //L0/2-5
  foreach_vertex() {
    coord cc = {x, y, z};
    phi[] = Point2Rec(&cc, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);
  }
  boundary({phi});
  fractions (phi, cs, fs);
  view (theta = -0.8, phi = 0.3, tx = -0.25); 
  boundary ({cs, fs});
  draw_vof("cs", "fs", edges = true);
  squares ("x + z", n = {0,1,0});
  save ("wall.png");
}
