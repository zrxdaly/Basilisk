#include "tracer.h"
#include "diffusion.h"
#include "vreman.h"

face vector Km[],Kh[];
(const) face vector Pr;
scalar Evis[]; //Cell Centered diffusivity
double molvis;
double Csmag;
scalar * tracers;

static inline void Evisprol(Point point,scalar s){
  foreach_child()
    Evis[]=bilinear(point,Evis)/4.;
}
static inline void Evisres(Point point,scalar s){
  double sum = 0.;
  foreach_child()
    sum += s[];
  s[] = sum/2.;
}

event defaults(i=0){
  if (dimension!=3) //Allow to run, but give a warning
    fprintf(stdout,"Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n",dimension);
  mu=Kh;
  Pr=unityf;
  molvis=0.;
  Csmag=0.12;
  Evis.prolongation=Evisprol;
  Evis.restriction=Evisres;

  #if TREE
  Evis.refine=no_restriction;
  Evis.coarsen=no_restriction;
  foreach_dimension(){
    Kh.x.refine=no_restriction;
    Km.x.coarsen=no_restriction;
  }
#endif
}

event Eddyvis(i++){
  eddyviscosity(Csmag,u,molvis,Evis);
  boundary({Evis});
  foreach_face(){
    Km.x[]=(Evis[]+Evis[-1])/2;
    Kh.x[]=(Pr.x[]*(Km.x[]-molvis))+molvis;
  }
  boundary_flux({Km,Kh}); 
}

event tracer_diffusion(i++){
  for (scalar s in tracers)
    diffusion(s,dt,Kh);
}
