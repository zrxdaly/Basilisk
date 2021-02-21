#include "grid/cartesian.h"
#include "run.h"

scalar A[];
scalar dA[];

double dt;

int main() {
  L0 = 1.;
  N = 256;
  run();
}

// boundary conditions

A[left] = 0.0 ;
A[top] = 0.0 ;
A[right] = 0.0 ;
A[bottom] = 0.0 ;

// initial conditions
event init (t = 0) {
  foreach()
    A[] =  1./0.1*(fabs(x*x+y*y)<0.05);
  boundary ({A});
}

// integration

event integration (i++) {
  foreach()
    dA[] =  (A[1,0] + A[-1,0] - 2. * A[])/Delta/Delta +
    (A[0,1] + A[0,-1] - 2. * A[])/Delta/Delta ;

 foreach()
    A[] = A[] + dt*dA[];
  boundary ({A});

}

// print
event print (i=10) {

  for (double y = 0 ; y <= L0; y += 0.01){
    printf("%g %g \n", 
      y, interpolate (A, 0, y));
}

}