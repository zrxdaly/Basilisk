// #include "grid/cartesian1D.h"
#include "run.h"

scalar T[], dT[], q[];
double dt;
#define EPS 0.1 

int main(){
    L0 = 20;
    X0 = -L0/2.;
    N = 512;
    DT = (L0/N)*(L0/N)/2 ;
    run();
    // foreach_cell()
    //     fprintf(ferr,"x = %g y = %g\n",x,y); 
}

T[left] = dirichlet(3.);
// T[right] = dirichlet(1);
// q[left] = dirichlet(100);
// q[right] = dirichlet(100);

event init(t = 0){
    foreach(){
        T[] = x>=0? 3: 0;
        // T[] = 1./EPS*(fabs(x)<EPS)/2;
    }
    boundary ({T});
}

event printdata (t = 0; t <= 2000 * DT; t += 250 * DT) {
  foreach()
    fprintf (stdout, "%g %g %g %g %g\n", x, T[], q[], dT[], t);
  fprintf (stdout, "\n\n");
}



// event integration (i++) {
//   double dt = DT;
//   dt = dtnext (dt);
//   foreach()
//     q[]= - (T[0,0] - T[-1,0])/Delta;
//   boundary ({q});
//   foreach()
//     dT[] = - (q[1,0]  - q[0,0])/Delta;
//   foreach()
//     T[] += dt*dT[];
//   boundary ({T});
// }

event integration (i++) {
  double dt = DT;
  dt = dtnext (dt);
  foreach()
    dT[] =  (T[1,0] - 2*T[0,0] + T[-1,0])/(Delta * Delta);
  foreach()
    T[] += dt*dT[];
  boundary ({T});
}

event adapt (i++) {
  adapt_wavelet ({T}, (double []){4e-1}, minlevel = 6, maxlevel = 8);
}

