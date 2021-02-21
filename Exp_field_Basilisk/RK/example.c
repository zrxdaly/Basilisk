#include "runge-kutta.h"

int main()
{
   init_grid(1);
   scalar u[];
   double dt = 1e-1;
   foreach()
      u[] = 1.;
   for(t=0; t<=2; t+=dt){
      foreach(){
         printf("%g %g\n", t, u[]);
      }

   }

}
