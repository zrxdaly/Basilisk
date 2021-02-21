    
    #include <stdio.h>
    #include <math.h>

int main ()
{
  int i,j,k;
  int N = 256 ;
  double L = 1. ;
  double x,y;
  double dx = L/N ;
  double dt = 0.00001;
  double A[N][N];
  double dA[N][N];

  
// boundary conditions
for (i = 0 ; i < N ; i++) A[i][0] = A[i][N-1] = 0. ;
for (j = 0 ; j < N ; j++) A[0][j] = A[N-1][j] = 0. ;

// initial conditions

  for (i = 0 ; i < N ; i++) 
    {
      for (j = 0 ; j < N ; j++)
        {
          x = i*dx - 0.5 ;
          y = j*dx - 0.5 ;
          A[i][j] = 1./0.1*((fabs(x*x+y*y) < 0.05)) ;
        }
    }

 for (j = 0 ; j < N ; j++)
        {
          printf("%f \n",A[(int)N/2][j]);
        }
 printf("\n\n");

  // time integration

  for (k = 0 ; k < 10 ; k++)
    {

  for (i = 1 ; i < N-1 ; i++) 
    {
       for (j = 1 ; j < N-1 ; j++)
        {
      dA[i][j] = (A[i+1][j] + A[i-1][j] - 2. * A[i][j])/dx/dx +
        (A[i][j+1] + A[i][j-1] - 2. * A[i][j])/dx/dx ;
    }
    }

  // update
  
 for (i = 0 ; i < N ; i++) 
    {
       for (j = 0 ; j < N ; j++)
        {
          A[i][j] =  A[i][j] + dt* dA[i][j] ;
        }
    }

    }

  // print solution (centerline)

       for (j = 0 ; j < N ; j++)
        {
                    printf("%f \n",A[(int)N/2][j]);
        }
    
}