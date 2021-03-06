#include "grid/bitree.h"
// #include "view.h"

#define FUN(x) (j1(x))
scalar near[]; 
int maxlevel = 19, baselevel = 3; 

int main() {
  L0 = 20;
  X0 = 1e-9;
  init_grid (1 << baselevel);
  // Iteratively refine the grid
  for (int l = baselevel; l < maxlevel; l++) {
    foreach_level (l) {
      near[] = 0; 
      if (FUN(x - Delta/2.)*FUN(x + Delta/2.) <= 0)
	      near[] = (double)true;
    }
    refine (near[] == (double)true && level <= l);
  }
  // Output the zeros:
  foreach_level (depth())
    if (FUN(x - Delta/2.)*FUN(x + Delta/2.) <= 0)
      printf ("%.8g\n", x);
  printf ("Error margin: +/- %.2g\n", L0/(1 << depth()));
  // Count the used cells (including parents)
  int n = 0;
  foreach_cell()
    n++;
  printf ("Number of Cells: %d\n"
	  "Equdistant-grid equivalent number: %d\n", n, 1 << depth());
}