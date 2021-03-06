#include "grid/quadtree.h"
#include "view.h"

int main() {
  L0 = 16;
  X0 = Y0 = 0;
  init_grid (4); // Initialize a 2 x 2 grid
  refine ((x > 12) && (y > 12) && (level < 3)); // Refine to top right corner
  unrefine ((x < 8) && (y < 8) && level >= 1); // Coarsen the bottom left corner
  printf ("#All Cells :\n");
  foreach_cell()
    printf ("%d %g %g %g\n", level, x, y, Delta);
  printf ("\n#leafs:\n");
  foreach()
    printf ("%d %g %g %g\n", level, x, y, Delta);
  int i = 1; 
  boundary (all);
  while (i <= depth()){
    fprintf (ferr, "\n#halo ghosts at level %d\n", i);
    foreach_halo (prolongation, i - 1){
      foreach_child(){
        fprintf (ferr, "%d\t%g\t%g\t%g\n", level, x, y, Delta);
      }
    }
    i++;
  }
  view (fov = 40, tx = -0.5, ty = -0.5);
  cells();
  save("grid.png");
}