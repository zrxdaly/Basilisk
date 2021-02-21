#ifndef BASILISK_HEADER_27
#define BASILISK_HEADER_27
#line 1 "/home/dai/software/basilisk/src/run.h"
/**
# Generic time loop

The `run()` function below implements a generic time loop which
executes events until termination.

The timestep `dt` can be accessed as a global variable. */

double dt = 1.;

#include "utils.h"

trace
void run (void)
{
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {

    /**
    We store the total number of cells advanced in time for computing
    speed statistics. */

    update_perf();
    iter = inext, t = tnext;
  }

  /**
  Time/speed statistics are written out on standard output. */

  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
}

#endif