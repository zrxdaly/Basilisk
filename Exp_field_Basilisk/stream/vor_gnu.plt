set term @PNG enhanced size 640,426
set output 'vorticity.png'
set size ratio -1
unset key
unset xtics
unset ytics
unset border
unset colorbox
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
                          0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
                      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
                      0.875 0.9333 0 0, 1 0.498 0 0 )

set multiplot layout 2,3 scale 1.6,1.6
splot 'omega-0'
splot 'omega-1'
splot 'omega-2'
splot 'omega-3'
splot 'omega-4'
splot 'omega-5'
unset multiplot