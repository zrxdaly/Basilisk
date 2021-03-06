set terminal png
set output "N_B.png"
set xr [0.005:0.05]

set xlabel 'Radius'
set ylabel 'Temperature'
set key bottom right
set grid
set size square
plot 'prof0' u 1:($2*27.3) w l lw 2 t 'U = 0 m/s', \
     'prof0.1' u 1:($2*27.3) w l lw 2 t 'U = 0.1 m/s', \
     'prof0.2' u 1:($2*27.3) w l lw 2 t 'U = 0.2 m/s', \
     'prof0.3' u 1:($2*27.3) w l lw 2 t 'U = 0.3 m/s'
