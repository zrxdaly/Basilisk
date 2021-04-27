set terminal png
set output "B_N2.png"
set xr [0.005:0.05]

set xlabel 'Radius'
set ylabel 'Temperature'
set key bottom right
set grid
set size square
plot 'prof_BN20' u 1:($2*27.3) w l lw 2 t 'U = 0 m/s', \
     'prof_BN20.1' u 1:($2*27.3) w l lw 2 t 'U = 0.1 m/s', \
     'prof_BN20.2' u 1:($2*27.3) w l lw 2 t 'U = 0.2 m/s', \
     'prof_BN20.3' u 1:($2*27.3) w l lw 2 t 'U = 0.3 m/s'