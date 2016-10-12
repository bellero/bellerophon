set term epslatex color

set key off

set title "Comparision of computation time"

set style line 1 lc rgb "black"
set style fill pattern 2
set boxwidth 0.5

set ylabel "Simulation time"
set yrange [0:*]

set output "../results/times.tex"

plot "../results/times.dat" using 1:3:xtic(2) with boxes ls 1

set output
