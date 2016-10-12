reset

set term epslatex color

set linestyle 1 lt 1 lw 1 lc 1
set linestyle 2 lt 1 lw 1 lc 3
set xlabel "Frequenz [Hz]"


set xrange [0:1.5]
set autoscale y
set output "../results/cl_spectra.tex"

plot "../results/spectra.dat" using 1:2 w l ls 1 title "Overset", "../results/spectra.dat" using 1:3 w l ls 2 title "Single"

#pause -1

set xrange [0:0.8]
set autoscale y
set output "../results/cd_spectra.tex"

plot "../results/spectra.dat" using 1:4 w l ls 1 title "Overset", "../results/spectra.dat" using 1:5 w l ls 2 title "Single"

#pause -1

set xlabel "Time [s]"
set xrange [100:200]
set output "../results/cl_timedomain.tex"

plot "../results/overset.dat" using 1:3 w l ls 1 title "Overset", "../results/single.dat" using 1:3 w l ls 2 title "Single"

set output "../results/cd_timedomain.tex"

plot "../results/overset.dat" using 1:2 w l ls 1 title "Overset", "../results/single.dat" using 1:2 w l ls 2 title "Single"

set output "/tmp/foo.tex"
