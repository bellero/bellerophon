reset

set term epslatex

set linestyle 1 lw 1 lc 1
set linestyle 2 lw 1 lc 2
set linestyle 3 lw 1 lc 1 pt 2
set linestyle 4 lw 1 lc 2 pt 2

set xlabel "Frequency [\\(Hz\\)]"
set xrange [0:1.5]

set ylabel "CL"
set yrange [2e-4:0.3]
#set yrange [1e-6:0.3]

#set y2range [1e-6:0.02]
set y2range [1e-6:0.1]
set y2tics
set y2label "CD"
set output "../results/cl_spectra.tex"

set logscale y
set logscale y2

plot "../results/donorCellTet_N20_spectra.dat" using 1:3 w l ls 3 title "Overset CL", \
     "../results/single_N20_spectra.dat" using 1:3 w l ls 4 title "Single CL", \
     "../results/donorCellTet_N20_spectra.dat" using 1:2 axis x1y2 w lp ls 1 pi 5 title "Overset CD", \
     "../results/single_N20_spectra.dat" using 1:2 axis x1y2 w lp ls 2 pi 5 title "Single CD"

set output "../results/foo.tex"

set xrange [0:1]
unset logscale y

a(x)=1
b(x)=2
c(x)=3
d(x)=4
e(x)=5
f(x)=6
g(x)=7
h(x)=8
i(x)=9
j(x)=10
k(x)=11
l(x)=12

plot a(x) w lp lc 1 pt 1 pi 5, \
     b(x) w lp lc 2 pt 2 pi 5, \
     c(x) w lp lc 3 pt 3 pi 5, \
     d(x) w lp lc 4 pt 4 pi 5, \
     e(x) w lp lc 5 pt 5 pi 5, \
     f(x) w lp lc 6 pt 6 pi 5, \
     g(x) w lp lc 7 pt 7 pi 5, \
     h(x) w lp lc 8 pt 8 pi 5, \
     i(x) w lp lc 9 pt 9 pi 5, \
     j(x) w lp lc 10 pt 10 pi 5, \
     k(x) w lp lc 11 pt 11 pi 5, \
     l(x) w lp lc 12 pt 12 pi 5

# #pause -1
# 
# set xrange [0:0.8]
# set autoscale y
# set output "../results/cd_spectra.tex"
# 
# plot "../results/spectra.dat" using 1:4 w l ls 1 title "Overset", "../results/spectra.dat" using 1:5 w l ls 2 title "Single"
# 
# #pause -1
# 
# set xlabel "Time [s]"
# set xrange [100:200]
# set output "../results/cl_timedomain.tex"
# 
# plot "../results/overset.dat" using 1:3 w l ls 1 title "Overset", "../results/single.dat" using 1:3 w l ls 2 title "Single"
# 
# set output "../results/cd_timedomain.tex"
# 
# plot "../results/overset.dat" using 1:2 w l ls 1 title "Overset", "../results/single.dat" using 1:2 w l ls 2 title "Single"

set output "/tmp/foo.tex"
