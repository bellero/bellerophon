set term epslatex color

set key outside right center

set title "Comparision of residuals for p and U"

set xlabel "Iteration"
set ylabel "Initial residual"

set xtics 1000

set logscale y

set output "../results/residuals.tex"

plot "../results/overset_pResiduals.dat" w l title "p overset", \
     "../results/single_pResiduals.dat" w l title "p single", \
     "../results/overset_UxResiduals.dat" w l title "Ux overset", \
     "../results/single_UxResiduals.dat" w l title "Ux single", \
     "../results/overset_UyResiduals.dat" w l title "Uy overset", \
     "../results/single_UyResiduals.dat" w l title "Uy single"


set output
