set term epslatex color

# Set various features of the contour plot
set pm3d
unset surface  # don't need surfaces
set view map
set key off
set cntrparam cubicspline  # smooth out the lines
set cntrparam levels incremental -2,0.2,2   # sets the num of contour lines
set pm3d interpolate 20,20 # interpolate the color

#splot 'daten.dat' using 1:2:3
#pause -1

set contour

# Export contour lines for pressure at overset grid
set table '../results/cont_overset_p.dat'
splot '../results/overset_p.dat' using 1:2:4
unset table

# Export contour lines for pressure at single grid
set table '../results/cont_single_p.dat'
splot '../results/single_p.dat' using 1:2:4
unset table

set cntrparam levels incremental 0,0.1,1   # sets the num of contour lines

# Export contour lines for pressure at overset grid
set table '../results/cont_overset_U.dat'
splot '../results/overset_U.dat' using 1:2:(($4**2+$5**2)**0.5)
unset table

# Export contour lines for pressure at overset grid
set table '../results/cont_single_U.dat'
splot '../results/single_U.dat' using 1:2:(($4**2+$5**2)**0.5)
unset table


# Unset the 3d mapping
unset pm3d

# Set a nice color palette
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")

set xrange [0:0.1]
set yrange [0:0.1]

set xlabel "x"
set ylabel "y"

# Plot velocity field as color plot
set output "../results/overset_U.tex"
set title "Velocity magnitude for overset grid"
plot '../results/overset_U.dat' using 1:2:(($4**2+$5**2)**0.5) w imag
#pause -1

set output "../results/single_U.tex"
set title "Velocity magnitude for single grid"
plot '../results/single_U.dat' using 1:2:(($4**2+$5**2)**0.5) w imag
#pause -1

set key outside
set key bottom center
# Plot contours of pressure and velocity
set output "../results/isobars.tex"
set title "Isobars"
l '<./labels.sh ../results/cont_overset_p.dat 0 15 0'
plot '<./labels.sh ../results/cont_overset_p.dat 1 15 0' w l lt -1 title "Overset grid", '<./labels.sh ../results/cont_single_p.dat 1 15 0' w l lt 1 title "Single grid"
#pause -1

unset label
set output "../results/isotachs.tex"
set title "Isotachs"
l '<./labels.sh ../results/cont_overset_U.dat 0 15 0'
plot '<./labels.sh ../results/cont_overset_U.dat 1 15 0' w l lt -1 title "Overset grid", '<./labels.sh ../results/cont_single_U.dat 1 15 0' w l lt 1 title "Single grid"
#pause -1

set output