# set xyplane 0
# set view equal xyz

# set pm3d depthorder border lw 2

# set style fill transparent solid 0.3


# splot fname u 1:2:3 with polygons lc "grey" pointtype 7 fill transparent solid 0.0 notitle

# # splot "data/atoms_test.dat" u 1:2:3:(0.4) with circles lc "black" fill solid 1. noborder fc "black" notitle,\
  
# # splot fname u 1:2:3 with polygons fc "grey" pointtype 7 notitle

set terminal qt size 900,900
# fname = "plots/bonds_test.txt"
fname = "plots/bonds.dat"
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

set size square
set view 62,118

set xyplane 0
set view equal
set pm3d depthorder border lw 2 
set style fill transparent solid 0.6
splot fname u 1:2:3 w lines lc "black" lw 1.5 notitle,\
 fname u 1:2:3 with polygons fc "grey" notitle
#  fname u 1:2:3 w points pointtype 7 lc "black" pointsize 2  notitle,\
