set terminal qt size 800,600

fname = "plots/pcf.dat"

set ylabel "pcf"
set xlabel "r"

# set yrange [-500:100]
r_sr = 3.582
M = 100
plot fname using ($1*2.5*r_sr/M):2 with boxes lc "black" notitle