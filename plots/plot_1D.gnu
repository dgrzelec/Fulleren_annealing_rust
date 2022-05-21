set terminal qt size 800,600

# fname = "plots/energy_tab.dat"
# # fname = "plots/r_tab.dat"
# set ylabel "E"
# set xlabel "iteration x 100"

# set yrange [-500:100]

# plot fname using 1:2 with lines lc "black" notitle


fname = "plots/EN_tab"

set ylabel "E/N"
set xlabel "N"


plot fname using ($1+30):2 with lines lc "black" notitle
