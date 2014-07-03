#! /usr/bin/env gnuplot
clear
reset
set key on
set border 3
set grid

set xlabel "MSE"
set ylabel "sample points"
set logscale y


set key opaque
set xtics 0,0.0025,0.01
#set ytics 0,100,300
set xrange [0.00025:0.01]
set yrange [1:300]

set ytics nomirror
set xtics nomirror

set grid  front lw 2


set parametric
const=0.005
set trange [0:500]

x = 3

# Output
set term pdfcairo enhanced font "Serif, 14"


set output "../graphics/mse_histogram_uniform.pdf"
plot\
   "mse_noimportance_hist.dat" u 2:3 t "Uniform" s bezier w filledcurves above x1 lw x lc rgb "forest-green"  fs transparent solid 0.5



set output "../graphics/mse_histogram_uniform_is.pdf"
plot\
   "mse_noimportance_hist.dat" u 2:3 t "Uniform" s bezier w filledcurves above x1 lw x lc rgb "forest-green"  fs transparent solid 0.5,\
   "mse_noadaptive_hist.dat" u 2:3 t "IS" s bezier w filledcurves above x1 lw x  lc rgb "blue" fs transparent solid 0.5



set term svg enhanced font "Serif, 14"
set output "../graphics/mse_histogram_is.svg"
plot\
   "mse_noadaptive_hist.dat" u 2:3 t "IS" s bezier w filledcurves above x1 lw x  lc rgb "blue" fs transparent solid 0.5


set output "../graphics/mse_histogram_is_as.svg"
plot\
   "mse_adaptive_hist.dat" u 2:3 t   "IS + AS" s bezier w filledcurves above x1 lw x lc rgb "red"  fs transparent solid 0.5,\
   "mse_noadaptive_hist.dat" u 2:3 t "IS" s bezier w filledcurves above x1 lw x  lc rgb "blue" fs transparent solid 0.5
