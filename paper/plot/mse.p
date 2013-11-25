#! /usr/bin/env gnuplot
clear
reset
set term wxt enhanced font "Serif, 14"
set key on
set border 3
set grid

set xlabel "MSE"
set logscale y


set key opaque
set xtics 0,0.0025,0.01
set ytics 0,100,300
set xrange [0.00025:0.01]
set yrange [1:300]

set ytics nomirror
set xtics nomirror

set grid  front lw 2


set parametric
const=0.005
set trange [0:500]

x = 3
#fs transparent pattern 4 bo
   
plot\
   "mse_noimportance_hist.dat" u 2:3 t "no importance sampling" s bezier w filledcurves above x1 lw x lc rgb "forest-green"  fs transparent solid 0.5,\
   "mse_noadaptive_hist.dat" u 2:3 t "importance sampling" s bezier w filledcurves above x1 lw x  lc rgb "blue" fs transparent solid 0.5




# Output
set term pngcairo enhanced font "Serif, 14"
set output "mse_importance.png"
replot
set term postscript
set output "mse_importance.ps"
replot
set term wxt enhanced font "Serif, 14"

plot\
   "mse_adaptive_hist.dat" u 2:3 t   "importance + adaptive sampling" s bezier w filledcurves above x1 lw x lc rgb "red"  fs transparent solid 0.5,\
   "mse_noadaptive_hist.dat" u 2:3 t "importance sampling" s bezier w filledcurves above x1 lw x  lc rgb "blue" fs transparent solid 0.5

# Output
set term pngcairo enhanced font "Serif, 14"
set output "mse_adaptive.png"
replot
set term postscript
set output "mse_adaptive.ps"
replot
set term wxt
