unset logscale x
unset logscale y
set key opaque
#set key center top# setze legende
set title ""
set xlabel "Sample points"
set ylabel "MSE"
set yrange [0:0.07]
set xrange [0:3210]
set xtics 0,400,3210

plot\
"expvec_noimportance.dat" w boxes t "MSE no importance sampling" lt rgb "blue",\
"expvec_noadaptive.dat" w boxes t "MSE not adaptive + importance" lt rgb "red",\
"expvec_adaptive.dat" w boxes t "MSE adaptive + importance" lt rgb "green",\
0.0015 t "avg. MSE (not)adaptive"

# Output
set term png
set output "mse.png"
replot
set term postscript
set output "mse.ps"
replot
set term x11
