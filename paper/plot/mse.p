unset logscale x
unset logscale y
set key opaque
#set key center top# setze legende
set title ""
set xlabel "Sample point"
set ylabel "Mean squared error"
set yrange [0:0.07]
set xrange [0:3210]
plot\
"expvec_noimportance.dat" w boxes t "MSE not importance" lt rgb "blue",\
"expvec_noadaptive.dat" w boxes t "MSE not adaptive" lt rgb "red",\
"expvec_adaptive.dat" w boxes t "MSE adaptive" lt rgb "green",\
0.0015 t "avg. MSE (not)adaptive"

# Output
set term png
set output "mse.png"
replot
set term postscript
set output "mse.ps"
replot
set term x11
