unset logscale x
unset logscale y
set key opaque
#set key center top# setze legende
set title ""
set xlabel "Sample point"
set ylabel "MSE"
#set format y "%.0t*10^%+03T"
set ytics nomirror
set yrange [0:0.07]
set xrange [320:480]
#set xrange [0:320]
set xtics 320,40,480
#set xtics 0,40,320
set grid noxtics front
set xtics nomirror

set style fill transparent solid 1 border -1

#plot\
#"expvec_noimportance.dat" w filledcurves  above x1 t "MSE no importance sampling" lt rgb "#2B83BA" ,\
#"expvec_noadaptive.dat" w filledcurves  above x1 t "MSE non adaptive + importance" lt rgb "#D7191C",\
#"expvec_adaptive.dat" w filledcurves above x1 t "MSE adaptive + importance" lt rgb "#ABDDA4",\
#0.005 t "MSE treshold" lt rgb "#ABDDA4" lw 2

plot\
"expvec_noimportance.dat" w boxes   t "MSE no importance sampling" lt rgb "#2B83BA" ,\
"expvec_noadaptive.dat" w boxes   t "MSE non adaptive + importance" lt rgb "#D7191C",\
"expvec_adaptive.dat" w boxes  t "MSE adaptive + importance" lt rgb "#ABDDA4",\
0.005 t "MSE treshold" w lines lt rgb "#ABDDA4" lw 2


# Output
set term pngcairo #size 800,400
set output "mse.png"
replot
set term postscript
set output "mse.ps"
replot
set term x11
