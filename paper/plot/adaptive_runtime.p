set xlabel "max. MSE"
set ylabel "runtime[s]"
set yrange [0:1000]
set xrange [0:1]
unset logscale x
unset logscale y

plot \
"adaptive_runtime.dat" u 1:2 w linespoints t "runtime no adaptive",\
"adaptive_runtime.dat" u 1:3 w linespoints t "runtime adaptive"

# Output
set term png
set output "adaptive_runtime.png"
replot
set term postscript
set output "adaptive_runtime.ps"
replot
set term x11
