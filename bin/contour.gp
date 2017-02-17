set contour base
unset surface
set table 'footprint.cnt'
splot 'footprint.out'
unset table
set size square
set xrange [0.66:0.71]
set yrange [0.66:0.71]
plot 'footprint.cnt' w p pt 5
pause -1 "ok"
