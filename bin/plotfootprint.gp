set size square
set xrange [0.66:0.71]
set yrange [0.66:0.71]
set trange [0.66:0.71]
#
n1 = 7.0
l1 = 10.0
m1 = 0.0
a1 = n1/l1
#
n2 = 7.0
l2 = 8.0
m2 = 2.0
a2 = n2/m2
b2 = -1*l2/m2
#
n3 = 7.0
l3 = 6.0
m3 = 4.0
a3 = n3/m3
b3 = -1*l3/m3
#
n4 = 7.0
l4 = 4.0
m4 = 6.0
a4 = n4/m4
b4 = -1*l4/m4
#
n5 = 7.0
l5 = 2.0
m5 = 8.0
a5 = n5/m5
b5 = -1*l5/m5
#
n6 = 7.0
l6 = 0.0
m6 = 10.0
a6 = n6/m6
b6 = -1*l6/m6
#
n7 = 8.0
l7 = 12.0
m7 = 0.0
a7 = n7/l7
#
n8 = 8.0
l8 = 10.0
m8 = 2.0
a8 = n8/m8
b8 = -1*l8/m8
#
n9 = 8.0
l9 = 8.0
m9 = 4.0
a9 = n9/m9
b9 = -1*l9/m9
#
n10 = 8.0
l10 = 6.0
m10 = 6.0
a10 = n10/m10
b10 = -1*l10/m10
#
n11 = 8.0
l11 = 4.0
m11 = 8.0
a11 = n11/m11
b11 = -1*l11/m11
#
n12 = 8.0
l12 = 2.0
m12 = 10.0
a12 = n12/m12
b12 = -1*l12/m12
#
n13 = 8.0
l13 = 0.0
m13 = 12.0
a13 = n13/m13
b13 = -1*l13/m13
#
set xtics 0.01
set ytics 0.01
set parametric
plot 'footprint.cnt' notitle w p pt 5,\
     'wp.dat' notitle w p pt 7 lc rgb "blue" ps 3,\
      a1,t notitle w l lt 1 lc rgb "green",\
      t,a2+b2*t notitle w l lt 1 lc rgb "green",\
      t,a3+b3*t notitle w l lt 1 lc rgb "green",\
      t,a4+b4*t notitle w l lt 1 lc rgb "green",\
      t,a5+b5*t notitle w l lt 1 lc rgb "green",\
      t,a6+b6*t notitle w l lt 1 lc rgb "green",\
      a7,t notitle w l lt 1 lc rgb "green",\
      t,a8+b8*t notitle w l lt 1 lc rgb "green",\
      t,a9+b9*t notitle w l lt 1 lc rgb "green",\
      t,a10+b10*t notitle w l lt 1 lc rgb "green",\
      t,a11+b11*t notitle w l lt 1 lc rgb "green",\
      t,a12+b12*t notitle w l lt 1 lc rgb "green",\
      t,a13+b13*t notitle w l lt 1 lc rgb "green"
pause -1 "ok"
set term png font times 14
set output "footprint.png"
replot
set term postscript eps color enhanced "Times-Roman" 21
set output "footprint.eps"
replot
