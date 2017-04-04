#!/bin/bash
echo "File: tev" > tmp.ltl
echo "Export: (hist)=>gauss" >> tmp.ltl
/home/shatilov/lifetrac/src/lifetrac tmp.ltl > tmp.hist.all

ns=`grep "Step:" tmp.hist.all | wc -l`
echo "Number of steps: " $ns

echo "f(x)=c*exp(-(x-a)**2/2/b**2)">tmp.gnuplot.x
echo "fit f(x) 'tmp.hist' u 1:2 via 'tmp.fitvals'">>tmp.gnuplot.x
echo "update 'tmp.fitvals'">>tmp.gnuplot.x
echo "f(x)=c*exp(-(x-a)**2/2/b**2)">tmp.gnuplot.y
echo "fit f(x) 'tmp.hist' u 1:5 via 'tmp.fitvals'">>tmp.gnuplot.y
echo "update 'tmp.fitvals'">>tmp.gnuplot.y
echo "f(x)=c*exp(-(x-a)**2/2/b**2)">tmp.gnuplot.z
echo "fit f(x) 'tmp.hist' u 1:8 via 'tmp.fitvals'">>tmp.gnuplot.z
echo "update 'tmp.fitvals'">>tmp.gnuplot.z

echo "# GsigX sigX GsigY sigY GsigZ sigZ" > gfit.out

i=1
while (($i<=$ns)) ; do
echo $i
if (($i<10)) ; then
sx=`grep -A 15 "Step:     $i" tmp.hist.all | grep "|sigm|"| cut -f 2 -d' '`
sy=`grep -A 15 "Step:     $i" tmp.hist.all | grep "|sigm|"| cut -f 6 -d' '`
sz=`grep -A 15 "Step:     $i" tmp.hist.all | grep "|sigm|"| cut -f 10 -d' '`
else if (($i<100)) ; then
sx=`grep -A 15 "Step:    $i" tmp.hist.all | grep "|sigm|"| cut -f 2 -d' '`
sy=`grep -A 15 "Step:    $i" tmp.hist.all | grep "|sigm|"| cut -f 6 -d' '`
sz=`grep -A 15 "Step:    $i" tmp.hist.all | grep "|sigm|"| cut -f 10 -d' '`
else if (($i<1000)) ; then
sx=`grep -A 15 "Step:   $i" tmp.hist.all | grep "|sigm|"| cut -f 2 -d' '`
sy=`grep -A 15 "Step:   $i" tmp.hist.all | grep "|sigm|"| cut -f 6 -d' '`
sz=`grep -A 15 "Step:   $i" tmp.hist.all | grep "|sigm|"| cut -f 10 -d' '`
fi
fi
fi

grep "S_$i " tmp.hist.all > tmp.hist
echo "a = 1.0" > tmp.fitvals
echo "b = 1.0" >>tmp.fitvals
echo "c = 1.0" >>tmp.fitvals
gnuplot tmp.gnuplot.x 2> tmp.gout
x=`grep b tmp.fitvals | cut -f2 -d'='`
echo "a = 1.0" > tmp.fitvals
echo "b = 1.0" >>tmp.fitvals
echo "c = 1.0" >>tmp.fitvals
gnuplot tmp.gnuplot.y 2> tmp.gout
y=`grep b tmp.fitvals | cut -f2 -d'='`
echo "a = 1.0" > tmp.fitvals
echo "b = 1.0" >>tmp.fitvals
echo "c = 1.0" >>tmp.fitvals
gnuplot tmp.gnuplot.z 2> tmp.gout
z=`grep b tmp.fitvals | cut -f2 -d'='`
echo $x $sx $y $sy $z $sz>> gfit.out
i=$(($i+1))
done

rm tmp.*
