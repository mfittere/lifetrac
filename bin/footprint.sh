#!/bin/bash
if [ $# != 3 ]; then
echo "usage: $0 file nturn npart"
exit
fi

ln -s $1 fort.10
echo $2 $3 > fort.77
footprint < fort.77 > footprint.out
gnuplot /home/valishev/lifetrac/bin/contour.gp
rm fort.77 fort.10
