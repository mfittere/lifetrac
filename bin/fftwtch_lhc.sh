#!/bin/bash
if [ $# != 5 ]; then
echo "usage: $0 file nshift nfourier xmin xmax"
exit
fi

ln -s $1 fort.10
echo $2 $3 > fort.77
fftwtch < fort.77 > fft.out
echo "set grid" > .fftwtch.gp
echo "set xrange [$4:$5]" >> .fftwtch.gp
echo "set logscale y" >> .fftwtch.gp
echo "plot \"fft.out\" u (\$2):(\$3) title \"x\" w histeps" >> .fftwtch.gp
echo "pause -1 \"press enter\" " >> .fftwtch.gp
echo "plot \"fft.out\" u (\$2):(\$4) title \"y\" w histeps" >> .fftwtch.gp
echo "pause -1 \"press enter\" " >> .fftwtch.gp
gnuplot .fftwtch.gp
rm fort.77 fort.10


