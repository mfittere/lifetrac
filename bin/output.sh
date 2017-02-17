#!/bin/bash -f
# this is a script for processing LIFETRAC output file .ltr
# to produce the following step-by-step data:
#     - normalized beam intensity: 'intensity.txt'
#     - normalized rate of losses: 'lossrate.txt'
#     - specific luminosity (L/L0*N/N0): 'luminosity.txt'
#     - horizontal and vertical emittances: 'emit.txt'
#     - longitudinal beam sigma: 'sigm.txt'
#
if [ $# != 2 ]; then
    if [ $# == 1 ]; then
         echo "Interaction Points for which the luminosity is monitored are:"
         grep -A 3 Luminosity $1 | head -3
         exit 0;
    fi
    echo "This is a script for processing LIFETRAC output"
    echo " "
    echo "invocation:"
    echo "     $0 ltr_file Watch_IP"
    echo " "
    echo "where Watch_IP is the name of the IP to print luminosity, e.g. \"IP_1077\"."
    echo " "
    echo "To find out which IPs are in your ltr file run: "
    echo "$0 ltr_file"
    exit 0;
fi

# path to specific executables
PATH=/home/valishev/lifetrac/bin:/home/shatilov/lifetrac/bin:$PATH
export PATH

grep 'Mean_values' $1 | cut -f2 -d'=' > inti.txt
normint < inti.txt > intensity.txt
lossrate2 < intensity.txt > lossrate.txt
rm inti.txt

grep ": $2" $1 | cut -f3 -d' ' > lumi.txt
speclumi > luminosity.txt
rm lumi.txt

emitgfit $1 > emit.txt
grep '|sigm|' $1 | cut -f10 -d' ' > sigm.txt
