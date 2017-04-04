#--- opening files
if( ($ARGV[0] eq "") || ($ARGV[1] eq "") ){
    print "usage: perl install_wires.pl madx_wire_params lifetracfile\n ";
    exit(0);
}
open(fpr1, $ARGV[0]) || die "Cannot open $ARGV[0] $!\n";
#open(fpr2, $ARGV[1]) || die "Cannot open $ARGV[1] $!\n";

printf "Reading wire parameters file (mad-x output)...\n";
$i=0;
while( <fpr1> ){
    @buf=split ;
    if( $buf[0] eq 'bbwire_l5.b2->charge' ){ $l5c=$buf[2]; $i=$i+1;};
    if( $buf[0] eq 'bbwire_r5.b2->charge' ){ $r5c=$buf[2]; $i=$i+1; };

    if( $buf[0] eq 'bbwire_l5.b2->xma' ){ $l5x=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r5.b2->xma' ){ $r5x=-$buf[2]*100; $i=$i+1; };

    if( $buf[0] eq 'bbwire_l5.b2->yma' ){ $l5y=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r5.b2->yma' ){ $r5y=-$buf[2]*100; $i=$i+1; };
}
close(fpr1);
printf "done.\n";
#
if($i != 6){ printf "Read wrong number of parameters %d, must be 6 (l+r x 3 param), exiting.\n",$i; exit(0);};
#
printf "Installing elements\n";
#
$n=`wc -l $ARGV[1]`;
$n1=`grep -n '# BBWIREL5:' $ARGV[1] |cut -f 1 -d ':'`;
$n2=`grep -n '# BBWIRER5:' $ARGV[1] |cut -f 1 -d ':'`;

system("cp $ARGV[1] tmp.1");

$ni=$n1+1;
$no=$n-$n1-2;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Current: $l5c\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");
$ni=$n1+5;
$no=$n-$n1-6;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Shift  (cm):         (x)=$l5x (y)=$l5y (comp)=1\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");

$ni=$n2+1;
$no=$n-$n2-2;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Current: $r5c\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");
$ni=$n2+5;
$no=$n-$n2-6;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Shift  (cm):         (x)=$r5x (y)=$r5y (comp)=1\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");

system("cp $ARGV[1] $ARGV[1].old");
system("cp tmp.1 $ARGV[1]");
system("rm tmp.1 tmp.2");
printf "done. Check $ARGV[1] and $ARGV[1].old\n";

