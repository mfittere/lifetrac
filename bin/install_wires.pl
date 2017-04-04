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
    if( $buf[0] eq 'bbwire_l5.b1->charge' ){ $l5c=$buf[2]; $i=$i+1;};
    if( $buf[0] eq 'bbwire_r5.b1->charge' ){ $r5c=$buf[2]; $i=$i+1; };
    if( $buf[0] eq 'bbwire_l1.b1->charge' ){ $l1c=$buf[2]; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r1.b1->charge' ){ $r1c=$buf[2]; $i=$i+1; };

    if( $buf[0] eq 'bbwire_l5.b1->xma' ){ $l5x=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r5.b1->xma' ){ $r5x=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_l1.b1->xma' ){ $l1x=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r1.b1->xma' ){ $r1x=-$buf[2]*100; $i=$i+1; };

    if( $buf[0] eq 'bbwire_l5.b1->yma' ){ $l5y=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r5.b1->yma' ){ $r5y=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_l1.b1->yma' ){ $l1y=-$buf[2]*100; $i=$i+1; };
    if( $buf[0] eq 'bbwire_r1.b1->yma' ){ $r1y=-$buf[2]*100; $i=$i+1; };
}
close(fpr1);
printf "done.\n";
#
if($i != 12){ printf "Read wrong number of parameters %d, must be 12, exiting.\n",$i; exit(0);};
#
printf "Installing elements\n";
#
$n=`wc -l $ARGV[1]`;
$n1=`grep -n '# BBWIREL5:' $ARGV[1] |cut -f 1 -d ':'`;
$n2=`grep -n '# BBWIRER5:' $ARGV[1] |cut -f 1 -d ':'`;
$n3=`grep -n '# BBWIREL1:' $ARGV[1] |cut -f 1 -d ':'`;
$n4=`grep -n '# BBWIRER1:' $ARGV[1] |cut -f 1 -d ':'`;

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

$ni=$n3+1;
$no=$n-$n3-2;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Current: $l1c\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");
$ni=$n3+5;
$no=$n-$n3-6;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Shift  (cm):         (x)=$l1x (y)=$l1y (comp)=1\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");

$ni=$n4+1;
$no=$n-$n4-2;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Current: $r1c\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");
$ni=$n4+5;
$no=$n-$n4-6;
system("head -$ni tmp.1 > tmp.2");
system("echo \"Shift  (cm):         (x)=$r1x (y)=$r1y (comp)=1\" >> tmp.2");
system("tail -$no tmp.1 >> tmp.2");
system("cp tmp.2 tmp.1");

system("cp $ARGV[1] $ARGV[1].old");
system("cp tmp.1 $ARGV[1]");
system("rm tmp.1 tmp.2");
printf "done. Check $ARGV[1] and $ARGV[1].old\n";

