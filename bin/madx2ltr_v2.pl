# This script converts MADX twiss output to Lifetrac lattice format
# Version for element-by-element tracking with thin mult - drift approach
# Version with proper CM units for lifetrac
# A.Valishev (valishev@fnal.gov), 5/29/2012
#
# Usage:
#        perl madx2ltr.pl madx.lattice madx.errors strong.optics lifetracfile
#
# madx.lattice file is supposed to have the default full twiss table output format
#
# Data are printed at 1) drifts; 2) thin mults; 3) rf cavities; 4) beam-beam markers
# :::: 2/14/2013 :::: 
# added new types of elements: 5) solenoids; 6) dipedge; 7) (h,v)kicker
#
#use,sequence=lhcb1;
#select,flag=twiss,clear;
#select,flag=twiss,class=drift;
#select,flag=twiss,class=beambeam,pattern=BB_PAR;
#select,flag=twiss,class=beambeam,pattern=BB_HO.*0;
#select,flag=twiss,class=multipole;
#select,flag=twiss,class=rfcavity;
#select,flag=twiss,class=solenoid;
#select,flag=twiss,class=dipedge;
#select,flag=twiss,class=kicker;
#select,flag=twiss,class=hkicker;
#select,flag=twiss,class=vkicker;
#twiss,file="out.lattice";
#
#use,sequence=lhcb2;
#select,flag=twiss,clear;
#select,flag=twiss,class=beambeam,pattern=BB_PAR;
#select,flag=twiss,class=beambeam,pattern=BB_HO.*0;
#twiss,file="out.strong";
#
#use,sequence=lhcb1;
#select,flag=error,clear;
#select,flag=error,class=multipole;
#esave;
#
  
# madx.errors contains multipole errors at the mults
# standard format: NAME K0L K0SL K1L K1SL K2L K2SL ...
# up to the 20th order
# !!!!!!!!!!!!!!!!!!! for now the 'dipole' error is not included !!!!!!!!!!!!!!!!!!
#
#
# strong optics is the madx twiss data for the strong beam at b-b interactions
# default full twiss table output
#
#
use Switch;
use List::Util qw[min max];
#
# Main IP marker for this lattice:
#$mainIP='IP1';
$mainIP='BBHO10';
# Secondary IP, if any
$mainIP2='BBHO50';
#
# Beam parameters:
# Particle type:
$particle='P';
#
#--- opening files
if( ($ARGV[0] eq "") || ($ARGV[1] eq "") || ($ARGV[2] eq "") || ($ARGV[3] eq "")){
    print "usage: perl madx2ltr.pl lattice errors strong.optics lifetracfile\n ";
    exit(0);
}
open(fpr1, $ARGV[0]) || die "Cannot open $ARGV[0] $!\n";
open(fpr2, $ARGV[1]) || die "Cannot open $ARGV[1] $!\n";
open(fpr3, $ARGV[2]) || die "Cannot open $ARGV[2] $!\n";
open(fpw, ">".$ARGV[3]) || die "Cannot open $ARGV[3] $!\n";

#------------------------------------------------------------------------------
printf "Reading machine lattice file...\n";
$n=0;
$fswitch=0;
$nip=0; $ndrift=0; $nmult=0; $nrf=0; $nsol=0; $ndpdg=0; $nkick=0;
while( <fpr1> ){
    @buf=split ;
    if( ($buf[0] eq '@') && ($buf[1] eq 'ENERGY') ){ $energy=$buf[3]; }
    if( ($buf[0] eq '@') && ($buf[1] eq 'GAMMA') ){ $gamma=$buf[3]; 
    $beta=sqrt(1-1/$gamma/$gamma); printf "gamma=%f beta=%lG\n",$gamma,$beta;}
    if( ($buf[0] eq '@') && ($buf[1] eq 'NPART') ){ $Np=$buf[3]; } 
    if( ($buf[0] eq '@') && ($buf[1] eq 'SIGE') ){ $sige=$buf[3]; }
    if( ($buf[0] eq '@') && ($buf[1] eq 'SIGT') ){ $blength=$buf[3]; } 
    if( ($buf[0] eq '@') && ($buf[1] eq 'EX') ){ $emitx=$buf[3]; }
    if( ($buf[0] eq '@') && ($buf[1] eq 'EY') ){ $emity=$buf[3]; }
    if( $buf[0] eq '*' ){
	printf "Input format:\n";
	$fswitch=1;
	$bufsize=scalar(@buf);
	for($j=1;$j<$bufsize;$j++){
	    switch ($buf[$j]) {
		case "NAME"    { $iname = $j-1; printf "NAME    at %d\n",$j; }
		case "KEYWORD" { $ikey  = $j-1; printf "KEYWORD at %d\n",$j; }
		case "S"       { $is    = $j-1; printf "S       at %d\n",$j; }
		case "X"       { $ix    = $j-1; printf "X       at %d\n",$j; }
		case "PX"      { $ipx   = $j-1; printf "PX      at %d\n",$j; }
		case "BETX"    { $ibetx = $j-1; printf "BETX    at %d\n",$j; }
		case "ALFX"    { $ialfx = $j-1; printf "ALFX    at %d\n",$j; }
		case "MUX"     { $imux  = $j-1; printf "MUX     at %d\n",$j; }
		case "DX"      { $idx   = $j-1; printf "DX      at %d\n",$j; }
		case "DPX"     { $idpx  = $j-1; printf "DPX     at %d\n",$j; }
		case "Y"       { $iy    = $j-1; printf "Y       at %d\n",$j; }
		case "PY"      { $ipy   = $j-1; printf "PY      at %d\n",$j; }
		case "BETY"    { $ibety = $j-1; printf "BETY    at %d\n",$j; }
		case "ALFY"    { $ialfy = $j-1; printf "ALFY    at %d\n",$j; }
		case "MUY"     { $imuy  = $j-1; printf "MUY     at %d\n",$j; }
		case "DY"      { $idy   = $j-1; printf "DY      at %d\n",$j; }
		case "DPY"     { $idpy  = $j-1; printf "DPY     at %d\n",$j; }
		case "L"       { $il    = $j-1; printf "L       at %d\n",$j; }
		case "LRAD"    { $ilrad = $j-1; printf "LRAD    at %d\n",$j; }
		case "KSI"     { $iksi  = $j-1; printf "KSI     at %d\n",$j; }
		case "H1"      { $ih1   = $j-1; printf "H1      at %d\n",$j; }
		case "E1"      { $ie1   = $j-1; printf "E1      at %d\n",$j; }
		case "FINT"    { $ifint = $j-1; printf "FINT    at %d\n",$j; }
		case "HGAP"    { $ihgap = $j-1; printf "HGAP    at %d\n",$j; }
		case "TILT"    { $itilt = $j-1; printf "TILT    at %d\n",$j; }
		case "VOLT"    { $ivolt = $j-1; printf "VOLT    at %d\n",$j; }
		case "LAG"     { $ilag  = $j-1; printf "LAG     at %d\n",$j; }
		case "FREQ"    { $ifreq = $j-1; printf "FREQ    at %d\n",$j; }
		case "HKICK"   { $ihkick= $j-1; printf "HKICK   at %d\n",$j; }
		case "VKICK"   { $ivkick= $j-1; printf "VKICK   at %d\n",$j; }
		case "K0L"     { $ikl[0] = $j-1; printf "K0L     at %d\n",$j; }
		case "K1L"     { $ikl[1] = $j-1; printf "K1L     at %d\n",$j; }
		case "K2L"     { $ikl[2] = $j-1; printf "K2L     at %d\n",$j; }
		case "K3L"     { $ikl[3] = $j-1; printf "K3L     at %d\n",$j; }
		case "K4L"     { $ikl[4] = $j-1; printf "K4L     at %d\n",$j; }
		case "K5L"     { $ikl[5] = $j-1; printf "K5L     at %d\n",$j; }
		case "K6L"     { $ikl[6] = $j-1; printf "K6L     at %d\n",$j; }
		case "K7L"     { $ikl[7] = $j-1; printf "K7L     at %d\n",$j; }
		case "K8L"     { $ikl[8] = $j-1; printf "K8L     at %d\n",$j; }
		case "K9L"     { $ikl[9] = $j-1; printf "K9L     at %d\n",$j; }
		case "K10L"    { $ikl[10]= $j-1; printf "K10L    at %d\n",$j; }
		case "K11L"    { $ikl[11]= $j-1; printf "K11L    at %d\n",$j; }
		case "K12L"    { $ikl[12]= $j-1; printf "K12L    at %d\n",$j; }
		case "K13L"    { $ikl[13]= $j-1; printf "K13L    at %d\n",$j; }
		case "K14L"    { $ikl[14]= $j-1; printf "K14L    at %d\n",$j; }
		case "K15L"    { $ikl[15]= $j-1; printf "K15L    at %d\n",$j; }
		case "K16L"    { $ikl[16]= $j-1; printf "K16L    at %d\n",$j; }
		case "K17L"    { $ikl[17]= $j-1; printf "K17L    at %d\n",$j; }
		case "K18L"    { $ikl[18]= $j-1; printf "K18L    at %d\n",$j; }
		case "K19L"    { $ikl[19]= $j-1; printf "K19L    at %d\n",$j; }
		case "K20L"    { $ikl[20]= $j-1; printf "K20L    at %d\n",$j; }
		case "K0SL"    { $iks[0] = $j-1; printf "K0SL    at %d\n",$j; }
		case "K1SL"    { $iks[1] = $j-1; printf "K1SL    at %d\n",$j; }
		case "K2SL"    { $iks[2] = $j-1; printf "K2SL    at %d\n",$j; }
		case "K3SL"    { $iks[3] = $j-1; printf "K3SL    at %d\n",$j; }
		case "K4SL"    { $iks[4] = $j-1; printf "K4SL    at %d\n",$j; }
		case "K5SL"    { $iks[5] = $j-1; printf "K5SL    at %d\n",$j; }
		case "K6SL"    { $iks[6] = $j-1; printf "K6SL    at %d\n",$j; }
		case "K7SL"    { $iks[7] = $j-1; printf "K7SL    at %d\n",$j; }
		case "K8SL"    { $iks[8] = $j-1; printf "K8SL    at %d\n",$j; }
		case "K9SL"    { $iks[9] = $j-1; printf "K9SL    at %d\n",$j; }
		case "K10SL"   { $iks[10]= $j-1; printf "K10SL   at %d\n",$j; }
		case "K11SL"   { $iks[11]= $j-1; printf "K11SL   at %d\n",$j; }
		case "K12SL"   { $iks[12]= $j-1; printf "K12SL   at %d\n",$j; }
		case "K13SL"   { $iks[13]= $j-1; printf "K13SL   at %d\n",$j; }
		case "K14SL"   { $iks[14]= $j-1; printf "K14SL   at %d\n",$j; }
		case "K15SL"   { $iks[15]= $j-1; printf "K15SL   at %d\n",$j; }
		case "K16SL"   { $iks[16]= $j-1; printf "K16SL   at %d\n",$j; }
		case "K17SL"   { $iks[17]= $j-1; printf "K17SL   at %d\n",$j; }
		case "K18SL"   { $iks[18]= $j-1; printf "K18SL   at %d\n",$j; }
		case "K19SL"   { $iks[19]= $j-1; printf "K19SL   at %d\n",$j; }
		case "K20SL"   { $iks[20]= $j-1; printf "K20SL   at %d\n",$j; }
	    }
	}
    }
    if( ($buf[0] ne '@') && ($buf[0] ne '$') && ($buf[0] ne '*') ){
	if( $fswitch == 0 ){ printf "No format line in input file, exiting \n"; exit 1;}
	$nelm1[$n] =$n; 
        $buf[$iname] =~ tr/\"//d;
        $buf[$iname] =~ tr/\_//d;
        $buf[$iname] =~ s/\.//g;
        $buf[$iname] =~ s/\$//g;
        $buf[$iname] =~ s/B1//g;
        $buf[$iname] =~ s/B2//g;
	$buf[$ikey]  =~ tr/\"//d;
	$type[$n]  =$buf[$ikey];
	switch ($type[$n]) {
	    case "DRIFT"     { $name1[$n]="$buf[$iname].$n"; $nameD[$ndrift]=$name1[$n]; 
			       $ldrift[$ndrift]=$buf[$il];
			       $ndrift=$ndrift+1; 
			     }
	    case "BEAMBEAM"  { $name1[$n]=$buf[$iname];        $nameI[$nip]= $name1[$n];  
			       $x1[$nip]     = $buf[$ix];    $px1[$nip]    = $buf[$ipx];
			       $betax1[$nip] = $buf[$ibetx]; $alfax1[$nip] = $buf[$ialfx];
			       $mux1[$nip]   = $buf[$imux];  $dx1[$nip]    = $buf[$idx];
			       $dpx1[$nip]   = $buf[$idpx];
			       $y1[$nip]     = $buf[$iy];       $py1[$nip] = $buf[$ipy];
			       $betay1[$nip] = $buf[$ibety]; $alfay1[$nip] = $buf[$ialfy];
			       $muy1[$nip]   = $buf[$imuy];     $dy1[$nip] = $buf[$idy];
			       $dpy1[$nip]   = $buf[$idpy]; $nip = $nip+1;
			     }
	    case "RFCAVITY"  { $name1[$n]="$buf[$iname].$n"; $nameR[$nrf] = $name1[$n]; 
			       $lrf[$nrf]  = $buf[$il];       $volt[$nrf] = $buf[$ivolt]; 
			       $freq[$nrf] = $buf[$ifreq];     $lag[$nrf] = $buf[$ilag]; 
                               if( $volt[$nrf] == 0 ){ $n=$n-1; }
			       else{ $nrf=$nrf+1; }
			     }
	    case "SOLENOID"  { $name1[$n]="$buf[$iname].$n"; $nameS[$nsol]=$name1[$n]; 
			       $ksi[$nsol]=$buf[$iksi]; $ks[$nsol]=$buf[$iksi]/$buf[$ilrad];
                               if( $ksi[$nsol] == 0 ){ $n=$n-1; }
			       else{ $nsol=$nsol+1; }
			     }
	    case "DIPEDGE"   { $name1[$n]="$buf[$iname].$n"; $nameDE[$ndpdg]=$name1[$n]; 
			       $h1[$ndpdg]=$buf[$ih1];           $e1[$ndpdg]=$buf[$ie1];
			       $fint[$ndpdg]=$buf[$ifint];     $hgap[$ndpdg]=$buf[$ihgap];
			       $ndpdg=$ndpdg+1; 
			     }
	    case "KICKER"    { $name1[$n]="$buf[$iname].$n"; $nameK[$nkick]=$name1[$n]; 
			       $hkick[$nkick]=$buf[$ihkick]; $vkick[$nkick]=$buf[$ivkick];
                               if( $hkick[$nkick]==0 && $vkick[$nkick]==0 ){ $n=$n-1; }
			       else{ $nkick=$nkick+1; }
			     }
	    case "HKICKER"   { $name1[$n]="$buf[$iname].$n"; $nameK[$nkick]=$name1[$n]; 
			       $hkick[$nkick]=$buf[$ihkick]; $vkick[$nkick]=$buf[$ivkick];
                               if( $hkick[$nkick]==0 && $vkick[$nkick]==0 ){ $n=$n-1; }
			       else{ $nkick=$nkick+1; }
			     }
	    case "VKICKER"   { $name1[$n]="$buf[$iname].$n"; $nameK[$nkick]=$name1[$n]; 
			       $hkick[$nkick]=$buf[$ihkick]; $vkick[$nkick]=$buf[$ivkick];
                               if( $hkick[$nkick]==0 && $vkick[$nkick]==0 ){ $n=$n-1; }
			       else{ $nkick=$nkick+1; }
       	                     }
	    case "MULTIPOLE" { $name1[$n]="$buf[$iname].$n"; $nameM[$nmult] = $name1[$n];
			       $lrad[$nmult] = $buf[$ilrad]; 
			       $xM[$nmult]   = $buf[$ix];    $yM[$nmult]    = $buf[$iy];
			       for($i=0;$i<21;$i++){ $knl[$nmult][$i]=$buf[$ikl[$i]]; }
			       $nn[$nmult]=20; $j=1;
			       do{ if($knl[$nmult][$nn[$nmult]]==0){
				   $nn[$nmult]=$nn[$nmult]-1;}else{$j=0;} 
				   if( $nn[$nmult]<0 ){$j=0;}
			       }while($j !=0 ); 
			       for($i=0;$i<21;$i++){ $ksl[$nmult][$i]=$buf[$iks[$i]]; }
			       $ns[$nmult]=20; $j=1;
			       do{ if($ksl[$nmult][$ns[$nmult]]==0){
				   $ns[$nmult]=$ns[$nmult]-1;}else{$j=0;}
				   if( $ns[$nmult]<0 ){$j=0;}
			       }while($j !=0 ); 
                               if( $nn[$nmult]<0 && $ns[$nmult]<0 ){ $n=$n-1; }
                               else{ $nmult=$nmult+1; }
			     }
	}
	$s1[$n]    =$buf[$is];
	$n=$n+1;
    }
}
$n1=$n;
printf "Number of Lines  read from lattice: %d \n", $n1;
printf "Number of IPs                     : %d \n", $nip;
printf "Number of DRIFTs                  : %d \n", $ndrift;
printf "Number of MULTs                   : %d \n", $nmult;
printf "Number of RFs                     : %d \n", $nrf;
printf "Number of SOLENOIDs               : %d \n", $nsol;
printf "Number of DIPEDGEs                : %d \n", $ndpdg;
printf "Number of correctors (KICKER)     : %d \n", $nkick;
close(fpr1);
#
#--- 
printf "Reading multipole errors file...\n";
$n=0;
while( <fpr2> ){
    @buf=split ;
    if( $buf[0] eq '*' ){
	printf "Input format:\n";
	$fswitch=1;
	$bufsize=scalar(@buf);
	for($j=1;$j<$bufsize;$j++){
	    switch ($buf[$j]) {
		case "NAME"    { $iname = $j-1; printf "NAME    at %d\n",$j; }
		case "K0L"     { $ikl[0] = $j-1; printf "K0L     at %d\n",$j; }
		case "K1L"     { $ikl[1] = $j-1; printf "K1L     at %d\n",$j; }
		case "K2L"     { $ikl[2] = $j-1; printf "K2L     at %d\n",$j; }
		case "K3L"     { $ikl[3] = $j-1; printf "K3L     at %d\n",$j; }
		case "K4L"     { $ikl[4] = $j-1; printf "K4L     at %d\n",$j; }
		case "K5L"     { $ikl[5] = $j-1; printf "K5L     at %d\n",$j; }
		case "K6L"     { $ikl[6] = $j-1; printf "K6L     at %d\n",$j; }
		case "K7L"     { $ikl[7] = $j-1; printf "K7L     at %d\n",$j; }
		case "K8L"     { $ikl[8] = $j-1; printf "K8L     at %d\n",$j; }
		case "K9L"     { $ikl[9] = $j-1; printf "K9L     at %d\n",$j; }
		case "K10L"    { $ikl[10]= $j-1; printf "K10L    at %d\n",$j; }
		case "K11L"    { $ikl[11]= $j-1; printf "K11L    at %d\n",$j; }
		case "K12L"    { $ikl[12]= $j-1; printf "K12L    at %d\n",$j; }
		case "K13L"    { $ikl[13]= $j-1; printf "K13L    at %d\n",$j; }
		case "K14L"    { $ikl[14]= $j-1; printf "K14L    at %d\n",$j; }
		case "K15L"    { $ikl[15]= $j-1; printf "K15L    at %d\n",$j; }
		case "K16L"    { $ikl[16]= $j-1; printf "K16L    at %d\n",$j; }
		case "K17L"    { $ikl[17]= $j-1; printf "K17L    at %d\n",$j; }
		case "K18L"    { $ikl[18]= $j-1; printf "K18L    at %d\n",$j; }
		case "K19L"    { $ikl[19]= $j-1; printf "K19L    at %d\n",$j; }
		case "K20L"    { $ikl[20]= $j-1; printf "K20L    at %d\n",$j; }
		case "K0SL"    { $iks[0] = $j-1; printf "K0SL    at %d\n",$j; }
		case "K1SL"    { $iks[1] = $j-1; printf "K1SL    at %d\n",$j; }
		case "K2SL"    { $iks[2] = $j-1; printf "K2SL    at %d\n",$j; }
		case "K3SL"    { $iks[3] = $j-1; printf "K3SL    at %d\n",$j; }
		case "K4SL"    { $iks[4] = $j-1; printf "K4SL    at %d\n",$j; }
		case "K5SL"    { $iks[5] = $j-1; printf "K5SL    at %d\n",$j; }
		case "K6SL"    { $iks[6] = $j-1; printf "K6SL    at %d\n",$j; }
		case "K7SL"    { $iks[7] = $j-1; printf "K7SL    at %d\n",$j; }
		case "K8SL"    { $iks[8] = $j-1; printf "K8SL    at %d\n",$j; }
		case "K9SL"    { $iks[9] = $j-1; printf "K9SL    at %d\n",$j; }
		case "K10SL"   { $iks[10]= $j-1; printf "K10SL   at %d\n",$j; }
		case "K11SL"   { $iks[11]= $j-1; printf "K11SL   at %d\n",$j; }
		case "K12SL"   { $iks[12]= $j-1; printf "K12SL   at %d\n",$j; }
		case "K13SL"   { $iks[13]= $j-1; printf "K13SL   at %d\n",$j; }
		case "K14SL"   { $iks[14]= $j-1; printf "K14SL   at %d\n",$j; }
		case "K15SL"   { $iks[15]= $j-1; printf "K15SL   at %d\n",$j; }
		case "K16SL"   { $iks[16]= $j-1; printf "K16SL   at %d\n",$j; }
		case "K17SL"   { $iks[17]= $j-1; printf "K17SL   at %d\n",$j; }
		case "K18SL"   { $iks[18]= $j-1; printf "K18SL   at %d\n",$j; }
		case "K19SL"   { $iks[19]= $j-1; printf "K19SL   at %d\n",$j; }
		case "K20SL"   { $iks[20]= $j-1; printf "K20SL   at %d\n",$j; }
	    }
	}
    }
    if( ($buf[0] ne '@') && ($buf[0] ne '$') && ($buf[0] ne '*') ){
	if( $fswitch == 0 ){ printf "No format line in input file, exiting \n"; exit 1;}
	$nmult1[$n] =$n; 
        $buf[$iname] =~ tr/\"//d;
        $buf[$iname] =~ tr/\_//d;
        $buf[$iname] =~ s/\.//g;
        $buf[$iname] =~ s/\$//g;
        $buf[$iname] =~ s/B1//g;
        $buf[$iname] =~ s/B2//g;
	$mult1[$n]=$buf[$iname];
#	$lrad[$nmult] = $buf[$ilrad]; 
	for($i=0;$i<21;$i++){ $knlm[$n][$i]=$buf[$ikl[$i]]; }
	  $nnm[$n]=20; $j=1;
	  do{ if($knlm[$n][$nnm[$n]]==0){$nnm[$n]=$nnm[$n]-1;}else{$j=0;} 
	      if( $nnm[$n]<0 ){$j=0;}
	  }while($j !=0 ); 
	for($i=0;$i<21;$i++){ $kslm[$n][$i]=$buf[$iks[$i]]; }
	  $nsm[$n]=20; $j=1;
	  do{ if($kslm[$n][$nsm[$n]]==0){$nsm[$n]=$nsm[$n]-1;}else{$j=0;}
	      if( $nsm[$n]<0 ){$j=0;}
	  }while($j !=0 ); 
         if( $nnm[$n]<0 && $nsm[$n]<0 ){ $n=$n-1; }
    $n=$n+1;
    }
}
$nm=$n;
printf "Number of multipole errors read: %d \n", $n;
close(fpr2);
#
#for($i=0;$i<$nm;$i++){
#    printf "%d %s\n",$i,$mult1[$i];
#    printf "KNL: ";
#    for($j=0;$j<21;$j++){printf "%f ",$knlm[$i][$j];}
#    printf "\n";
#    printf "KSL: ";
#    for($j=0;$j<21;$j++){printf "%f ",$kslm[$i][$j];}
#    printf "\n";
#}
#--- 
printf "Reading strong optics file...\n";
$n=0;
$fswitch=0;
while( <fpr3> ){
    @buf=split ;
    if( $buf[0] eq '*' ){
	printf "Input format:\n";
	$fswitch=1;
	$bufsize=scalar(@buf);
	for($j=1;$j<$bufsize;$j++){
	    switch ($buf[$j]) {
		case "NAME"    { $iname = $j-1; printf "NAME    at %d\n",$j; }
		case "KEYWORD" { $ikey  = $j-1; printf "KEYWORD at %d\n",$j; }
		case "S"       { $is    = $j-1; printf "S       at %d\n",$j; }
		case "X"       { $ix    = $j-1; printf "X       at %d\n",$j; }
		case "PX"      { $ipx   = $j-1; printf "PX      at %d\n",$j; }
		case "BETX"    { $ibetx = $j-1; printf "BETX    at %d\n",$j; }
		case "ALFX"    { $ialfx = $j-1; printf "ALFX    at %d\n",$j; }
		case "MUX"     { $imux  = $j-1; printf "MUX     at %d\n",$j; }
		case "DX"      { $idx   = $j-1; printf "DX      at %d\n",$j; }
		case "DPX"     { $idpx  = $j-1; printf "DPX     at %d\n",$j; }
		case "Y"       { $iy    = $j-1; printf "Y       at %d\n",$j; }
		case "PY"      { $ipy   = $j-1; printf "PY      at %d\n",$j; }
		case "BETY"    { $ibety = $j-1; printf "BETY    at %d\n",$j; }
		case "ALFY"    { $ialfy = $j-1; printf "ALFY    at %d\n",$j; }
		case "MUY"     { $imuy  = $j-1; printf "MUY     at %d\n",$j; }
		case "DY"      { $idy   = $j-1; printf "DY      at %d\n",$j; }
		case "DPY"     { $idpy  = $j-1; printf "DPY     at %d\n",$j; }
		case "L"       { $il    = $j-1; printf "L       at %d\n",$j; }
		case "LRAD"    { $ilrad = $j-1; printf "LRAD    at %d\n",$j; }
		case "VOLT"    { $ivolt = $j-1; printf "VOLT    at %d\n",$j; }
		case "LAG"     { $ilag  = $j-1; printf "LAG     at %d\n",$j; }
		case "FREQ"    { $ifreq = $j-1; printf "FREQ    at %d\n",$j; }
		case "K0L"     { $ikl[0] = $j-1; printf "K0L     at %d\n",$j; }
		case "K1L"     { $ikl[1] = $j-1; printf "K1L     at %d\n",$j; }
		case "K2L"     { $ikl[2] = $j-1; printf "K2L     at %d\n",$j; }
		case "K3L"     { $ikl[3] = $j-1; printf "K3L     at %d\n",$j; }
		case "K4L"     { $ikl[4] = $j-1; printf "K4L     at %d\n",$j; }
		case "K5L"     { $ikl[5] = $j-1; printf "K5L     at %d\n",$j; }
		case "K6L"     { $ikl[6] = $j-1; printf "K6L     at %d\n",$j; }
		case "K7L"     { $ikl[7] = $j-1; printf "K7L     at %d\n",$j; }
		case "K8L"     { $ikl[8] = $j-1; printf "K8L     at %d\n",$j; }
		case "K9L"     { $ikl[9] = $j-1; printf "K9L     at %d\n",$j; }
		case "K10L"    { $ikl[10]= $j-1; printf "K10L    at %d\n",$j; }
		case "K11L"    { $ikl[11]= $j-1; printf "K11L    at %d\n",$j; }
		case "K12L"    { $ikl[12]= $j-1; printf "K12L    at %d\n",$j; }
		case "K13L"    { $ikl[13]= $j-1; printf "K13L    at %d\n",$j; }
		case "K14L"    { $ikl[14]= $j-1; printf "K14L    at %d\n",$j; }
		case "K15L"    { $ikl[15]= $j-1; printf "K15L    at %d\n",$j; }
		case "K16L"    { $ikl[16]= $j-1; printf "K16L    at %d\n",$j; }
		case "K17L"    { $ikl[17]= $j-1; printf "K17L    at %d\n",$j; }
		case "K18L"    { $ikl[18]= $j-1; printf "K18L    at %d\n",$j; }
		case "K19L"    { $ikl[19]= $j-1; printf "K19L    at %d\n",$j; }
		case "K20L"    { $ikl[20]= $j-1; printf "K20L    at %d\n",$j; }
		case "K0SL"    { $iks[0] = $j-1; printf "K0SL    at %d\n",$j; }
		case "K1SL"    { $iks[1] = $j-1; printf "K1SL    at %d\n",$j; }
		case "K2SL"    { $iks[2] = $j-1; printf "K2SL    at %d\n",$j; }
		case "K3SL"    { $iks[3] = $j-1; printf "K3SL    at %d\n",$j; }
		case "K4SL"    { $iks[4] = $j-1; printf "K4SL    at %d\n",$j; }
		case "K5SL"    { $iks[5] = $j-1; printf "K5SL    at %d\n",$j; }
		case "K6SL"    { $iks[6] = $j-1; printf "K6SL    at %d\n",$j; }
		case "K7SL"    { $iks[7] = $j-1; printf "K7SL    at %d\n",$j; }
		case "K8SL"    { $iks[8] = $j-1; printf "K8SL    at %d\n",$j; }
		case "K9SL"    { $iks[9] = $j-1; printf "K9SL    at %d\n",$j; }
		case "K10SL"   { $iks[10]= $j-1; printf "K10SL   at %d\n",$j; }
		case "K11SL"   { $iks[11]= $j-1; printf "K11SL   at %d\n",$j; }
		case "K12SL"   { $iks[12]= $j-1; printf "K12SL   at %d\n",$j; }
		case "K13SL"   { $iks[13]= $j-1; printf "K13SL   at %d\n",$j; }
		case "K14SL"   { $iks[14]= $j-1; printf "K14SL   at %d\n",$j; }
		case "K15SL"   { $iks[15]= $j-1; printf "K15SL   at %d\n",$j; }
		case "K16SL"   { $iks[16]= $j-1; printf "K16SL   at %d\n",$j; }
		case "K17SL"   { $iks[17]= $j-1; printf "K17SL   at %d\n",$j; }
		case "K18SL"   { $iks[18]= $j-1; printf "K18SL   at %d\n",$j; }
		case "K19SL"   { $iks[19]= $j-1; printf "K19SL   at %d\n",$j; }
		case "K20SL"   { $iks[20]= $j-1; printf "K20SL   at %d\n",$j; }
	    }
	}
    }
    if( ($buf[0] ne '@') && ($buf[0] ne '$') && ($buf[0] ne '*') ){
	if( $fswitch == 0 ){ printf "No format line in input file, exiting \n"; exit 1;}
	$nelm2[$n] =$n; 
        $buf[$iname] =~ tr/\"//d;
        $buf[$iname] =~ tr/\_//d;
        $buf[$iname] =~ s/\.//g;
        $buf[$iname] =~ s/\$//g;
        $buf[$iname] =~ s/B1//g;
        $buf[$iname] =~ s/B2//g;
	$buf[$ikey]  =~ tr/\"//d;
	$type[$n]  =$buf[$ikey];
	switch ($type[$n]) {
	    case "BEAMBEAM"  { $name2[$n]=$buf[$iname];
			       $x2[$n]     = $buf[$ix];    $px2[$n]    = $buf[$ipx];
			       $betax2[$n] = $buf[$ibetx]; $alfax2[$n] = $buf[$ialfx];
			       $mux2[$n]   = $buf[$imux];  $dx2[$n]    = $buf[$idx];
			       $dpx2[$n]   = $buf[$idpx];
			       $y2[$n]     = $buf[$iy];    $py2[$n]    = $buf[$ipy];
			       $betay2[$n] = $buf[$ibety]; $alfay2[$n] = $buf[$ialfy];
			       $muy2[$n]   = $buf[$imuy];     $dy2[$n] = $buf[$idy];
			       $dpy2[$n]   = $buf[$idpy];
			     }
	}
	$n=$n+1;
    }
}
$n2=$n;
printf "Number of IPs read from strong optics file: %d \n", $n2;
close(fpr3);
for($j=0;$j<$n2;$j++){ printf "%s \n",$name2[$j]; }
#
#
if($nip != $n2 ) { printf "Numbers of IPs do not match. Exiting.\n"; exit(0);}
#
#--- Creating and writing output ---------------------------------
#
$kp=1.0E+2;
$km=1.0E-2;
printf fpw "File: lhc\n";
printf fpw "\n___________________________Working_Parameters___________________________\n";
printf fpw "Transport:  LIN\n";
printf fpw "Levels:  (X,Y)\n";
printf fpw "Gamma_weak:  %f                 (type)=%s \n",$gamma,$particle;
printf fpw "Emitt_str  (cm*rad):  (x)=%lG      (y)=%lG  \n",$emitx*$kp,$emity*$kp;
printf fpw "Emitt_weak (cm*rad):  (x)=%lG      (y)=%lG  \n",$emitx*$kp,$emity*$kp;
printf fpw "Sigma_str  (cm, ..):  (z)=%f        (dE/E)=%lG\n", $blength*$kp, $sige;
printf fpw "Sigma_weak (cm, ..):  (z)=%f        (dE/E)=%lG\n", $blength*$kp, $sige;
printf fpw "Aperture  (sigm):     (x)=20.           (y)=20.           (z)=20.\n";
printf fpw "Seed:            (comb_1)=89787    (comb_2)=325      (comb_3)=493\n";
printf fpw "Boundary:          (part)=100      (step)=1024\n";

#--- creating the structure list ------------------------
printf fpw "\nStructure:\n";
for($i=0;$i<$n1;$i++){
    printf fpw "%s ", $name1[$i];
    if( $name1[$i] =~ $mainIP){ printf fpw "Watch_point "; }
    if( ($i>0) && ($i % 4 == 0) ){ printf fpw "\n"; }
}
printf fpw "\nEnd_structure\n";
#
printf fpw "\n# Watch_point: WATCH\n";

#--- IPs ----------------------------------------
printf "Name_I Name_S Xsep Ysep Px Py\n";
for($i=0;$i<$nip;$i++){
    if($nameI[$i] =~ $mainIP){ 
	printf fpw "\n# ".$nameI[$i].": IP_BASE\n";
	printf fpw "Xi_bs:    (p)=%lG\n",-$Np; 
    }else{ 
	printf fpw "\n# ".$nameI[$i].": IP\n";
    }
    printf fpw "Latt_str: L_".$nameI[$i]."_s\n";
    printf fpw "Current: 1\n";
    printf fpw "Norm: Off\n";
    if($nameI[$i] eq $mainIP || $nameI[$i] eq $mainIP2 ){
	printf fpw "Lumi: On\n";
	printf fpw "Slices:  12\n";
    }else{
	printf fpw "Lumi: Off\n";
	printf fpw "Slices:  1\n";
    }
    for($j=0;$j<$n2;$j++){ if($name2[$j] eq $nameI[$i]){ $jj=$j; } }
#--- separation and angle --------
    $x[$i]=($x1[$i]-$x2[$jj])*$kp;
    $y[$i]=($y1[$i]-$y2[$jj])*$kp;
    $px[$i]=($px1[$i]-$px2[$jj]);
    $py[$i]=($py1[$i]-$py2[$jj]);
    if($nameI[$i] eq $mainIP || $nameI[$i] eq $mainIP2 ){
	printf fpw "Shift  (cm):         (x)=%G (y)=%G\n",$x[$i],$y[$i];
    }else{
	printf fpw "Shift  (cm):         (x)=%G (y)=%G (comp)=1\n",$x[$i],$y[$i];
    };
    printf fpw "Angle (rad):         (x)=%G (y)=%G\n",$px[$i],$py[$i];
#--- strong lattice parameters ------
    printf fpw "\n# L_".$nameI[$i]."_s: LATT\n";
    printf fpw "Beta   (cm):         (x)=%f (y)=%f\n",$betax2[$jj]*$kp,$betay2[$jj]*$kp;
    printf fpw "Alpha:               (x)=%f (y)=%f\n",$alfax2[$jj],$alfay2[$jj];
    printf fpw "Disp   (cm):         (x)=%f (y)=%f\n",$dx2[$jj]*$kp,$dy2[$jj]*$kp;
    printf fpw "Disp_drv:            (x)=%f (y)=%f\n",$dpx2[$jj],$dpy2[$jj];

    printf "%s %s %G %G %G %G\n", $nameI[$i], $name2[$jj],$x[$i],$y[$i],$px[$i],$py[$i];
}
#--- DRIFTs -------------------------------------
$circumf=0;
for($i=0;$i<$ndrift;$i++){
    printf fpw "\n# ".$nameD[$i].": DRIFT\n";
    printf fpw "Drift (cm): %lG (dz)=1\n", $kp*$ldrift[$i];
    $circumf=$circumf+$ldrift[$i];
}
printf "Machine length = %12.10f m \n",$circumf;
#
#--- RF cavities ---------------------------------
for($i=0;$i<$nrf;$i++){
    printf fpw "\n# ".$nameR[$i].": EXT_RFCV\n";
    printf fpw "Value_1: %G \n", $volt[$i]/$energy/1.0E3*$beta;
    printf fpw "Value_2: %G \n", 2.99792458E10*$beta/$freq[$i]/1.0E6;
}
#
#--- SOLENOIDs  ---------------------------------
for($i=0;$i<$nsol;$i++){
    printf fpw "\n# ".$nameS[$i].": EXT_SOLE\n";
    printf fpw "Value_1: %G \n", $ks[$i]*$km;
    printf fpw "Value_2: %G \n", $ksi[$i];
}
#
#--- DIPEDGEs  ---------------------------------
for($i=0;$i<$ndpdg;$i++){
    printf fpw "\n# ".$nameDE[$i].": EXT_DPDG\n";
    printf fpw "Value_1: %G \n", $h1[$i]*$km;
    printf fpw "Value_2: %G \n", $e1[$i];
    printf fpw "Value_3: %G \n", $fint[$i];
    printf fpw "Value_4: %G \n", $hgap[$i]*$kp;
}
#
#--- KICKERs  ---------------------------------
for($i=0;$i<$nkick;$i++){
    printf fpw "\n# ".$nameK[$i].": EXT_CORR\n";
    printf fpw "Value_1: %G \n", $hkick[$i];
    printf fpw "Value_2: %G \n", $vkick[$i];
}
#
#--- MULTs ---------------------------------------
for($i=0;$i<$nmult;$i++){
    $err=-1;
    $sname=substr($nameM[$i],0,index($nameM[$i],'.'));
    for($j=0;$j<$nm;$j++){
	if( $mult1[$j] eq $sname ){
	    printf "attaching multipole errors %s to %s\n",$mult1[$j],$nameM[$i];
	    $err=$j;
	}
    }
    printf fpw "\n# ".$nameM[$i].": MULT\n";
    $zero=0;
    printf fpw "Shift: (x)=%G (y)=%G\n", $xM[$i]*$kp*$zero, $yM[$i]*$kp*$zero;
    printf fpw "Length: %G \n", $lrad[$i]*$kp;
    if( $nn[$i] >= 0 || $nnm[$err] >=0 ){printf fpw "KNL: ";}
    for($j=0;$j<max($nn[$i]+1,$nnm[$err]+1);$j++){
	printf fpw "%G \n",$knl[$i][$j]*$km**$j+$knlm[$err][$j]*$km**$j;
    }
    if( $ns[$i] >= 0 || $nsm[$err] >=0 ){printf fpw "KSL: ";}
    for($j=0;$j<max($ns[$i]+1,$nsm[$err]+1);$j++){
	printf fpw "%G \n",$ksl[$i][$j]*$km**$j+$kslm[$err][$j]*$km**$j;
    }
#    if( ($ns[$i] < 0) && ($nn[$i] < 0) ){printf fpw "KNL: 0\n";}
}
#
printf fpw "\n_______________________End_Working_Parameters___________________________\n";
#
close(fpw);
