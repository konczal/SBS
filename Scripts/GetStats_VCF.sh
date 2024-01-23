VCF=$1
ID=$2

MIS=`grep -v "^#" $VCF | grep "missense_variant" | wc -l`
SYN=`grep -v "^#" $VCF | grep -v "missense_variant" | grep "synonymous_" | wc -l`
STOP_G=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" | grep "stop_gain"| wc -l`
START_L=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" |  grep -v "stop_gain" | grep "start_lost" | wc -l`
STOP_L=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" |  grep -v "stop_gain" | grep -v "start_lost" | grep "stop_lost" | wc -l`
START_R=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" |  grep -v "stop_gain" | grep -v "start_lost" | grep -v "stop_lost" | grep "start_retained" | wc -l`
STOP_R=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" |  grep -v "stop_gain" | grep -v "start_lost" | grep -v "stop_lost" | grep -v "start_retained" | grep "stop_retained" | wc -l`
OTHER=`grep -v "^#" $VCF | grep -v "missense_variant" | grep -v "synonymous_" |  grep -v "stop_gain" | grep -v "start_lost" | grep -v "stop_lost" | grep -v "start_retained" | grep -v "stop_retained" | wc -l`


echo $ID "\tmisense_variants\t" $MIS
echo $ID "\tmynonymous_variants\t" $SYN
echo $ID "\tstop_gained_variants\t" $STOP_G
echo $ID "\tstart_lost_variants\t" $START_L
echo $ID "\tstop_lost_variants\t" $STOP_L
echo $ID "\tstart_retained_variants\t" $START_R
echo $ID "\tstop_retained_variants\t" $STOP_R
echo $ID "\tother_variants\t" $OTHER
