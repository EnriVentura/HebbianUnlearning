#!/bin/bash

if [ $# -ne 3 ]; then
echo Usage: N NORM_TYPE NSAMPLES
exit 1 
fi
declare -a AlphaArray=(1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1) 
declare -a strenghtNArray=(0.01)
declare -a D_maxstrenghtArray=(5)

N=$1
NORM_TYPE=$2
NSAMPLES=$3

FNAME=scr_instructions$(date '+%d:%m:%Y:%H:%M:%S').txt
LOGNAME=scr_instructions$(date '+%d:%m:%Y:%H:%M:%S').log
 
for D_maxstrenght in ${D_maxstrenghtArray[@]}; do
for strenghtN in ${strenghtNArray[@]}; do
for a in ${AlphaArray[@]}; do
    P=$(awk -v alpha="$a" -v N="$N" 'BEGIN {print int(alpha*N)}')
echo ./scr_unlearning${NORM_TYPE}_N${N}P${P}a${a}strenghtN${strenghtN}D_maxstrenght${D_maxstrenght}NSAMPLES${NSAMPLES}.exe $N $a $strenghtN ${D_maxstrenght} ${NORM_TYPE} ${NSAMPLES} >> ./${FNAME}
cp hebbian_unlearning_v2.exe scr_unlearning${NORM_TYPE}_N${N}P${P}a${a}strenghtN${strenghtN}D_maxstrenght${D_maxstrenght}NSAMPLES${NSAMPLES}.exe
done
done
done
	
cat ./${FNAME} | xargs -n 7 -P 5 time >> ./${LOGNAME} 2>&1 &
#

