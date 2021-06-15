#!/bin/bash

if [ $# -ne 1 ]; then
echo Usage:  normalization_tipe  
exit 1 
fi

normalization_tipe=$1
D_maxstrenght=1
n_samples=50
delta_SEED=1000

#declare -a NArray=(500 400 300) 
#declare -a AlphaArray=(0.65 0.64 0.63 0.62 0.61 0.6 0.59 0.58 0.57 0.56 0.55 0.54 0.53 0.52 0.51 0.5 0.45 0.4 0.35 0.3) 
declare -a NArray=(20) 
declare -a AlphaArray=(0.1)

#rm scr_instructions.txt
FNAME=scr_instructions$(date '+%d:%m:%Y:%H:%M:%S').txt
LOGNAME=scr_instructions$(date '+%d:%m:%Y:%H:%M:%S').log

for N in ${NArray[@]}; do
strenghtN=$(gawk -v N=$N 'BEGIN{print 0.01*N}')
    for a in ${AlphaArray[@]}; do
    
echo ./scr_hu_N${N}a${a}S${n_samples}.exe $N $a $strenghtN $D_maxstrenght $normalization_tipe $n_samples $delta_SEED>> ./${FNAME}
	cp hebbian_unlearning_v2.4.exe scr_hu_N${N}a${a}S${n_samples}.exe
	delta_SEED=delta_SEED+500
	done
done
	
cat ./${FNAME} | xargs -n 10 -P 5 time >> ./${LOGNAME} 2>&1 &
#

