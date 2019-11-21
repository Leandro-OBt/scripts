#!/bin/bash

rm temporary_* log_* 2>/dev/null
rm -r frames_min_nma/
mkdir frames_min_nma eigenfreqs 2>/dev/null
mv *eigenfreq*xvg eigenfreqs/

for grp in complex receptor ligand ; do
	input_trr=''
	input_t=''
	n=0
	while [ $n -le 100 ]; do
		echo -en "$grp"."$n""\t"
		echo 'System System' | gmx_d rms -f "$grp"_min_nma.$n.trr -s "$grp"_min_nma.$n.tpr -o temporary.$n.xvg >/dev/null 2>/dev/null
		last_step=$(cat temporary.$n.xvg | awk '{printf "%d\n", $1}' | sort -n | tail -n 1)
		rm temporary.$n.xvg
		echo -en "$last_step""\t"

		echo 'System' | gmx_d trjconv -f "$grp"_min_nma.$n.trr -s "$grp"_min_nma.$n.tpr -o frames_min_nma/"$grp"_min_nma.$n.trr -b $last_step -e $last_step >/dev/null 2>/dev/null
		echo 'ok'

		input_trr="$input_trr"' '"$grp"'_min_nma.'"$n"'.trr'
		input_t="$input_t""$n""\n"

		let n=$n+1
	done
	echo -en 'joining'"\t"
	cd frames_min_nma/
	echo -e $input_t | gmx_d trjcat -f $input_trr -o ../"$grp"_min_nma.trr -settime #>/dev/null 2>/dev/null
	cd ../
	echo 'ok'
done
#rm -r frames_min_nma/
