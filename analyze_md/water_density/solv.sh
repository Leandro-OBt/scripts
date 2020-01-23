#!/bin/bash

grp="OW"
input_filename_prefix="frame_""$grp"
output_filename_prefix="water-density_K"

flag_frames_pdb=1
flag_franes_dat=1
flag_dg=1

if [ $flag_frames_pdb -eq 1 ]; then
	rep=1
	while [ $rep -le 3 ]; do
		echo "$grp" | gmx trjconv -f ../md"$rep"_fit.xtc -s ../../md"$rep".tpr -n index.ndx -o temporary_frame_"$grp"."$rep"..pdb -sep -b 100000 -e 200000 -dt 100
		let rep=$rep+1
	done
fi





if [ $flag_franes_dat -eq 1 ]; then
	frame=1
	for filename in temporary_frame_"$grp"* ; do
		grep "ATOM " "$filename" | awk '{print $6,$7,$8}' > "$input_filename_prefix""_""$frame"".dat"
		
		let frame=$frame+1
			
	done
fi




if [ $flag_dg -eq 1 ]; then
	total_atm=$(wc -l "$input_filename_prefix""_1.dat" | awk '{print $1}')
	total_frame=$(find . -name "$input_filename_prefix""_*.dat" | wc -l | awk '{print $1}')
	/home/leandro/repos/scripts/analyze_md/water_density/solv "$input_filename_prefix" $total_atm $total_frame "$output_filename_prefix"
fi





#for stdev in $(awk '{print $2}' "$output_filename_prefix"".pdb" | sort -n | uniq); do
#	grep " ""$stdev"" " "$output_filename_prefix"".pdb" > "$output_filename_prefix""_""$stdev""stdev.pdb"
#done
