#!/bin/bash

prefix=$1
ni=$2
nf=$3

path_frames="frames"
path_mdps="../mdp"
path_gmx="/usr/local/bin/"
outout_eigenfreqs_filename="$prefix""_eigenfreq.dat"
temperature=309


rm "$outout_eigenfreqs_filename" 2>/dev/null

let total_frames=$nf-$ni
let total_frames=$total_frames+1

total_atoms=$(head -n 2 "$path_frames""/""$prefix"".0.gro" | tail -n 1 | awk '{print $1}')
let total_modes=$total_atoms*3

n=$ni

while [ $n -le $nf ]; do
	if [ ! -e "$prefix""_min_nma"."$n"".gro" ]; then
	
		echo "System" | "$path_gmx""/"gmx_d editconf -f "$path_frames""/""$prefix"".""$n"".gro" -o "$prefix""_box.""$n"".gro" -c -d 2.0 -princ
		
		"$path_gmx""/"gmx_d grompp -f "$path_mdps"/"min_nma"".mdp" -c "$prefix""_box.""$n"".gro" -p "$prefix"".top" -o "$prefix""_min_nma"."$n"".tpr"
		"$path_gmx""/"gmx_d mdrun -s "$prefix""_min_nma"."$n"".tpr" -v -deffnm "$prefix""_min_nma"."$n" -nt 1
		
		# Sometimes the mininimization stops with segmentation fault.
		# When this happens, the .gro file is not generated. Thus, generate it using the trr file
		if [ ! -e "$prefix""_min_nma"."$n"".gro" ]; then
			echo 0 0 | gmx rms -f "$prefix""_min_nma"."$n"".trr" -s "$prefix""_min_nma"."$n"".tpr" -o "temporary_rmsd_""$prefix""_min_nma"."$n"".xvg"
			last_step=$(tail -n 1 "temporary_rmsd_""$prefix""_min_nma"."$n"".xvg" | awk '{printf "%d\n",$1}')
			echo 0 | gmx trjconv -f "$prefix""_min_nma"."$n"".trr" -s "$prefix""_min_nma"."$n"".tpr" -o "$prefix""_min_nma"."$n"".gro" -b "$last_step" -e "$last_step"
		fi
	fi
	
	
	
	
	
	if [ ! -e "$prefix""_eigenfreq.""$n"".xvg" ]; then
		"$path_gmx""/"gmx_d grompp -f "$path_mdps"/"nma.mdp" -c "$prefix""_min_nma"."$n"".gro" -t "$prefix""_min_nma"."$n"".trr" -p "$prefix"".top" -o "$prefix""_nma.""$n"".tpr"
		"$path_gmx""/"gmx_d mdrun -s "$prefix""_nma.""$n"".tpr" -v -deffnm "$prefix""_nma.""$n" -nt 1

		"$path_gmx""/"gmx_d nmeig -f "$prefix""_nma.""$n"".mtx" -s "$prefix""_nma.""$n"".tpr" -of "$prefix""_eigenfreq.""$n"".xvg" -ol "$prefix""_eigenval.""$n"".xvg" -v "$prefix""_eigenvec.""$n"".trr" -T "$temperature" -first 1 -last $total_modes
	
		maxforce_min=$(grep "Maximum force" "$prefix""_min_nma.""$n"".log" | awk '{print $4}')
		maxforce_nma=$(grep "Maximum force" "$prefix""_nma.""$n"".log" | awk '{print $3}')
		echo -e "$n""\t""$maxforce_min""\t""$maxforce_nma" >> "$prefix""_min_nma_forces.log"
		
		rm mdout.mdp \#*
		rm "$prefix""_box.""$n"".gro" "$prefix""_nma.""$n"* "$prefix""_eigenvec.""$n"".trr" "$prefix""_eigenval.""$n"".xvg"
#		rm "$prefix""_min_nma.""$n"*
	
		for i in $(grep -v "#" "$prefix""_eigenfreq.""$n"".xvg" | grep -v "@" | awk '{print $2}') ; do
			echo -n -e "$i""\t" >> "$outout_eigenfreqs_filename"
		done
		echo "" >> "$outout_eigenfreqs_filename"
#		rm "$prefix""_eigenfreq.""$n"".xvg"
	fi
	
	
	let n=$n+1
done

#./entropy_vibrational "$prefix" "$total_atoms" "$total_frames"

