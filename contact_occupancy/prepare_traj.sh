#!/bin/bash

destination="dcsign_3ca"
destination=$1
path="/media/leandro/lob_hd1/_OK_/ctl/calcium/data/""$destination""/analysis"
here=$(pwd)

grp_protein="Protein"
grp_ligand="CA"

mkdir "$destination"
cd "$destination"

#echo -e '"Protein" | "CA"'"\n"'q'"\n" | make_ndx -f "$path"/md1_non-water.tpr -o index.ndx
#echo -e '"Protein" | "NA"'"\n"'q'"\n" | make_ndx -f "$path"/md1_non-water.tpr -o index.ndx

mkdir "frames"

n=1
out=0
while [ $n -le 3 ]; do
	
	mkdir frames/"$n"
	
	echo "$grp_protein" | gmx trjconv -f "$path"/md"$n"_non-water.xtc -s "$path"/md"$n"_non-water.tpr -o prot..pdb -sep -b 100000 -e 200000 -dt 100 #-n index.ndx
	echo "$grp_ligand" | gmx trjconv -f "$path"/md"$n"_non-water.xtc -s "$path"/md"$n"_non-water.tpr -o lig..pdb -sep -b 100000 -e 200000 -dt 100 #-n index.ndx
	
	
	total_frames=$(ls -lh lig.*.pdb | wc -l | awk '{print $1}')
	
	frame=0
	out=0	#####
	while [ $frame -lt $total_frames ]; do
		sed 's/HETATM/ATOM/g' prot."$frame".pdb | grep "ATOM" | cut -c 23-26,30-55 > frames/"$n"/frame."$out".pdb
		sed 's/HETATM/ATOM/g' lig."$frame".pdb | grep "ATOM" | cut -c 23-26,30-55 >> frames/"$n"/frame."$out".pdb
		
		let out=$out+1
		let frame=$frame+1
	done
	
	if [ $n -eq 1 ]; then	# If this is the first repetition, get some information from the first pdb file
		atm_protein=$(sed 's/HETATM/ATOM/g' prot.0.pdb | grep "ATOM" | wc -l | awk '{print $1}')
		atm_ligand=$(sed 's/HETATM/ATOM/g' lig.0.pdb | grep "ATOM" | wc -l | awk '{print $1}')
		res_protein=$(sed 's/HETATM/ATOM/g' prot.0.pdb | grep "ATOM" | awk '{print $3}' | grep ^"N"$ | wc -l | awk '{print $1}')
	fi
	
	rm prot.*.pdb lig.*.pdb
	
	let n=$n+1
done

let total_frames=$out-1


echo -e "$total_frames\n$atm_protein\n$atm_ligand\n$res_protein" > "parameters.dat"

cd "$here"
