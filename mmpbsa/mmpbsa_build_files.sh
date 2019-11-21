#!/bin/bash

# -------------------
#
#     This script will create the files needed to run mmpbsa.sh
#     The main activity is to make a single trajectory which contains only the relevant frames for the MM/PBSA calculation from a folder with many complete trajectories
#     In order to do this:
#          - Leave only the trajectories that will be used in a given directory
#          - Run this script
#     If time_e and time_b must be different between trajectories, then an input list approach must be used
#     The key parameter is number_of_frames. This is the number of frames you want to use in the MM/PBSA calculation, i.e. the number of frames that the final trajectory will have.
#
#
# Leandro Olivera Bortot
# 05 Oct 2016
#
# -------------------


input_trj_folder=$1		# ../
input_tpr=$2			# ../md_non-water.1.tpr
#time_dt_mmpbsa=$3
number_of_frames=$3		# 100
time_b=$4				# 10000
time_e=$5				# 20000
receptor=$6			
ligand=$7				
index=$8				# if this file doesnt exist, it will be created

trj_name="md_non-water"	# Used to find the files in the directory.
output_name="mmpbsa"	# Used to name the output files

path_gmx='/home/'"$USER"'/programs/gmx-5.1.2/bin/'
path_gmx='/usr/local/bin'

let time_dt_trj=$time_e-$time_b
string_settime=""
string_trjname=""

n=1
time=0
time_total=0

for trj in $(find "$input_trj_folder"/ -maxdepth 1 -name "$trj_name"".*.xtc" | sed 's/\//\ /g' | awk '{print $NF}' ) ; do



	echo 0 | "$path_gmx"/gmx trjconv -f "$input_trj_folder"/"$trj" -s "$input_tpr" -o temporary_trj."$n".xtc -b $time_b -e $time_e >/dev/null 2>/dev/null
	
	string_settime="$string_settime""$time"'\n'
	string_trjname="$string_trjname""temporary_trj.""$n"".xtc "
	
	let time=$time+$time_dt_trj
	let time_total=$time_total+$time_dt_trj
	
	let n=$n+1

done

echo -e "$string_settime" | "$path_gmx"/gmx trjcat -f $string_trjname -o temporary_trj.xtc -settime >/dev/null 2>/dev/null
rm temporary_trj.*.xtc

let time_dt_mmpbsa=$time_total/$number_of_frames
echo 'time total = '$time_total
echo 'time dt    = '$time_dt_mmpbsa
if [ ! -e "$index" ]; then
	echo -e '"'"$receptor"'" | "'"$ligand"'"'"\n"'q'"\n" | "$path_gmx"/gmx make_ndx -f "$input_tpr" -o "$index" >/dev/null 2>/dev/null
fi
complex="$receptor"'_'"$ligand"
echo "$complex" | "$path_gmx"/gmx trjconv -f temporary_trj.xtc -s "$input_tpr" -n "$index" -dt "$time_dt_mmpbsa" -o "$output_name".pdb >/dev/null 2>/dev/null
#echo "$complex" | "$path_gmx"/gmx convert-tpr -s "$input_tpr" -n "$index"  -o "$output_name".tpr >/dev/null 2>/dev/null
#echo -e 'q'"\n" | "$path_gmx"/gmx make_ndx -f "$output_name".pdb -o "$output_name".ndx >/dev/null 2>/dev/null

#echo 'System' | "$path_gmx"/gmx trjconv -f "$output_name".xtc -s "$output_name".tpr -n "$output_name".ndx -o "$output_name".pdb #>/dev/null 2>/dev/null

rm temporary_trj.xtc




#for dir in dcsign_gl*/; do cd $dir ; mkdir mmpbsa ; cd eq_smd ; echo -e '"Protein-H" | "Other"'"\n"'q'"\n" | ~/programs/gmx-5.1.2/bin/gmx make_ndx -f ../md.1.tpr -o temporary.ndx ; rm mmpbsa.pdb ; n=1 ; for i in md-frame*gro ; do echo 'Protein-H_Other' | ~/programs/gmx-5.1.2/bin/gmx trjconv -s $i -f $i -n temporary.ndx -o temporary.pdb ; echo 'MODEL '"$n" >> mmpbsa.pdb ; grep 'ATOM' temporary.pdb >> mmpbsa.pdb ; echo 'ENDMDL' >> mmpbsa.pdb ; rm temporary.pdb ; let n=$n+1 ; done ; echo 'Protein-H_Other' | ~/programs/gmx-5.1.2/bin/gmx convert-tpr -s ../md.1.tpr -n temporary.ndx -o ../mmpbsa/mmpbsa.tpr ; mv mmpbsa.pdb ../mmpbsa/ ; rm temporary.ndx ; cd ../mmpbsa/ ; echo 'System' | ~/programs/gmx-5.1.2/bin/gmx trjconv -s mmpbsa.tpr -f mmpbsa.pdb -o mmpbsa.xtc ; echo -e '"Protein-H" | "Other"'"\n"'q'"\n" | ~/programs/gmx-5.1.2/bin/gmx make_ndx -f mmpbsa.tpr -o mmpbsa.ndx ; cd ../ ; cd ../ ; done
