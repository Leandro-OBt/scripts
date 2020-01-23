#!/bin/bash

glycosylation=$1

time_b=25000
time_e=50000
time_dt=10

if [ $glycosylation -eq 5 ]; then
	index_res='/media/'"$USER"'/lob_hd/_OK_/ctl/glycosylation/fit/dcsign/index_res_glcnac2man5.ndx'
else
	index_res='/media/'"$USER"'/lob_hd/_OK_/ctl/glycosylation/fit/build_9/index_res_glcnac2man9.ndx'
fi
dcsign_renum='/media/'"$USER"'/lob_hd/_OK_/ctl/glycosylation/fit/dcsign/dcsign_renumbered.pdb'
output_directory='hydrogen_bond_profile'


find unbind/ -name "*.xtc" | awk -F '.' '{print $2}' > temporary_unbind

string_trj=''

mkdir "$output_directory" 2>/dev/null
rm "$output_directory"/temporary* 2>/dev/null

for rep in $(echo '1 2 3 4 5'); do
	
	a=$(grep $rep temporary_unbind)
	if [ -z $a ]; then
		string_trj="$string_trj"' '"$rep"
	fi
	
done
rm temporary_unbind


if [ ! -z "$string_trj" ]; then



total_grp=$(grep '\[ ' "$index_res" | wc -l | awk '{print $1-1}')

for rep in $(echo "$string_trj") ; do 
	
	echo 'Protein Other' | gmx hbond -f md_non-water."$rep".xtc -s md_non-water."$rep".tpr -num "$output_directory"/hb."$rep".xvg -xvg none -b $time_b -e $time_e -dt $time_dt >/dev/null 2>/dev/null
	
	points_total=$(wc -l "$output_directory"/hb."$rep".xvg | awk '{print $1}')
	points_zero=$(grep '   0   ' "$output_directory"/hb."$rep".xvg | wc -l | awk '{print $1}')
	
	hb_num=$(gmx analyze -f "$output_directory"/hb."$rep".xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\n", $2}')
	hb_occ=$(echo $points_total $points_zero | awk '{printf "%.4f\n", 1.000000-($2/$1)}')
	
	echo "$hb_num" >> "$output_directory"/temporary_hb_num.xvg
	echo "$hb_occ" >> "$output_directory"/temporary_hb_occ.xvg
	
	grp=1
	while [ $grp -le $total_grp ]; do
		
		if [ ! -e "$output_directory"/hb_r"$grp".xvg ]; then
			echo '0 '"$grp" | gmx hbond -f md_non-water."$rep".xtc -s md_non-water."$rep".tpr -num "$output_directory"/hb_r"$grp"."$rep".xvg -xvg none -b $time_b -e $time_e -dt $time_dt -n "$index_res" >/dev/null 2>/dev/null
		fi
		
		points_total=$(wc -l "$output_directory"/hb_r"$grp"."$rep".xvg | awk '{print $1}')
		points_zero=$(grep '   0   ' "$output_directory"/hb_r"$grp"."$rep".xvg | wc -l | awk '{print $1}')
		
		hb_num=$(gmx analyze -f "$output_directory"/hb_r"$grp"."$rep".xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\n", $2}')
		hb_occ=$(echo $points_total $points_zero | awk '{printf "%.4f\n", 1.000000-($2/$1)}')
		
		#echo -e "$rep""\t""$grp""\t""$hb_num""\t""$hb_occ"
		
		echo "$hb_num" >> "$output_directory"/temporary_hb_r"$grp"_num.xvg
		echo "$hb_occ" >> "$output_directory"/temporary_hb_r"$grp"_occ.xvg
		
		let grp=$grp+1
		
	done
	
done





echo -e "residue""\t""hb_num_avg""\t""hb_num_err" > "$output_directory"/hb_r_num.xvg
echo -e "residue""\t""hb_occ_avg""\t""hb_occ_err" > "$output_directory"/hb_r_occ.xvg
echo "$total_grp"                                 > "$output_directory"/hb_r_num.dat
echo "$total_grp"                                 > "$output_directory"/hb_r_occ.dat


hb_num_avg_err=$(gmx analyze -f "$output_directory"/temporary_hb_num.xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\t%.4f\n", $2, $4}')
hb_occ_avg_err=$(gmx analyze -f "$output_directory"/temporary_hb_occ.xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\t%.4f\n", $2, $4}')

echo -e "$hb_num_avg_err" >> "$output_directory"/hb_num.xvg
echo -e "$hb_occ_avg_err" >> "$output_directory"/hb_occ.xvg


grp=1
while [ $grp -le $total_grp ]; do
	
	res_name=$(grep '\[ ' "$index_res"  | tail -n+2 | head -n $grp | tail -n 1 | awk -F ' ' '{print $2}' | awk -F '_' '{print $2}')
	res_num=$(grep  '\[ ' "$index_res"  | tail -n+2 | head -n $grp | tail -n 1 | awk -F ' ' '{print $2}' | awk -F '_' '{print $3}')
	
	hb_num_avg_err=$(gmx analyze -f "$output_directory"/temporary_hb_r"$grp"_num.xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\t%.4f\n", $2, $4}')
	hb_occ_avg_err=$(gmx analyze -f "$output_directory"/temporary_hb_r"$grp"_occ.xvg 2>/dev/null | grep 'SS1' | awk '{printf "%.4f\t%.4f\n", $2, $4}')
	
	echo -e "$res_name"'-'"$res_num""\t""$hb_num_avg_err" >> "$output_directory"/hb_r_num.xvg
	echo -e "$res_name"'-'"$res_num""\t""$hb_occ_avg_err" >> "$output_directory"/hb_r_occ.xvg
	
	let grp=$grp+1
done

rm "$output_directory"/temporary*
rm "$output_directory"/\#* 2>/dev/null

tail -n+2 "$output_directory"/hb_r_num.xvg | awk 'BEGIN{n=1}{print n++,$2}' >> "$output_directory"/hb_r_num.dat
tail -n+2 "$output_directory"/hb_r_occ.xvg | awk 'BEGIN{n=1}{print n++,$2}' >> "$output_directory"/hb_r_occ.dat

gmx editconf -f "$dcsign_renum" -bf "$output_directory"/hb_r_num.dat -o temporary.pdb >/dev/null 2>/dev/null
gmx editconf -f temporary.pdb -resnr 252 -o "$output_directory"/hb_r_num.pdb >/dev/null 2>/dev/null
rm temporary.pdb

gmx editconf -f "$dcsign_renum" -bf "$output_directory"/hb_r_occ.dat -o temporary.pdb >/dev/null 2>/dev/null
gmx editconf -f temporary.pdb -resnr 252 -o "$output_directory"/hb_r_occ.pdb >/dev/null 2>/dev/null
rm temporary.pdb

else

	echo -n 'no trajectories ... '

fi

