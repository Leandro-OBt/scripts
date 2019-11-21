# Script para execução de MM/PBSA com AmberTools a partir de trajetórias do gromacs
#
#
# 04 Oct 2016


input_xtc=$1			# Trajectory containing ONLY the frames that will be used for the MM/PBSA calculation
input_tpr=$2			# Corresponding tpr
input_ndx=$3			# Index file. Must contain grp_receptor and grp_ligand, e.g. Protein-H_chain1, Other_chain1
grp_receptor=$4		# For proteins, the group must have no hydrogens, e.g. Protein-H
grp_ligand=$5			# grp name for the ligand
prefix=$6				# For temporary and output files
internal_dielectric=$7	# internal dielectric constant
input_charges_mol2=$8	# .mol2 with the correct charges (only for gaff)
input_frcmod=$9		# .frcmod file (only for gaff)

path_gmx='/home/'"$USER"'/programs/gmx-5.1.2/bin/'
path_gmx='/usr/local/bin/'
path_ambertools='/home/'"$USER"'/programs/amber16/bin'
path_ambertools='/home/'"$USER"'/programs/amber18/bin'

output_filename_energies="$prefix""_mmpbsa_energies.xvg"
output_filename_statistics="$prefix""_mmpbsa_results.dat"
output_filename_decomp="$prefix""_mmpbsa_decomp.dat"

ligand_protein=0	# currently, if this is =1, glycam should be = 1
glycam=0	# If =1, GLYCAM_06j-1. Otherwise, uses GAFF
	branched_ligand=0
		# A1, A2, B1, B2
		#branch_total=1
		#branch_resn_all="0fA"
		#branch_resi_all=4
		#branch_atm_all="H2O"
		# core2
		#branch_total=1
		#branch_resn_all="0LB"
		#branch_resi_all=3
		#branch_atm_all="H2O"
		# core4
		#branch_total=1
		#branch_resn_all="0LB"
		#branch_resi_all=4
		#branch_atm_all="H2O"
		
		# GlcNAc2Man5
		branch_total=4
		branch_resn_all="VMB VMA 0MA 0MA"
		branch_resi_all="4 5 6 7"
		branch_atm_all="O3 O3 H2O H2O"
		
#options

save_pdb_frames=0
create_amber_files=0
create_mmpbsa_script=1
run_calculations=1
	gbsa=0	# If =1, runs MM/GBSA instead of MM/PBSA
statistics=0
	bootstrap_resamplings=50
	keep_raw_files=1
	
# dumb parallelization (NOT IMPLEMENTED YET)
total_cores=6



echo
echo -n "Starting MM/"
if  [ $gbsa -eq 0 ]; then 
	echo -n "PB"
else
	echo -n "GB"
fi
echo "SA calculations ..."


if [ $save_pdb_frames -eq 1 ]; then
	time_begin=$(date +%s)
	
	# First, save all frames separately as pdb
	# !!!!!!!!!! The frame must contain the complex. Proteins must have no hydrogens!!!!
	echo -en "\t""Saving frames in separate pdb files ... "
	echo "$grp_receptor" | "$path_gmx"/gmx trjconv -f "$input_xtc" -s "$input_tpr" -n "$input_ndx" -o "$prefix"_temporary_rec_frame..pdb -sep >/dev/null 2>/dev/null
	echo "$grp_ligand" | "$path_gmx"/gmx trjconv -f "$input_xtc" -s "$input_tpr" -n "$input_ndx" -o "$prefix"_temporary_lig_frame..pdb -sep >/dev/null 2>/dev/null

	total_frames=$(find . -maxdepth 1 -name "$prefix""_temporary_rec_frame.*.pdb" | wc -l | awk '{print $1}')
	frame=0

	rm "$prefix"_complex.*.pdb 2>/dev/null
	while [ $frame -lt $total_frames ]; do
		# Firts, uses the temporary_rec_frame files to look for chains and correctly put "TER" 
		sed 's\OC1\O  \g' "$prefix"_temporary_rec_frame."$frame".pdb | sed 's\OC2\OXT\g' | sed 's\CD  ILE\CD1 ILE\g' | sed 's\HETATM\ATOM  \g' | grep "ATOM" > "$prefix"_temporary_frame.pdb
		awk -v ch=1 -v prefix="$prefix" '{print > prefix"_temporary_chain_"ch".pdb"}/OXT/{ch++;}' "$prefix"_temporary_frame.pdb 
		rm "$prefix"_temporary_frame.pdb
		
		total_chains=$(find . -maxdepth 1 -name "$prefix""_temporary_chain_*.pdb" | wc -l | awk '{print $1}')
		chain=1
		while [ $chain -le $total_chains ]; do
			cat "$prefix"_temporary_chain_"$chain".pdb >> "$prefix"_rec."$frame".pdb
			echo "TER" >> "$prefix"_rec."$frame".pdb
			
			let chain=$chain+1
		done
		echo "END" >> "$prefix"_rec."$frame".pdb
		
		if [ $ligand_protein -eq 1 ]; then
			
			sed 's\OC1\O  \g' "$prefix"_temporary_lig_frame."$frame".pdb | sed 's\OC2\OXT\g' | sed 's\CD  ILE\CD1 ILE\g' | sed 's\HETATM\ATOM  \g' | grep "ATOM" > "$prefix"_temporary_frame.pdb
			awk -v ch=1 -v prefix="$prefix" '{print > prefix"_temporary_chain_"ch".pdb"}/OXT/{ch++;}' "$prefix"_temporary_frame.pdb 
			rm "$prefix"_temporary_frame.pdb
			
			rm "$prefix"_lig."$frame".pdb 2>/dev/null
			total_chains=$(find . -maxdepth 1 -name "$prefix""_temporary_chain_*.pdb" | wc -l | awk '{print $1}')
			chain=1
			while [ $chain -le $total_chains ]; do
				cat "$prefix"_temporary_chain_"$chain".pdb >> "$prefix"_lig."$frame".pdb
				echo "TER" >> "$prefix"_lig."$frame".pdb
			
				let chain=$chain+1
			done
			echo "END" >> "$prefix"_lig."$frame".pdb
			
		else
			grep "ATOM" "$prefix"_temporary_lig_frame."$frame".pdb > "$prefix"_temporary_lig_frame_.pdb
			if [ $branched_ligand -eq 1 ]; then
				
				branch=1
				while [ $branch -le $branch_total ]; do
				
					branch_resn=$(echo "$branch_resn_all" | awk '{print $'"$branch"'}')
					branch_resi=$(echo "$branch_resi_all" | awk '{print $'"$branch"'}')
					branch_atm=$(echo  "$branch_atm_all"  | awk '{print $'"$branch"'}')
					
					line=$(grep "$branch_resn" "$prefix"_temporary_lig_frame_.pdb | grep "$branch_atm" | grep " ""$branch_resi"" ")
					line_num=$(grep "$line" "$prefix"_temporary_lig_frame_.pdb -n | sed 's\:\ \g' | tail -n 1 | awk '{print $1}')
					
					head -n $line_num "$prefix"_temporary_lig_frame_.pdb > "$prefix"_temporary_lig_frame_ter.pdb
					echo "TER" >> "$prefix"_temporary_lig_frame_ter.pdb
					let line_num=$line_num+1
					
					tail -n+$line_num "$prefix"_temporary_lig_frame_.pdb >> "$prefix"_temporary_lig_frame_ter.pdb
					
					mv "$prefix"_temporary_lig_frame_ter.pdb "$prefix"_temporary_lig_frame_.pdb
					
					let branch=$branch+1
				done
			fi
			
			echo "TER" >> "$prefix"_temporary_lig_frame_.pdb
			cp "$prefix"_temporary_lig_frame_.pdb "$prefix"_lig."$frame".pdb
			echo "END" >> "$prefix"_lig."$frame".pdb
			
			if [ $glycam -eq 0 ]; then
				antechamber -i "$prefix"_lig."$frame".pdb -fi pdb -o "$prefix"_temporary_lig.mol2 -fo mol2 -pf yes >/dev/null 2>/dev/null
				antechamber -i "$prefix"_temporary_lig.mol2 -fi mol2 -o "$prefix"_lig."$frame".mol2 -fo mol2 -ao crg -a "$input_charges_mol2" -fa mol2 -pf yes -dr no >/dev/null 2>/dev/null
				rm "$prefix"_lig."$frame".pdb
			fi
		fi
		
		
		rm "$prefix"_temporary_chain_*.pdb "$prefix"_temporary_lig_frame_.pdb "$prefix"_temporary_lig.mol2 2>/dev/null
		
		let frame=$frame+1
	done
	rm "$prefix"_temporary_lig_frame.*.pdb "$prefix"_temporary_rec_frame.*.pdb
	echo "done"
	
	time_end=$(date +%s)
	let time=$time_end-$time_begin
	echo -e "pdb frames""\t""$time" >> "$prefix"_time.log
fi



if [ $create_amber_files -eq 1 ]; then
	time_begin=$(date +%s)
	
	# Build tleap scripts
	echo -en "\t""Creating AMBER files ... "
	
	total_frames=$(find . -maxdepth 1 -name "$prefix""_rec.*.pdb" | wc -l | awk '{print $1}')
	frame=0 
	while [ $frame -lt $total_frames ]; do
		echo 'source oldff/leaprc.ff99SBildn' > "$prefix"_temporary_tleap_script
		if [ $glycam -eq 1 ]; then
			echo 'source leaprc.GLYCAM_06j-1' >> "$prefix"_temporary_tleap_script
		else
			echo 'source leaprc.gaff' >> "$prefix"_temporary_tleap_script
			if [ ! -z "$input_frcmod" ]; then
				echo 'loadamberparams '"$input_frcmod" >> "$prefix"_temporary_tleap_script
			fi
		fi
		echo 'set default PBRadii mbondi2' >> "$prefix"_temporary_tleap_script
		
		if [ $glycam -eq 1 ] || [ $ligand_protein -eq 1 ] ; then
			echo 'lig = loadpdb '"$prefix"'_lig.'"$frame"'.pdb' >> "$prefix"_temporary_tleap_script
		else
			echo 'lig = loadmol2 '"$prefix"'_lig.'"$frame"'.mol2' >> "$prefix"_temporary_tleap_script
		fi
		echo 'saveamberparm lig '"$prefix"'_lig.prmtop '"$prefix"'_lig.'"$frame"'.inpcrd' >> "$prefix"_temporary_tleap_script
		echo 'rec = loadpdb '"$prefix"'_rec.'"$frame"'.pdb' >> "$prefix"_temporary_tleap_script
		echo 'saveamberparm rec '"$prefix"'_rec.prmtop '"$prefix"'_rec.'"$frame"'.inpcrd' >> "$prefix"_temporary_tleap_script
		echo 'complex = combine { rec lig }' >> "$prefix"_temporary_tleap_script
		echo 'saveamberparm complex '"$prefix"'_complex.prmtop '"$prefix"'_complex.'"$frame"'.inpcrd' >> "$prefix"_temporary_tleap_script
		echo 'quit' >> "$prefix"_temporary_tleap_script
		
		"$path_ambertools"/tleap -f "$prefix"_temporary_tleap_script >/dev/null 2>/dev/null
		
#		rm "$prefix"_rec."$frame".pdb
#		if [ $glycam -eq 1 ]; then
#			rm "$prefix"_lig."$frame".pdb
#		else
#			rm "$prefix"_lig."$frame".mol2
#		fi
		
		let frame=$frame+1
	done
	echo "done"
	
	rm "$prefix"_temporary_tleap_script leap.log
	
	time_end=$(date +%s)
	let time=$time_end-$time_begin
	echo -e "Amber files""\t""$time" >> "$prefix"_time.log
fi





if [ $create_mmpbsa_script -eq 1 ]; then
	time_begin=$(date +%s)
	
	# Build the mmpbsa input file
	echo -en "\t""Creating the input script ... "
	echo 'Input file for running MM/PBSA with AmberTools' > "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '&general'                                      >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '   entropy = 0,'                               >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '   keep_files = 0,'                            >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '   verbose = 1,'                               >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '/'                                             >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	if [ $gbsa -eq 1 ]; then
		echo '&gb'                                   >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   igb=5,'                             >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   saltcon=0.150,'                     >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '/'                                     >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	else
		echo '&pb'                                   >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   exdi = 80.0,'                       >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   indi = '"$internal_dielectric"','   >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   istrng=0.150, '                     >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   radiopt=0,'                         >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '   inp=1,'                             >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
		echo '/'                                     >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	fi
	echo '&decomp'                                       >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '  idecomp=1,'                                  >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '  dec_verbose=0,'                              >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo '/'                                             >> "$prefix"_mmpbsa_eint-"$internal_dielectric".in
	echo "done"
	
	time_end=$(date +%s)
	let time=$time_end-$time_begin
	echo -e "Script""\t""$time" >> "$prefix"_time.log
fi




if [ $run_calculations -eq 1 ]; then
	time_begin=$(date +%s)
	
	total_frames=$(find . -maxdepth 1 -name "$prefix""_complex.*.inpcrd" | wc -l | awk '{print $1}')
	echo -en "\t""Running calculations for ""$total_frames"" frames ... "
	
	frame=0 
	while [ $frame -lt $total_frames ]; do
		
		"$path_ambertools"/MMPBSA.py -i "$prefix"_mmpbsa_eint-"$internal_dielectric".in -cp "$prefix"_complex.prmtop -rp "$prefix"_rec.prmtop -lp "$prefix"_lig.prmtop -y "$prefix"_complex."$frame".inpcrd -prefix "$prefix""_temporary_mmpbsa_eint-""$internal_dielectric""_" -eo "$prefix""_eint-""$internal_dielectric""_mmpbsa_energies.""$frame"".csv" -deo "$prefix""_eint-""$internal_dielectric""_mmpbsa_decomp.""$frame"".csv" -o "$prefix""_eint-""$internal_dielectric""_FINAL_RESULTS_MMPBSA.""$frame"".dat" -do  "$prefix""_eint-""$internal_dielectric""_FINAL_DECOMP_MMPBSA.""$frame"".dat" >/dev/null 2>/dev/null
		
#		rm "$prefix"_complex."$frame".inpcrd "$prefix"_FINAL_DECOMP_MMPBSA.dat "$prefix"_FINAL_RESULTS_MMPBSA.dat
		
		let frame=$frame+1
	done
	echo "done"

#	rm "$prefix"_mmpbsa.in "$prefix"_rec.prmtop "$prefix"_lig.prmtop "$prefix"_complex.prmtop
	
	time_end=$(date +%s)
	let time=$time_end-$time_begin
	echo -e "Calculations""\t""$time" >> "$prefix"_time.log
fi








if [ $statistics -eq 1 ]; then
	time_begin=$(date +%s)
	

	# Bootstrap Statistics with R
	# to install R: sudo apt-get install r-base
	echo -en "\t""Bootstrap statistics with ""$bootstrap_resamplings"" resamplings ... "
	
	
	echo '@ s0 legend "VdW"' > "$output_filename_energies"
	echo '@ s1 legend "Elec"' >> "$output_filename_energies"
	echo '@ s2 legend "PB"' >> "$output_filename_energies"
	echo '@ s3 legend "SA"' >> "$output_filename_energies"
	echo '@ s4 legend "MM"' >> "$output_filename_energies"
	echo '@ s5 legend "PBSA"' >> "$output_filename_energies"
	echo '@ s6 legend "MM/PBSA"' >> "$output_filename_energies"
	
	
	
	# Prepare input files
	
	first_frame=$(find . -maxdepth 1 -name "$prefix""_mmpbsa_decomp.*.csv" | sed 's/\./\ /g' | awk '{print $(NF-1)}' | sort -n | head -n 1)
	total_residues=$(cat "$prefix"_mmpbsa_decomp."$first_frame".csv | grep ^'1,' | grep ',R ' | wc -l | awk '{print $1}')
	
	rm "$prefix"_temporary_decomp_res_* 2>/dev/null
	
	for frame in $(find . -maxdepth 1 -name "$prefix""_mmpbsa_decomp.*.csv" | sed 's/\./\ /g' | awk '{print $(NF-1)}' | sort -n); do
		
		grep '0,' "$prefix"_mmpbsa_energies."$frame".csv | tail -n 1 | sed 's\,\ \g' | awk -v frame=$frame '{printf "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", frame,$2,$3,$4,$5,$2+$3,$4+$5,$2+$3+$4+$5}' >> "$output_filename_energies"
		
		grep ^'1' "$prefix"_mmpbsa_decomp."$frame".csv | grep ',R ' | sed 's\,\ \g' | awk '{print $2"-"$3,$7,$8,$9,$10,$11,$12}' > "$prefix"_temporary_decomp."$frame".csv
		grep ^'1' "$prefix"_mmpbsa_decomp."$frame".csv | grep ',L ' | sed 's\,\ \g' | awk '{print $2"-"$3,$7,$8,$9,$10,$11,$12}' >> "$prefix"_temporary_decomp."$frame".csv
		
		
		
		# Split the contribution files into residue-specific files
		
		residue=1
		while [ $residue -le $total_residues ]; do
			
			head -n $residue "$prefix"_temporary_decomp."$frame".csv | tail -n 1 | awk -v frame=$frame '{printf "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", frame,$2,$3,$4,$5,$6,$2+$3+$4+$5+$6}' >> "$prefix"_temporary_decomp_res_"$residue"
			
			let residue=$residue+1
		done
		
	done



	# Statistics for the global values
	
	
	rm "$prefix"_temporary_values 2>/dev/null
	i=2
	while [ $i -le 8 ]; do
		grep -v "@" "$output_filename_energies" | awk '{print $'"$i"'}' > "$prefix"_data.dat

		echo 'list <- read.table("'"$prefix"'_data.dat")' > "$prefix"_temporary_script_R.r
		echo 'data <- as.numeric( unlist(list) )' >> "$prefix"_temporary_script_R.r
		echo 'resamples <- lapply(1:'"$bootstrap_resamplings"', function(i) sample(data, replace = T))' >> "$prefix"_temporary_script_R.r
		echo 'averages <- sapply(resamples, mean)' >> "$prefix"_temporary_script_R.r
		echo 'write( mean(data), file="'"$prefix"'_avg.dat" )' >> "$prefix"_temporary_script_R.r
		echo 'write( sqrt( var(averages) ) , file="'"$prefix"'_stderr.dat" )' >> "$prefix"_temporary_script_R.r
		echo 'q(save="no")' >> "$prefix"_temporary_script_R.r
		
		rm .Rdata 2>/dev/null
		Rscript --vanilla "$prefix"_temporary_script_R.r
		rm .Rdata 2>/dev/null
		
		avg=$(cat "$prefix"_avg.dat)
		stderr=$(cat "$prefix"_stderr.dat)
		
		echo "$avg"" ""$stderr" | awk '{printf "%.4f\t%.4f\n", $1,$2}' >> "$prefix"_temporary_values
		
		let i=$i+1
	done

	echo -e "vdw\t\nelec\t\npb\t\nsa\t\nmm\t\npbsa\t\nmmpbsa\t" > "$prefix"_temporary_header

	echo -e "contrib\t\tavg[kcal/mol]\tstderr[kcal/mol]" > "$output_filename_statistics"
	paste "$prefix"_temporary_header "$prefix"_temporary_values >> "$output_filename_statistics"
	rm "$prefix"_temporary_header "$prefix"_temporary_values




	# Statistics for the residue-specific values


	echo -e "residue\tinternal\tvdw[kcal/mol]\telec[kcal/mol]\tpb[kcal/mol]\tsa[kcal/mol]\ttotal[kcal/mol]" > "$output_filename_decomp"

	residue=1
	while [ $residue -le $total_residues ]; do
		
		id=$(head -n $residue "$prefix"_temporary_decomp."$first_frame".csv | tail -n 1 | awk '{print $1}')
		echo -en "$id""\t" >> "$output_filename_decomp"
		
		i=2
		while [ $i -le 7 ]; do
			awk '{print $'"$i"'}' "$prefix"_temporary_decomp_res_"$residue" > "$prefix"_data.dat
			
			echo 'list <- read.table("'"$prefix"'_data.dat")' > "$prefix"_temporary_script_R.r
			echo 'data <- as.numeric( unlist(list) )' >> "$prefix"_temporary_script_R.r
			echo 'resamples <- lapply(1:'"$bootstrap_resamplings"', function(i) sample(data, replace = T))' >> "$prefix"_temporary_script_R.r
			echo 'averages <- sapply(resamples, mean)' >> "$prefix"_temporary_script_R.r
			echo 'write( mean(data), file="'"$prefix"'_avg.dat" )' >> "$prefix"_temporary_script_R.r
			echo 'write( sqrt( var(averages) ) , file="'"$prefix"'_stderr.dat" )' >> "$prefix"_temporary_script_R.r
			
			rm .Rdata 2>/dev/null
			#R CMD BATCH
			Rscript --vanilla "$prefix"_temporary_script_R.r
			rm .Rdata 2>/dev/null
			
			avg=$(cat "$prefix"_avg.dat)
			stderr=$(cat "$prefix"_stderr.dat)
			
			echo "$avg"" ""$stderr""\t" | awk '{printf "%.4f\t%.4f\t", $1,$2}' >> "$output_filename_decomp"
			
			let i=$i+1
		done
		echo ""  >> "$output_filename_decomp"
		
		let residue=$residue+1
	done

	rm "$prefix"_data.dat "$prefix"_avg.dat "$prefix"_stderr.dat "$prefix"_temporary_script_R.r "$prefix"_temporary_script_R.r.Rout  "$prefix"_temporary_decomp_res_* "$prefix"_temporary_decomp.*.csv
	rm leap.log 2>/dev/null
	
	
	if [ $keep_raw_files -eq 0 ]; then
		rm "$prefix"_mmpbsa_energies.*.csv "$prefix"_mmpbsa_decomp.*.csv
	fi

	echo "done"
	echo
	
	time_end=$(date +%s)
	let time=$time_end-$time_begin
	echo -e "Statistics""\t""$time" >> "$prefix"_time.log

fi




