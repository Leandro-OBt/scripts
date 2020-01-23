#!/bin/bash


t_b=100000
t_e=200000
t_dt=100

#t_b=5000
#t_e=10000
#t_dt=1

#t_b=50000
#t_e=100000
#t_dt=50

r=1
while [ $r -le 5 ]; do
#	echo 'Solute' | gmx covar -f ../md"$r"_non-water.xtc -s ../md"$r"_non-water.tpr -n ../index.ndx -o eigenval.$r.xvg -v eigenvec.$r.trr -av average.$r.pdb -l covar.$r.log -last -1 -nofit -noref -b $t_b -e $t_e -dt $t_dt
#	gmx anaeig -v eigenvec.$r.trr -temp 309 -nevskip 6 -entropy 2>/dev/null | awk '{print $9,$10"."$11,$6}' > S_conf.$r.xvg

	let r=$r+1
done


for method in Quasiharmonic Schlitter ; do
	cat S_conf.*.xvg | grep "$method" > tmp.xvg
	gmx analyze -f tmp.xvg 2>/dev/null | grep 'SS1' | awk '{printf "%f\t%f\tJ/(mol.K)\n",$2,$4}' > 'Sconf_avg-err_'"$method"'.dat'
	rm tmp.xvg
done
