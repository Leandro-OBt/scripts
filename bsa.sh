
r=1
while [ $r -le 5 ]; do
gmx sasa -f ../../md_non-water.$r.xtc -s ../../md_non-water.$r.pdb -n index.ndx -xvg none -dt 100 -b 25000 -e 100000 -o asa_chA.$r.xvg  -surface 'chA'     -output 'chA'
gmx sasa -f ../../md_non-water.$r.xtc -s ../../md_non-water.$r.pdb -n index.ndx -xvg none -dt 100 -b 25000 -e 100000 -o asa_chB.$r.xvg  -surface 'chB'     -output 'chB'
gmx sasa -f ../../md_non-water.$r.xtc -s ../../md_non-water.$r.pdb -n index.ndx -xvg none -dt 100 -b 25000 -e 100000 -o asa_chAB.$r.xvg -surface 'Protein' -output 'Protein'
paste asa_chA.$r.xvg asa_chB.1.xvg asa_chAB.$r.xvg | awk '{print $1,($3+$6)-$9}' > bsa_chA-chB.$r.xvg
let r=$r+1
done


