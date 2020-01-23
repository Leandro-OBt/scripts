#!/bin/bash


#pdb2gmx -f "new.pdb" -o "new.gro" -p "new.top" -ignh -ff "amber99sb-ildn" -water "none" >/dev/null 2>/dev/null 
gmx mdrun -s "md.tpr" -rerun "new.pdb" -deffnm "new" -nt 1 >/dev/null 2>/dev/null 
echo "Potential" | gmx energy -f "new.edr" -o "new.xvg" >/dev/null 2>/dev/null 
tail -n 1 "new.xvg" | awk '{print $2}' > "potential.dat"
rm "new."*


