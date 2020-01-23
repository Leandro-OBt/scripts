#!/bin/bash

echo "C-alpha Protein" | gmx rms -f "new.pdb" -s "old.pdb" -o "rmsd_ca-ca_previous.xvg" >/dev/null 2>/dev/null 
tail -n 1 "rmsd_ca-ca_previous.xvg" | awk '{print $2*10}' > "rmsd_ca-ca_previous.dat"
rm "rmsd_ca-ca_previous.xvg" \#* "new.pdb" "old.pdb" 2>/dev/null