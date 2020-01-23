#!/bin/bash

echo "C-alpha C-alpha" | gmx rms -f "new.pdb" -s "1l2y_native.pdb" -o "rmsd_ca-ca_native.xvg" >/dev/null 2>/dev/null 
tail -n 1 "rmsd_ca-ca_native.xvg" | awk '{print $2*10}' > "rmsd_ca-ca_native.dat"
rm "rmsd_ca-ca_native.xvg" \#* "new.pdb" 2>/dev/null

