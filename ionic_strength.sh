#!/bin/bash

#
#
#	Script for calculating the number of cations and anions necessary to reach a given ionic strenght
# for a given number of water molecules in the simulation box
#
#		Leandro Oliveira Bortot
#		04 Oct 2016
#

input_gro=$1	# gro file to which the ions will be added
is=$2		# final ionic strenght
z_cat=$3		# absolute charge of the cation
z_an=$4		# absolute charge of the anion
lipid14=$5	# if =1, use amber names for water
if [ -z $lipid14 ]; then
	lipid14=0
fi

# Default values is they are left blank
if [ -z $is ]; then is=0.150 ; fi
if [ -z $z_cat ]; then z_cat=1 ; fi
if [ -z $z_an  ]; then z_an=1  ; fi

if [ $lipid14 -eq 1 ]; then
	sol=$(grep "WAT " "$input_gro" | grep ' O' | wc -l)
else
	sol=$(grep " OW" "$input_gro" | grep "SOL " | wc -l)
fi

c_cat=$(echo "(2*$is)/(($z_cat*$z_cat)+($z_an*$z_cat) )" | bc -l | awk '{printf "%.3f", $1}')
c_an=$(echo "$c_cat*($z_cat/$z_an)" | bc -l | awk '{printf "%.3f", $1}')

x_cat=$(echo "$c_cat*($sol/55.555)" | bc -l | awk '{printf "%d", $1}')
x_an=$(echo  "$c_an*($sol/55.555)"  | bc -l | awk '{printf "%d", $1}')


echo -e "$x_cat\t$x_an"

