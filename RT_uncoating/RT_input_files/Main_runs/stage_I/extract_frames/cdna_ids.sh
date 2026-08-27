#!/bin/bash

i=1
nbeads=3000
ncap=60444

[ -f cdna_ids.dat ] && rm cdna_ids.dat
touch cdna_ids.dat

# This make a list of the cDNA atom (bead) IDs
# Input is the last LAMMPS data file from stage I (ie. formation of RNA-DNA hybrid)

while [ $i -le $nbeads ]
do

        echo "Entering iteration for bead $i:"

        bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
        angle_line=$(grep -n "Angles" start.data | cut -d: -f1)

	i_rna=$((i + ncap))

	var=$(awk "NR>$bond_line && NR<$angle_line && / 521 ${i_rna} / {print \$4}" start.data)
	echo ${var} >> cdna_ids.dat

	i=$[$i+1]

done
