#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node=128
#SBATCH --job-name=s3S2

module unload intel-oneapi
module load intel-mpi

nodes=1
ntasks=128
cores=$((nodes * ntasks))

lmp=/ocean/projects/mca94p017p/kuntalg/softwares/LAMMPS/lammps-21Jul20/src/lmp_mpi

i=1
nbeads=3000
ncap=60444
rna_mol=1315
rna2_mol=1316
cap_bonds=681966

while [ $i -le $nbeads ]
do

	echo " "
	echo "Entering iteration (for cleavage) for bead $i:"
	echo " "	

#       Relabelling bond types for cleavage for both strands

	# FOR POLYMER 1	
	num_bonds=$(grep 'bonds' start.data | awk '{print $1}')

        bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
        angle_line=$(grep -n "Angles" start.data | cut -d: -f1)

	i_rna=$((i + ncap))
	ii=${i_rna}

	awk -v min="$bond_line" -v max="$angle_line" -v ii="$ii" 'NR>=min && NR<=max && $2==521 && $3==ii {$2=524} 1' start.data > temp.data
	mv temp.data start.data

	mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
        cat dump1.lammpstrj >> traj.lammpstrj

#       Checking if cleavage is successful (ie, RNP-DNA bonds have broken) else run again

	num_new_bonds=$(grep 'bonds' start.data | awk '{print $1}')

	nb=$((num_bonds - num_new_bonds))

	count=0
#       Continuing runs to ensure cleavage
        while [ $nb -ne 1 ]; do

                if [ ${count} -gt 100 ]; then
                        echo " "
			echo "Entering continued run for bead $i (for chain 1):"
                        echo " "
                fi

                mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
                cat dump1.lammpstrj >> traj.lammpstrj

		num_new_bonds=$(grep 'bonds' start.data | awk '{print $1}')

		nb=$((num_bonds - num_new_bonds))

                count=$((${count}+1))

        done

	i=$[$i+1]
done
