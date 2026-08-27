#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node=128
#SBATCH --job-name=s3SI

module unload intel-oneapi
module load intel-mpi

nodes=1
ntasks=128
cores=$((nodes * ntasks))

lmp=/ocean/projects/mca94p017p/kuntalg/softwares/LAMMPS/lammps-21Jul20/src/lmp_mpi

rm *lammpstrj* *restart* start.data log.lammps

i=1
nbeads=3000
ncap=60444
rna_mol=1315
rna2_mol=1316
cap_bonds=681966

cp cap_rnp_dntp_ang.data start.data

mpirun -np ${cores} ${lmp} -in inCapsid_init -var SEED ${RANDOM}

sed -i -e 's/^/ /' start.data   # Adds a space at the beginning of every line for sed purposes

while [ $i -le $nbeads ]
do

	echo " "
	echo "Entering iteration for bead $i:"
	echo " "	

#       Relabelling CG RNA beads so that dNTP can bind to it
        atom_line=$(grep -n "Atoms" start.data | cut -d: -f1)
        bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)

	i_rna=$((i + ncap))
#	i_rna2=$((i_rna + 3000))   # For res-id 1316: the second RNP polymer
        sed -i "${atom_line},${bond_line}s/ ${i_rna} ${rna_mol} 93 / ${i_rna} ${rna_mol} 97 /g" start.data
#        sed -i "${atom_line},${bond_line}s/ ${i_rna2} ${rna2_mol} 93 / ${i_rna2} ${rna2_mol} 97 /g" start.data

        mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
        cat dump1.lammpstrj >> traj.lammpstrj

#	Checking if dNTP binding is successful
        sed -i -e 's/^/ /' start.data
        wc=$(grep -i " ${i_rna} ${rna_mol} 95 " start.data | wc -l)
#        wc2=$(grep -i " ${i_rna2} ${rna2_mol} 95 " start.data | wc -l)
        sed -i 's/^[ \t]*//' start.data

	count=0
#	Continuing runs to ensure dNTP binds	
	while [ $wc -eq 0 ]; do

		if [ ${count} -gt 100 ]; then
			echo " "
			echo "Entering continued run for bead $i of chain 1:"
			echo " "
		fi
	
	        mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
        	cat dump1.lammpstrj >> traj.lammpstrj

                sed -i -e 's/^/ /' start.data
                wc=$(grep -i " ${i_rna} ${rna_mol} 95 " start.data | wc -l)
                sed -i 's/^[ \t]*//' start.data

		count=$((${count}+1))

	done

#	count=0
#	Continuing runs to ensure dNTP binds	
#	while [ $wc2 -eq 0 ]; do

#		if [ ${count} -gt 100 ]; then
#			echo " "
#			echo "Entering continued run for bead $i of chain 2:"
#			echo " "
#		fi
	
#		mpirun -np ${cores} ${lmp} -in inCapsid -var SEED ${RANDOM}
#		cat dump1.lammpstrj >> traj.lammpstrj

#                sed -i -e 's/^/ /' start.data
#                wc2=$(grep -i " ${i_rna2} ${rna2_mol} 95 " start.data | wc -l)
#                sed -i 's/^[ \t]*//' start.data

#		count=$((${count}+1))

#	done


#	dNTP-dNTP binding
        if [ $i -ge 2 ]; then
		bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
		angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
		num_bonds=$(sed -n 's/^\s*\([0-9]\+\)\s\+bonds$/\1/p' start.data)
		rna_bonds=$((num_bonds - cap_bonds))		

		fourth_column_values=()
	
		sed -i '$!N;/\n.*Angles/!P;D' start.data  # Removes blank line before "Angles" line

		# FOR POLYMER 1	
		ii=${i_rna}
		jj=$[$ii-1]

		# Use awk to process the file
		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==521 && $3==jj {print $4}' start.data)
		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==521 && $3==ii {print $4}' start.data)

		new_bond_no=$[${num_bonds}+1]

		if (( fourth_1 < fourth_2 )); then
		  sed -i "/Angles/i\\${new_bond_no} 522 ${fourth_1} ${fourth_2}\n" start.data
		else
		  sed -i "/Angles/i\\${new_bond_no} 522 ${fourth_2} ${fourth_1}\n" start.data
		fi
	
#		sed -i '$!N;/\n.*Angles/!P;D' start.data  # Removes blank line before "Angles" line
	
#		# FOR POLYMER 2
#		ii=${i_rna2}
#		jj=$[$ii-1]

#		# Use awk to process the file
#		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==521 && $3==jj {print $4}' start.data)
#		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==521 && $3==ii {print $4}' start.data)

#		new_bond_no=$[${new_bond_no}+1]

#		if (( fourth_1 < fourth_2 )); then
#		  sed -i "/Angles/i\\${new_bond_no} 522 ${fourth_1} ${fourth_2}\n" start.data
#		else
#		  sed -i "/Angles/i\\${new_bond_no} 522 ${fourth_2} ${fourth_1}\n" start.data
#		fi

		sed -i "s/${num_bonds} bonds/${new_bond_no} bonds/g" start.data
	fi

#       dNTP-dNTP-dNTP angular stiffness
        if [ $i -ge 3 ]; then
		bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
		angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
		num_angles=$(sed -n 's/^\s*\([0-9]\+\)\s\+angles$/\1/p' start.data)
		rna_angles=${num_angles}     # cap_angles = 0	

		fourth_column_values=()
	
		# FOR POLYMER 1	
		ii=${i_rna}
		jj=$[$ii-1]
		kk=$[$ii-2]

		# Use awk to process the file
		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v kk="$kk" 'NR>=min && NR<=max && $2==521 && $3==kk {print $4}' start.data)
		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==521 && $3==jj {print $4}' start.data)
		fourth_3=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==521 && $3==ii {print $4}' start.data)

		new_angle_no=$[${num_angles}+1]

		sed -i '$a\'"${new_angle_no} 4 ${fourth_1} ${fourth_2} ${fourth_3}" start.data


#		# FOR POLYMER 2
#		ii=${i_rna2}
#		jj=$[$ii-1]
#		kk=$[$ii-2]

#		# Use awk to process the file
#		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v kk="$kk" 'NR>=min && NR<=max && $2==521 && $3==kk {print $4}' start.data)
#		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==521 && $3==jj {print $4}' start.data)
#		fourth_3=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==521 && $3==ii {print $4}' start.data)

#		new_angle_no=$[${new_angle_no}+1]

#		sed -i '$a\'"${new_angle_no} 4 ${fourth_1} ${fourth_2} ${fourth_3}" start.data


		sed -i "s/${num_angles} angles/${new_angle_no} angles/g" start.data
	fi


#	Stiffening of the RNA-RNA bonded interaction due to dNTP addition
	dntp=$(grep -i " 0 96 " start.data | wc -l)
	bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
	angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
	last_line=$(cat start.data | wc -l)

#	if [ $dntp -ge 4 ]; then   # 2 for each polymer  # 4 for both
	if [ $dntp -ge 2 ]; then   # 2 for each polymer  # 4 for both


		sed -i 's/$/ /' start.data

		# FOR POLYMER 1
		ii=${i_rna}
		jj=$[$ii-1]
		sed -i "${bond_line},${angle_line}s/ 520 $jj $ii / 523 $jj $ii /g" start.data

#		# FOR POLYMER 2
#		ii=${i_rna2}
#		jj=$[$ii-1]
#		sed -i "${bond_line},${angle_line}s/ 520 $jj $ii / 523 $jj $ii /g" start.data

		sed -i 's/ $//' start.data

	fi

#	Stiffening of the RNA-RNA angular interaction due to dNTP addition
	dntp=$(grep -i " 0 96 " start.data | wc -l)
	angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
	last_line=$(cat start.data | wc -l)

#	if [ $dntp -ge 6 ]; then   # 3 for each polymer
	if [ $dntp -ge 3 ]; then   # 3 for each polymer

		# FOR POLYMER 1
		ii=${i_rna}
		jj=$[$ii-1]
		kk=$[$ii-2]
		sed -i "${angle_line},${last_line}s/ 1 $kk $jj $ii/ 2 $kk $jj $ii/g" start.data

#		# FOR POLYMER 2
#		ii=${i_rna2}
#		jj=$[$ii-1]
#		kk=$[$ii-2]
#		sed -i "${angle_line},${last_line}s/ 1 $kk $jj $ii/ 2 $kk $jj $ii/g" start.data


	fi

	sed -i -e 's/^/ /' start.data

	i=$[$i+1]
done
