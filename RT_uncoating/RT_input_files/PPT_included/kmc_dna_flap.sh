#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node=128
#SBATCH --job-name=flap_1

module unload intel-oneapi
module load intel-mpi

nodes=1
ntasks=128
cores=$((nodes * ntasks))

lmp=/ocean/projects/mca94p017p/kuntalg/softwares/LAMMPS/lammps-21Jul20/src/lmp_mpi

rm *lammpstrj* *restart* log.lammps start.data

#i=1
global_count=100    #Not starting from bead #1 to mimic PPT (~2900 from the other end)
nbeads=1600
ncap=60444
rna_mol=1315
rna2_mol=1316
cap_bonds=681966

cp start_cdna.data start.data

job_id=$(for f in slurm-*.out; do   echo "$f" | grep -oP '(?<=slurm-)\d+(?=\.out)'; done)

mpirun -np ${cores} ${lmp} -in inCapsid_init -var SEED ${RANDOM}

sed -i -e 's/^/ /' start.data   # Adds a space at the beginning of every line for sed purposes

flap_count=1500

while [ ${global_count} -le $nbeads ]
do

	if (( ${flap_count} > 1533 )); then
	    i=${global_count}
	    route=A
	else
	    rand=$(awk 'BEGIN { srand(); print rand() }')
	    if (( $(echo "$rand > 0.4" | bc -l) )); then
		i=${global_count}
		route=A
	    else
		i=${flap_count}
		route=B
	    fi
	fi

	echo " "
	echo "Entering iteration for bead $i:"
	echo " "	

	if [ "$(grep -i -n 'missing' slurm* | wc -l)" -ne 0 ]; then
		scancel ${job_id}
		exit 1
	fi

	if [ "$(grep -i -n 'continued' slurm* | wc -l)" -ne 0 ]; then
		scancel ${job_id}
		exit 1
	fi


#       Relabelling CG RNA beads so that dNTP can bind to it
        atom_line=$(grep -n "Atoms" start.data | cut -d: -f1)
        bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)

#	cdna_ids.dat stores the atom ids of the cDNA beads (generated after stage I)
	i_cdna=$(head -n $i cdna_ids.dat | tail -n 1)

        sed -i "${atom_line},${bond_line}s/ ${i_cdna} 0 96 / ${i_cdna} 0 98 /g" start.data

	total_count=$((global_count - 100 + flap_count - 1500))
	if (( ${total_count} > 1100 )); then
	    input_Cap="inCapsid_high"
	else
	    input_Cap="inCapsid_low"
	fi

        mpirun -np ${cores} ${lmp} -in ${input_Cap} -var SEED ${RANDOM}
        cat dump1.lammpstrj >> traj.lammpstrj

#	Checking if dNTP binding is successful
        sed -i -e 's/^/ /' start.data
        wc=$(grep -i " ${i_cdna} 0 99 " start.data | wc -l)

        sed -i 's/^[ \t]*//' start.data

	count=0
#	Continuing runs to ensure dNTP binds	
	while [ $wc -eq 0 ]; do

		if [ ${count} -gt 100 ]; then
			echo " "
			echo "Entering continued run for bead $i of chain 1:"
			echo " "
		fi

		if [ "$(grep -i -n 'missing' slurm* | wc -l)" -ne 0 ]; then
			scancel ${job_id}
			exit 1
		fi

		if [ "$(grep -i -n 'continued' slurm* | wc -l)" -ne 0 ]; then
			scancel ${job_id}
			exit 1
		fi

		if (( $i > 1100 )); then
		    input_Cap="inCapsid_high"
		else
		    input_Cap="inCapsid_low"
		fi

		mpirun -np ${cores} ${lmp} -in ${input_Cap} -var SEED ${RANDOM}
#        	cat dump1.lammpstrj >> traj.lammpstrj

                sed -i -e 's/^/ /' start.data
                wc=$(grep -i " ${i_cdna} 0 99 " start.data | wc -l)
                sed -i 's/^[ \t]*//' start.data

		count=$((${count}+1))

	done


#	dNTP-dNTP binding
	if [[ "$route" == "A" && ${global_count} -ge 102 ]] || [[ "$route" == "B" && ${flap_count} -ge 1502 ]]; then
#	if [ ${global_count} -ge 1502 ]; then
		bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
		angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
		num_bonds=$(sed -n 's/^\s*\([0-9]\+\)\s\+bonds$/\1/p' start.data)
		rna_bonds=$((num_bonds - cap_bonds))		

		fourth_column_values=()
	
		sed -i '$!N;/\n.*Angles/!P;D' start.data  # Removes blank line before "Angles" line

		# FOR POLYMER 1	
		ii=${i_cdna}
		j=$[$i-1]
		jj=$(head -n $j cdna_ids.dat | tail -n 1)

		# Use awk to process the file
		# This additional check is necessary as the template strand for this code (ie, cDNA) is itself made from dNTPs

		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==525 && $3==jj {print $4}' start.data)
		[[ -z "${fourth_1}" ]] && fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==525 && $4==jj {print $3}' start.data)

		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==525 && $3==ii {print $4}' start.data)
		[[ -z "${fourth_2}" ]] && fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==525 && $4==ii {print $3}' start.data)

		new_bond_no=$[${num_bonds}+1]

		if (( fourth_1 < fourth_2 )); then
		  sed -i "/Angles/i\\${new_bond_no} 526 ${fourth_1} ${fourth_2}\n" start.data
		else
		  sed -i "/Angles/i\\${new_bond_no} 526 ${fourth_2} ${fourth_1}\n" start.data
		fi
	
		sed -i "s/${num_bonds} bonds/${new_bond_no} bonds/g" start.data
	fi

#       dNTP-dNTP-dNTP angular stiffness
	if [[ "$route" == "A" && ${global_count} -ge 103 ]] || [[ "$route" == "B" && ${flap_count} -ge 1503 ]]; then
#	if [ ${global_count} -ge 1503 ]; then
		bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
		angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
		num_angles=$(sed -n 's/^\s*\([0-9]\+\)\s\+angles$/\1/p' start.data)
		rna_angles=${num_angles}     # cap_angles = 0	

		fourth_column_values=()
	
		# FOR POLYMER 1	
		ii=${i_cdna}
		j=$[$i-1]
		jj=$(head -n $j cdna_ids.dat | tail -n 1)
		k=$[$i-2]
		kk=$(head -n $k cdna_ids.dat | tail -n 1)

		# Use awk to process the file
		fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v kk="$kk" 'NR>=min && NR<=max && $2==525 && $3==kk {print $4}' start.data)
		[[ -z "${fourth_1}" ]] && fourth_1=$(awk -v min="${bond_line}" -v max="${angle_line}" -v kk="$kk" 'NR>=min && NR<=max && $2==525 && $4==kk {print $3}' start.data)

		fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==525 && $3==jj {print $4}' start.data)
		[[ -z "${fourth_2}" ]] && fourth_2=$(awk -v min="${bond_line}" -v max="${angle_line}" -v jj="$jj" 'NR>=min && NR<=max && $2==525 && $4==jj {print $3}' start.data)

		fourth_3=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==525 && $3==ii {print $4}' start.data)
		[[ -z "${fourth_3}" ]] && fourth_3=$(awk -v min="${bond_line}" -v max="${angle_line}" -v ii="$ii" 'NR>=min && NR<=max && $2==525 && $4==ii {print $3}' start.data)

		new_angle_no=$[${num_angles}+1]

		sed -i '$a\'"${new_angle_no} 6 ${fourth_1} ${fourth_2} ${fourth_3}" start.data

		sed -i "s/${num_angles} angles/${new_angle_no} angles/g" start.data
	fi


#	Stiffening of the RNA-RNA bonded interaction due to dNTP addition
	dntp=$(grep -i " 0 100 " start.data | wc -l)
	bond_line=$(grep -n "Bonds" start.data | cut -d: -f1)
	angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
	last_line=$(cat start.data | wc -l)

	if [ $dntp -ge 2 ]; then   # 2 for each polymer  # 4 for both


		sed -i 's/$/ /' start.data

		# FOR POLYMER 1
		ii=${i_cdna}
		j=$[$i-1]
		jj=$(head -n $j cdna_ids.dat | tail -n 1)
		sed -i "${bond_line},${angle_line}s/ 522 $jj $ii / 527 $jj $ii /g" start.data
		sed -i "${bond_line},${angle_line}s/ 522 $ii $jj / 527 $ii $jj /g" start.data

		sed -i 's/ $//' start.data

	fi

#	Stiffening of the RNA-RNA angular interaction due to dNTP addition
	dntp=$(grep -i " 0 100 " start.data | wc -l)
	angle_line=$(grep -n "Angles" start.data | cut -d: -f1)
	last_line=$(cat start.data | wc -l)

	if [ $dntp -ge 3 ]; then   # 3 for each polymer

		# FOR POLYMER 1
		ii=${i_cdna}
		j=$[$i-1]
		jj=$(head -n $j cdna_ids.dat | tail -n 1)
		k=$[$i-2]
		kk=$(head -n $k cdna_ids.dat | tail -n 1)

		sed -i "${angle_line},${last_line}s/ 4 $kk $jj $ii/ 7 $kk $jj $ii/g" start.data
		sed -i "${angle_line},${last_line}s/ 4 $ii $jj $kk/ 7 $ii $jj $kk/g" start.data

	fi

	sed -i -e 's/^/ /' start.data

	[ "$route" = "B" ] && flap_count=$[${flap_count}+1]

	[ "$route" = "A" ] && global_count=$[${global_count}+1]
done
