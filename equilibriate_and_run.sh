#!/bin/bash

MDP="/home/finnl92/thesis/extension/free_energy/MDP"
home="$(pwd)"

for (( i=2; i<20; i++ ))
do

	lambda=$i
	
	# Make a separate dir for each lambda value
	mkdir -p lambda_$lambda
	cp run_job.sh verify_equilibriation.py lambda_$lambda
	cd lambda_$lambda

	# Energy Minimization
	gmx194 grompp -f $MDP/em_steep_$lambda.mdp -c $home/initial.gro -p $home/topol.top -o em_$lambda.tpr -maxwarn 1 -nobackup
	./run_job.sh gmx194 mdrun -deffnm em_$lambda -nt 16 -pin on -nobackup

	# NVT Equilibriation
	gmx194 grompp -f $MDP/nvt_$lambda.mdp -c em_$lambda.gro -r em_$lambda.gro -p $home/topol.top -o nvt_$lambda.tpr -maxwarn 2 -nobackup
	./run_job.sh gmx194 mdrun -deffnm nvt_$lambda -nt 16 -pin on -nobackup
	echo "Temperature" | gmx194 energy -f nvt_$lambda.edr -o temperature_$lambda.xvg -nobackup

	# NPT Equilibriation 1
	gmx194 grompp -f $MDP/npt1_$lambda.mdp -c nvt_$lambda.gro -r nvt_$lambda.gro -t nvt_$lambda.cpt -p $home/topol.top -o npt1_$lambda.tpr -maxwarn 2 -nobackup
	./run_job.sh gmx194 mdrun -deffnm npt1_$lambda -nt 16 -pin on -nobackup
	echo "Temperature" | gmx194 energy -f npt1_$lambda.edr -o temperature1_$lambda.xvg -nobackup
	echo "Pressure" | gmx194 energy -f npt1_$lambda.edr -o pressure1_$lambda.xvg -nobackup
	echo "Density" | gmx194 energy -f npt1_$lambda.edr -o density1_$lambda.xvg -nobackup

	# NPT Equilibriation 2
	gmx194 grompp -f $MDP/npt2_$lambda.mdp -c npt1_$lambda.gro -r npt1_$lambda.gro -t npt1_$lambda.cpt -p $home/topol.top -o npt2_$lambda.tpr -maxwarn 1 -nobackup
	./run_job.sh gmx194 mdrun -deffnm npt2_$lambda -nt 16 -pin on -nobackup
	echo "Temperature" | gmx194 energy -f npt2_$lambda.edr -o temperature2_$lambda.xvg -nobackup
	echo "Pressure" | gmx194 energy -f npt2_$lambda.edr -o pressure2_$lambda.xvg -nobackup
	echo "Density" | gmx194 energy -f npt2_$lambda.edr -o density2_$lambda.xvg -nobackup

	# Production run
	gmx194 grompp -f $MDP/md_$lambda.mdp -c npt2_$lambda.gro -r npt2_$lambda.gro -t npt2_$lambda.cpt -p $home/topol.top -o md_$lambda.tpr -maxwarn 1 -nobackup
	./run_job.sh gmx194 mdrun -deffnm md_$lambda -nt 16 -pin on -nobackup	

	python verify_equilibriation.py $lambda
	
	cd $home
    
done
