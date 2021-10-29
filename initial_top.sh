#!/bin/bash

# Default solvent is TIP3P
solvent="spc216"

# Parse the command line
while getopts n:s: flag
do
    case "${flag}" in 
        n) name=${OPTARG};;
	s) solvent=${OPTARG};;
    esac
done

# If name.pdb not in directrory, exit
usage() { echo "Usage: ./initial_top.sh -n NAME -s solvent"; exit 1;}

if [ -z "${name}.pdb" ] || [ -z "$name" ]; then
    usage
fi

# Set up forcefield files
cp -r ~/thesis/ff/amber14sb_sp.ff .
cp amber14sb_sp.ff/residuetypes.dat . 

# Generate topology using .pdb of the molecule
gmx194 pdb2gmx -ff amber14sb_sp -f ${name}.pdb -o ${name}.gro -water tip3p -nobackup

# Define box with ligand at center, > 1.2 nm to edge
gmx194 editconf -f ${name}.gro -o ${name}_box.gro -c -d 1.2 -bt cubic -nobackup

# Solvate (select option "-s dcm.gro" for DCM, the default is TIP3P)
gmx194 solvate -cp ${name}_box.gro -cs $solvent -o ${name}_solv.gro -p topol.top -nobackup	

# Add sodium ions
#gmx194 grompp -f em_steep.mdp -c ${name}_solv.gro -p topol.top -o ions.tpr -maxwarn 1 -nobackup
#echo "SOL" | gmx194 genion -s ions.tpr -o initial.gro -p topol.top -np 1 -pname NA -pq 1 -nobackup

cp ${name}_solv.gro initial.gro

