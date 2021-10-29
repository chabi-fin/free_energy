#!/bin/bash

# create job directory
CALCDIR=$(basename $(pwd))
FULLPATH=$(pwd)
mkdir -p /scratch/finnl92/tmp/$CALCDIR
TDIR=/scratch/finnl92/tmp/$CALCDIR

# copy input data to job directory
cp -r * $TDIR
cd $TDIR

lambda=$1
gmxs mdrun -deffnm md_$lambda -nt 16 -pin on -nobackup

# Move job to home directory
cp -r $TDIR/* $FULLPATH
cd $FULLPATH
rm -rf $TDIR

exit 0
