#!/bin/bash

# create job directory
CALCDIR=$(basename $(pwd))
FULLPATH=$(pwd)
mkdir -p /scratch/finnl92/tmp/$CALCDIR
TDIR=/scratch/finnl92/tmp/$CALCDIR

# copy input data to job directory
cp -r * $TDIR
cd $TDIR

#JOB INPUT GOES HERE
$@

# Move job to home directory
cp -r $TDIR/* $FULLPATH
cd $FULLPATH
rm -rf $TDIR

exit 0
