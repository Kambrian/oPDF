#!/bin/bash

export ESTIMATOR=$1
dir=exe/mockfit_mc$ESTIMATOR
mkdir -p $dir
make lib -B
cp dynio.py libdyn.so mockML.py job.bsub $dir
cd $dir
bsub <job.bsub