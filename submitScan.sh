#!/bin/bash

make -B scan ESTIMATOR=$1 
exe=scanHuge$1
mv scan exe/$exe
cp job.tmpl exe/job.bsub
cd exe
echo ./$exe >>job.bsub
bsub <job.bsub