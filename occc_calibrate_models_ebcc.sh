#!/bin/sh

#$ -N log_weights
#$ -S /bin/bash
#$ -cwd
#$ -j y

#$ -t 1-395:1

echo $SGE_TASK_ID

time R --vanilla --quiet --args $SGE_TASK_ID < occc_run_models_ebcc.r > log_weights$SGE_TASK_ID.log

