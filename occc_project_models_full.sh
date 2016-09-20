#!/bin/sh

#$ -N proj_mod
#$ -S /bin/bash
#$ -cwd
#$ -j y

#$ -t 1-124:1

echo $SGE_TASK_ID

time R --vanilla --quiet --args $SGE_TASK_ID < project_afe_on_full_climate_data.r > proj_mod$SGE_TASK_ID.log

