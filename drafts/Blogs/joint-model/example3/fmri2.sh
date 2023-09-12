#!/bin/bash
#PBS -l select=1:ncpus=16:mem=10gb
#PBS -l walltime=20:00:00  
#PBS -k oe

source /etc/profile.d/modules.sh  
module load R/3.6.3  

cd /home/rji870/blogs
  
Rscript --no-save --no-restore fmri2.R

exit 0