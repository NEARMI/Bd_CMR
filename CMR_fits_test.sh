#!/bin/bash

export IFS=","

cat which_model_fit.csv | while read a b c d; do 

job_file="CMR_${a}${b}_pop.job"

echo "#!/bin/bash
#SBATCH  -p normal
#SBATCH --account=nearmi
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=3      
#SBATCH --time=${c}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkain@usgs.gov
#SBATCH --output=CMR_${a}${b}pop.log

ml R/4.1.1
Rscript top_level_script.R "$a" "$b" 1" > $job_file

    sbatch $job_file

done








