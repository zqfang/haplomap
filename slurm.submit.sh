#!/bin/bash
#SBATCH -J snakeflow
#SBATCH --time=120:00:00
#SBATCH --qos long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p normal
#SBATCH --mem=1g
####SBATCH --mail-type=FAIL
####SBATCH --mail-user=fztsing@126.com


# activate conda enviroment
source activate base

# run jobs
snakeflow="/scratch/users/fangzq/HBCGM/workflows"

ls /scratch/users/fangzq/mpd_clean.split.*.txt | while read tags 
do
    echo -e "split: ${tags}"
    snakemake --nolock -k -p -j 666 -s ${snakeflow}/haplomap.smk \
        --configfile ${snakeflow}/config.yaml \
        --cluster "sbatch --time={cluster.time_min} -p {cluster.partition} --mem={cluster.mem_mb} -c {cluster.cpus} " \
        --cluster-config ${snakeflow}/slurm_config.yaml \
        --config TRAIT_IDS=${tags}
        
    echo -e "finished ${tags}"
done