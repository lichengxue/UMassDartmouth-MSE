#!/bin/bash
#SBATCH -p cpu 
#SBATCH -n 1
#SBATCH --mem=200G 
#SBATCH -t 01:00:00 
#SBATCH --job-name=Proj_analysis
#SBATCH --output=outputs.out
#SBATCH --error=errors.out

cd /work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/Projections

#output file directory 
outdir="Outputs"
mkdir -p "$outdir"

#load conda
module load conda/latest
conda activate CONDATEST 

#run R script
Rscript ProjectionsAnalysis_Cluster.R \
  --outdir $outdir