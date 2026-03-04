#!/bin/bash
#SBATCH -p cpu 
#SBATCH -n 1
#SBATCH --mem=50G 
#SBATCH -t 01:00:00 
#SBATCH --array=0-799%800 
#SBATCH --job-name=whamMSE_5yrSENS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

cd work/pi_gfay_umassd_edu/Wulfing/CEFI_Draft2/5Yr/SensitivityTest

#set seeds 
seeds=(124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224)

#match seeds with block + sim number 
sim_id=$SLURM_ARRAY_TASK_ID
block=$((sim_id / 100 + 1))
index_in_block=$((sim_id % 100))
seed=${seeds[$index_in_block]}

#index running message files/logs by block_seed_sim
mkdir -p logs
logfile="logs/block${block}_seed${seed}_sim${sim_id}.out"
exec > "$logfile" 2>&1
scontrol update JobId=$SLURM_JOB_ID JobName=blk${block}_seed${seed}_sim${sim_id}

#output file directory 
outdir="Outputs/block${block}"
mkdir -p "$outdir"

#load conda
module load conda/latest
conda activate CONDATEST 

echo "Simulation $sim_id using seed $seed in block $block"

#run R script
Rscript 5Yr_Cluster.R \
  --seed $seed \
  --sim_id $SLURM_ARRAY_TASK_ID \
  --block $block \
  --outdir $outdir