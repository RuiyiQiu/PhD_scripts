#!/bin/bash
#SBATCH --job-name=LDhat_CostaRica
#SBATCH -p normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

CHR=$1

# Run LDhat interval
BASE=$2
LDhat='/netscratch/dep_mercier/grp_marques/bin/LDhat'

cd $BASE
mkdir ${CHR}
cd $BASE/${CHR}
$LDhat/interval -seq $BASE/${CHR}.ldhat.sites -loc $BASE/${CHR}.ldhat.locs -lk $BASE/CostaRica_new_lk.txt -its 10000000 -samp 5000 -bpen 5 -theta 0.01

# Run LDhat stat
$LDhat/stat -input rates.txt -burn 20 > ${CHR}.rec_summary.txt
