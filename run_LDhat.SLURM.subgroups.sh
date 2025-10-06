#!/bin/bash
#SBATCH --job-name=LDhat_sub
#SBATCH -p normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

POP=$1
CHR=$2

# Run LDhat interval
BASE=$3
LDhat='/netscratch/dep_mercier/grp_marques/bin/LDhat'

cd $BASE
mkdir ${POP}_${CHR}
cd $BASE/${POP}_${CHR}
$LDhat/interval -seq $BASE/${POP}_${CHR}.ldhat.sites -loc $BASE/${POP}_${CHR}.ldhat.locs -lk $BASE/${POP}_new_lk.txt -its 10000000 -samp 5000 -bpen 5 -theta 0.01

# Run LDhat stat
$LDhat/stat -input rates.txt -burn 20 > ${POP}_${CHR}.rec_summary.txt
