WD="/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/01_mapping"
REF="/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/00_ref/Rciliata.h1.fasta"
READS_DIR="/biodata/dep_mercier/grp_marques/rqiu/Rciliata_Costa_Rica_Andres_WGS/nervosa"
cd $WD

bwa index $REF

cat sample.list | while read name; do
    r1="$READS_DIR/${name}/${name}_1.fq.gz"
    r2="$READS_DIR/${name}/${name}_2.fq.gz"
    
    bwa mem -t 20 $REF $r1 $r2 | \
    samtools view -Sb - | \
    samtools sort -@ 20 -o $WD/${name}.bam
    echo "Processed: $name (sorted BAM)"
done
bash $WD/stats/mapping_stats.sh
