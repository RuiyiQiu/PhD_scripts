WD='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/02_SNP'
REF="/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/00_ref/Rciliata.h1.fasta"
REF_DIR='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/00_ref/'
BAM_DIR='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/01_mapping'

cd $REF_DIR
samtools faidx $REF

cd $WD
cat callSNP.list | while read name; do
    bam=$name.bam
    bcftools mpileup -a AD,DP,SP -Ob -f $REF $BAM_DIR/$bam > $name.bcf
    bcftools call -f GQ,GP -mO z -o ./$name.vcf.gz $name.bcf
    rm $name.bcf
    bcftools index $name.vcf.gz
done
# merge vcf
bcftools merge -o merge.vcf.gz -O z *.vcf.gz
