# nohup bash mask_TRC_1.sh > mask_TRC_1.log 2>&1 &

WD='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat'

REPEAT='/netscratch/dep_mercier/grp_marques/rqiu/Rciliata_6041A/RepeatExplorer/Rciliata.h1.chr.TRC_1.sort.bed'
# generation: 
# grep '_h1' Rciliata.TRC_1.bed > Rciliata.h1.chr.TRC_1.bed
# sort -k1,1 -k2,2n Rciliata.h1.chr.TRC_1.bed > Rciliata.h1.chr.TRC_1.sort.bed

CHR_SIZE='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/chr_size.txt'
IN_VCF='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/CostaRica.thin1000.vcf.gz'
VCF=$WD/CostaRica.flt.noTRC1.vcf.gz
POPMAP='/netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/04_mask_repeat/00_test/CostaRica.subgroup.txt'

BIN=500000
BIN_N='500kb'

cd $WD

# Remove SNPs in TRC_1 from vcf
bcftools view -T ^$REPEAT -Oz -o CostaRica.flt.noTRC1.vcf.gz $IN_VCF &
bcftools index CostaRica.flt.noTRC1.vcf.gz

# generate nonrepeat.bed
bedtools complement -i $REPEAT -g $CHR_SIZE > nonrepeats.noheader.bed
## add header
sed -i 'CHR\tSTART\tEND' nonrepeats.noheader.bed
awk 'BEGIN { print "CHR\tSTART\tEND" } { print }' nonrepeats.noheader.bed > nonrepeats.bed

# file configuration
mkdir -p 01_LDhat
mkdir -p 02_pi
mkdir -p 03_Fst

# subgroup
bcftools query -l $VCF > samples.list # print vcf sample names
while read -r sample; do
    id=$(basename $sample | cut -d. -f1)   # extract sample name from bam file directory
    pop=$(grep -w $id $POPMAP | cut -f2)
    echo -e "${sample}\t${pop}"
done < samples.list > $WD/pop.info
awk '$2=="South"{print $0}' pop.info > South.list
awk '$2=="North"{print $0}' pop.info > North.list
#conda activate py-popgen
vcftools --gzvcf $VCF \
    --keep South.list \
    --recode --stdout | bgzip -c > South.vcf.gz &
vcftools --gzvcf $VCF \
    --keep North.list \
    --recode --stdout | bgzip -c > North.vcf.gz &


# LDhat
cd $WD/01_LDhat
mkdir -p $WD/01_LDhat/all
mkdir -p $WD/01_LDhat/subgroups
## for all CostaRica
cp /netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/04_LDhat/CostaRica_new_lk.txt $WD/01_LDhat/all/CostaRica_new_lk.txt
for CHR in $(cat $WD/chr.list)
do
  # split vcf file by chr, and convert to LDhat input format
  vcftools --gzvcf $VCF --chr $CHR --ldhat-geno --out all/$CHR
  # submit SLURM job for each chr
  sbatch $WD/run_LDhat.SLURM.all.sh $CHR $WD/01_LDhat/all
done
## for subgroups
cp /netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/04_LDhat/subpop/South_new_lk.txt $WD/01_LDhat/subgroups/South_new_lk.txt
cp /netscratch/dep_mercier/grp_marques/rqiu/Rnervosa_Costa_Rica_WGS_analysis/mapping_c_h1_final/03_pop_analyses/04_LDhat/subpop/North_new_lk.txt $WD/01_LDhat/subgroups/North_new_lk.txt
for POP in North South; do
    for CHR in $(cat $WD/chr.list); do
        vcftools --gzvcf $VCF \
                 --keep $WD/${POP}.list \
                 --chr $CHR \
                 --ldhat-geno \
                 --out subgroups/${POP}_${CHR}
        sbatch $WD/run_LDhat.SLURM.subgroups.sh $POP $CHR $WD/01_LDhat/subgroups
    done
done

## clean: only keep res.txt
rm *.out
rm */*/type_table.txt
rm */*/new_lk.txt
rm */*/bounds.txt
rm */*/rates.txt
# process res.txt into bed format
# Rscript process_LDhat.R <BIN> <.locs> <LDhat_res.txt> <chr> <OUT>
for CHR in $(cat $WD/chr.list)
do
  Rscript $WD/process_LDhat.R $BIN $WD/01_LDhat/all/${CHR}.ldhat.locs $WD/01_LDhat/all/${CHR}/res.txt ${CHR} CostaRica_${CHR}.$BIN_N.windowed.bed
done
for POP in North South; do
    for CHR in $(cat $WD/chr.list); do
        Rscript $WD/process_LDhat.R $BIN $WD/01_LDhat/subgroups/${POP}_${CHR}.ldhat.locs $WD/01_LDhat/subgroups/${POP}_${CHR}/res.txt ${CHR} ${POP}_${CHR}.$BIN_N.windowed.bed
    done
done
## merge results in bed format
cat CostaRica* >CostaRica.$BIN_N.windowed.ldhat.bed
cat North* >North.$BIN_N.windowed.ldhat.bed
cat South* >South.$BIN_N.windowed.ldhat.bed
# Remove windows overlapping repeats
bedtools intersect -v -a CostaRica.$BIN_N.windowed.ldhat.bed -b $REPEAT > CostaRica.$BIN_N.masked.ldhat.bed
bedtools intersect -v -a North.$BIN_N.windowed.ldhat.bed -b $REPEAT > North.$BIN_N.masked.ldhat.bed
bedtools intersect -v -a South.$BIN_N.windowed.ldhat.bed -b $REPEAT > South.$BIN_N.masked.ldhat.bed


# pi
cd $WD/02_pi
vcftools --gzvcf $VCF \
  --bed $WD/nonrepeats.bed \
  --window-pi $BIN \
  --window-pi-step $BIN \
  --out CostaRica.$BIN_N
# for subgroups
vcftools --gzvcf $WD/South.vcf.gz \
  --bed $WD/nonrepeats.bed \
  --window-pi $BIN \
  --window-pi-step $BIN \
  --out South.$BIN_N
vcftools --gzvcf $WD/North.vcf.gz \
  --bed $WD/nonrepeats.bed \
  --window-pi $BIN \
  --window-pi-step $BIN \
  --out North.$BIN_N
for NAME in CostaRica South North
do
  awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' $NAME.$BIN_N.windowed.pi > $NAME.$BIN_N.windowed.pi.bed
  bedtools intersect -v -a $NAME.$BIN_N.windowed.pi.bed -b $REPEAT > $NAME.$BIN_N.masked.pi.bed
done


# Fst
cd $WD/03_Fst
vcftools --gzvcf $VCF \
  --bed $WD/nonrepeats.bed \
  --weir-fst-pop $WD/South.list \
  --weir-fst-pop $WD/North.list \
  --fst-window-size $BIN \
  --fst-window-step $BIN \
  --out CostaRica.$BIN_N
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' CostaRica.$BIN_N.windowed.weir.fst > CostaRica.$BIN_N.windowed.fst.bed  # weighted FST
bedtools intersect -v -a CostaRica.$BIN_N.windowed.fst.bed -b $REPEAT > CostaRica.$BIN_N.masked.fst.bed
# SNP number
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4}' CostaRica.$BIN_N.windowed.weir.fst > CostaRica.$BIN_N.windowed.SNP.bed
bedtools intersect -v -a CostaRica.$BIN_N.windowed.SNP.bed -b $REPEAT > CostaRica.$BIN_N.masked.SNP.bed


# masked regions marked as NA in .bedGraph 
cd $WD
mkdir 04_summary
bedtools makewindows -g $CHR_SIZE -w $BIN > Rciliata.h1.chr.$BIN_N.windows
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 01_LDhat/CostaRica.$BIN_N.masked.ldhat.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/rho_all.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 01_LDhat/South.$BIN_N.masked.ldhat.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/rho_s.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 01_LDhat/North.$BIN_N.masked.ldhat.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/rho_n.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 02_pi/CostaRica.$BIN_N.masked.pi.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/pi_all.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 02_pi/South.$BIN_N.masked.pi.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/pi_s.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 02_pi/North.$BIN_N.masked.pi.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/pi_n.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 03_Fst/CostaRica.$BIN_N.masked.fst.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/Fst.$BIN_N.intersect.bedGraph
bedtools intersect -a Rciliata.h1.chr.$BIN_N.windows -b 03_Fst/CostaRica.$BIN_N.masked.SNP.bed -loj \
  | awk '{ val=$7; if(val=="."||val=="") val="NA"; print $1, $2, $3, val }' OFS="\t" \
  > 04_summary/SNP.$BIN_N.intersect.bedGraph

# then make .list as signals for multi tracks plotting
mkdir -p $WD/04_summary/signals
cd $WD/04_summary
ls *.bedGraph | while read FILE
do
  cut -f4 $FILE > $WD/04_summary/signals/$(basename $FILE ".$BIN_N.intersect.bedGraph").list # use double quotes to expand variables!
done



