#!/bin/bash

#REF="/common/genomics-core/reference/BWA/hg38_WES"
REF="/common/genomics-core/data/Temp/Di_reference/GRCm38_WGS/"
PICARD="/common/genomics-core/apps/"

module load R
module load samtools

#echo "##############Faidx is started
#/hpc/apps/samtools/1.6/bin/samtools faidx $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna"
#/common/genomics-core/anaconda2/bin/samtools faidx $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna
#echo "Subject: Faidx is done for REF" | sendmail -v yizhou.wang@cshs.org

#echo "##############Create SequenceDictionary is started"
#/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar CreateSequenceDictionary R=$REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna O=$REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.dict
#echo "Subject: Create SequenceDictionary is done for REF" | sendmail -v yizhou.wang@cshs.org

echo "##############Mapping is started"
/common/genomics-core/anaconda2/bin/bwa mem -M -t 10 -R '@RG\tID:$1\tLB:$1\tPL:ILLUMINA\tPM:HISEQ\tSM:$1' $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna $2 $3 > $1_aligned_reads.sam


echo "##############Sorting BAM is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar SortSam INPUT=$1_aligned_reads.sam OUTPUT=$1_sorted_reads.bam SORT_ORDER=coordinate TMP_DIR=`pwd`/tmp


echo "##############Collect Alignment & Insert Size Metrics is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics R=$REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna I=$1_sorted_reads.bam O=$1_alignment_metrics.txt


echo "##############Collect Alignment & Insert Size Metrics is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar CollectInsertSizeMetrics INPUT=$1_sorted_reads.bam OUTPUT=$1_insert_metrics.txt HISTOGRAM_FILE=$1_insert_size_histogram.pdf


echo "##############Depth calculation is started"
/hpc/apps/samtools/1.6/bin/samtools depth -a $1_sorted_reads.bam > $1_depth_out.txt


echo "##############Sorting BAM is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar MarkDuplicates INPUT=$1_sorted_reads.bam OUTPUT=$1_dedup_reads.bam METRICS_FILE=$1_metrics.txt


echo "##############BuildBamINdex is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/picard.jar BuildBamIndex INPUT=$1_dedup_reads.bam


echo "##############RealignerTargetCreator is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_dedup_reads.bam -o $1_realignment_targets.list

echo "##############Realign Indels is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T IndelRealigner -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_dedup_reads.bam -targetIntervals $1_realignment_targets.list -o $1_realigned_reads.bam


echo "##############Call Variants is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_realigned_reads.bam -o $1_raw_variants.vcf


echo "##############Extract SNPs is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T SelectVariants -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_variants.vcf -selectType SNP -o $1_raw_snps.vcf


echo "##############Extract Index is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T SelectVariants -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_variants.vcf -selectType INDEL -o $1_raw_indels.vcf


echo "##############Filter SNP is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantFiltration -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" --filterName "basic_snp_filter" -o $1_filtered_snps.vcf


echo "##############Filter indel is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantFiltration -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "basic_indel_filter" -o $1_filtered_indels.vcf


echo "##############Base Quality Scroe Recalibration is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_realigned_reads.bam -knownSites $1_filtered_snps.vcf -knownSites $1_filtered_indels.vcf -o $1_recal_data.table


echo "##############Post Base Quality Score Recalibration is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_realigned_reads.bam -knownSites $1_filtered_snps.vcf -knownSites $1_filtered_indels.vcf -BQSR $1_recal_data.table -o $1_post_recal_data.table


echo "##############Analyze Covariates is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -before $1_recal_data.table -after $1_post_recal_data.table -plots $1_recalibration_plots.pdf


echo "##############Apply BQSR is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T PrintReads -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_realigned_reads.bam -BQSR $1_recal_data.table -o $1_recal_reads.bam


echo "##############Call Variants is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -I $1_recal_reads.bam -o $1_raw_variants_recal.vcf


echo "##############Extract SNPs & Indels is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T SelectVariants -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_variants_recal.vcf -selectType SNP -o $1_raw_snps_recal.vcf


echo "##############Extract SNPs & Indels is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T SelectVariants -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_variants_recal.vcf -selectType INDEL -o $1_raw_indels_recal.vcf


echo "##############Filter SNP is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantFiltration -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $1_filtered_snps_final.vcf

echo "##############Filter Indels is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantFiltration -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -V $1_raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $1_filtered_indels_final.vcf


echo "##############SNP anno is started"
/common/genomics-core/anaconda2/bin/java -jar ~/genomics/apps/snpEff_latest_core/snpEff/snpEff.jar -csvStats $1_snps.snpEff.summary.csv -v -s $1_snps.snpEff.summary.html GRCh38.86 $1_filtered_snps_final.vcf > $1_filtered_snps_final.ann.vcf

echo "##############indel anno is started"
/common/genomics-core/anaconda2/bin/java -jar ~/genomics/apps/snpEff_latest_core/snpEff/snpEff.jar -csvStats $1_indels.snpEff.summary.csv -v -s $1_indels.snpEff.summary.html GRCh38.86 $1_filtered_indels_final.vcf > $1_filtered_indels_final.ann.vcf

echo "##############GATK eval for SNP  is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantEval --eval $1_filtered_snps_final.vcf -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -o $1.snps.eval


echo "##############GATK eval for Indels  is started"
/common/genomics-core/anaconda2/bin/java -jar $PICARD/GenomeAnalysisTK.jar -T VariantEval --eval $1_filtered_indels_final.vcf -R $REF/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna -o $1.indels.eval
echo "Subject: Analysis is done for $1" | sendmail -v di.wu@cshs.org
#echo "Subject: Analysis is done for $1" | sendmail -v yizhou.wang@cshs.org

#bedtools genomecov -bga -ibam recal_reads.bam > genomecov.bedgraph
