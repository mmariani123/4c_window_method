#!/usr/bin/env bash

##Mike Mariani, Frietze Lab, UVM, 2021

##Derived from Kim et al. 2020 nature paper, EBV and 4c-seq (I set mapq filter >=10 for alignment) 
##Also can check out general chipseq tutorial:
##https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html

##Galaxy submit command:
##qsub -cwd -pe threads 16 -j yes -o /slipstream/home/mmariani/projects/hsv1_4c/logs /slipstream/home/mmariani/scripts/hsv1_4c_scripts/four_c_window_method_pipeline.bash

################################ Set initial variables and perform alignment ###############################

##Initial variables:
##input_dir="~/path/to/trimmed/fastq/folder"
##output_dir="~/path/to/desired/outupt/folder"
##bowtie2_ref="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/GRCh38"

##Round 1 HSV-1 sequencing
##input_dir="/slipstream/home/mmariani/scripts/4c_window_method/data/trimmed_fastqs"
##output_dir="/slipstream/home/mmariani/scripts/4c_window_method/output"
##bowtie2_ref="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/GRCh38"

##Combined, round1 + round2:
input_dir="/slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs"
output_dir="/slipstream/home/mmariani/projects/hsv1_4c/output/round_1_and_round_2_combined"
bowtie2_ref="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/GRCh38"

#####################################################################################################

##Perform read trimming:
##cd $input_dir
##if files are in individual folders:
##cp $(find ./  name=".fastq.gz" | grep ".gz") ./
##you can then delete the subfolders as long as you remember which file is which.
##for i in *.fastq 
##do 
##awk '{if(NR%4==2 || NR%4==0){print substr($0,21);}else{print;}}' $i > $(basename $i ".fastq")".trimmed.fastq"
##done

#####################################################################################################

##Combined multiple (2 here) rounds of trimmed files if necessary:
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_A_1_N701_S10_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_A_1_N701_S1_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_A_1.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_A_2_N707_S16_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_A_2_N707_S7_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_A_2.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_B_1_N702_S11_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_B_1_N702_S2_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_B_1.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_B_2_N708_S17_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_B_2_N708_S8_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_B_2.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_D_1_N703_S12_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_D_1_N703_S3_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_D_1.comibned.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_D_2_N709_S18_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_D_2_N709_S9_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_D_2.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_E_1_N704_S13_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_E_1_N704_S4_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_E_1.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_E_2_N710_S19_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_E_2_N710_S10_L001_R1_001.trimmed.fastq > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_E_2.combined.fastq
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_F_1_N705_S14_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_F_1_N705_S5_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_F_1.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_F_2_N711_S20_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_F_2_N711_S11_L001_R1_001.trimmed.fastq > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_F_2.combined.fastq
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_G_1_N706_S15_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_G_1_N706_S6_L001_R1_001.trimmed.fastq  > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_G_1.combined.fastq 
##cat /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_1/trimmed_fastqs/HSV1_G_2_N712_S21_L002_R1_001.trimmed.fastq /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/round_2/Frietze_SK_15538_021121/trimmed_fastqs/HSV1_G_2_N712_S12_L001_R1_001.trimmed.fastq > /slipstream/home/mmariani/projects/4c_data_and_logs/data/4c_seq_illumina/hsv1_sK/combined/trimmed_fastqs/HSV1_G_2.combined.fastq

#####################################################################################################

##Perform alignment with bowtie2:

cd $input_dir

##for i in *.trimmed.fastq
##do
##
##base_name=$(basename $i ".trimmed.fastq")
##
##/slipstream/home/mmariani/programs/miniconda3/bin/bowtie2 \
##-p 16 \
##--met-file $output_dir"/"$base_name".bt2.metrics.txt" \
##-x $bowtie2_ref \
##-U $i \
##| samtools view -@ 15 -bS -F 4 -q 10 - \
##| samtools sort -@ 15 -o $output_dir"/"$base_name".sorted.mapped.bam" -
##
####Index sorted, mapped bam and collect alignment stats: 
##samtools index $output_dir"/"$base_name".sorted.mapped.bam"
##samtools flagstat $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.flagstat"
##samtools stats $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.stats"
##samtools idxstats $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.idxstats"
##
##done
##
###################################### perform peak calling ###############################
##
####Make sure that you have the necessary R libraries installed
####for the R script and set the path to the R libraries folder:
####e.g.,
##export R_LIBS="/slipstream/home/mmariani/R/x86_64-pc-linux-gnu-library/3.6"
##
####Now use my Rscript that calculates the p-value based
####of Kim et al. windowed method:
##
####Required variables:
####suffix="_L002_R1_001.sorted.mapped.bam$" ##Can change the suffix of your files as needed
####read_length=60 ##read length after trimming viewpoint
####wm_output_dir="~/path/to/desired/output/folder/for/window/method" ##set output folder
##
##suffix=".combined.sorted.mapped.bam$" ##Can change the suffix of your files as needed
##read_length=65 ##read length after trimming viewpoint
##wm_output_dir=$output_dir ##set output folder
##
####Run the newest version of script (version 2)
####Rscript "~/path/to/scripts/folder/read_poisson_pvalues_mm_v2.R" \
####"~/path/to/aligned/bams/from/above/folder" \
####$suffix \
####$read_length \
####$wm_output_dir
##
##Rscript "/slipstream/home/mmariani/scripts/hsv1_4c_scripts/read_poisson_pvalues_mm_v2.R" \
##$wm_output_dir \
##$suffix \
##$read_length \
##$wm_output_dir
##
################################ Consolidate peaks called from above ###############################
##
####Here we use a function that is part of macs2
####that will combine smaller peaks that are 
####within a certain distance of larger peaks
##
##cd $wm_output_dir
##
##for i in *.sorted.mapped.bdg
##do
##
####Use the default parameters from Kim et al. 
####because we are looking at another episomal 
####herpesvirus (like EBV): HSV-1 in this case.
##/slipstream/home/mmariani/programs/miniconda3/bin/macs2 bdgpeakcall \
##-i $i \
##-c 5 \
##-l 20000 \
##-g 10000 \
##-o $(basename $i ".bdg")".narrowPeak"
##
####Understand what the above parameters are doing.
####For viral integration detection in HHV-6A I used
####larger peak size parameters than I use for other
####herpesviruses such as EBV, VZV, HSV-1 - in which
####case I would start with the default parameters from 
####the Kim et al. paper. 
##
##done
##
##for i in *.narrowPeak 
##do 
##
##grep -v "track" $i | awk '{print  $1"\t"$2"\t"$3;}' > $(basename $i)".bed" 
##
##done
##
####################### output circos plots ####################################

circos_output_dir=$output_dir
circos_input_dir=$output_dir

Rscript "/slipstream/home/mmariani/scripts/hsv1_4c_scripts/process_circos.R" \
$circos_input_dir \
".combined.sorted.mapped.narrowPeak.bed$" \
$circos_output_dir \
"/slipstream/home/mmariani/scripts/hsv1_4c_scripts/output_circos_mm.R" \
"/slipstream/home/mmariani/scripts/4c_window_method/cytoband.ucsc.hg38.and.hsv1.txt"

####################### Run final multiqc on output files #######################

/slipstream/home/mmariani/programs/miniconda3/bin/multiqc -f -o $output_dir
