#!/usr/bin/env bash

##Mike Mariani, Frietze Lab, UVM, 2021

##Derived from Kim et al. 2020 nature paper, EBV and 4c-seq (I set mapq filter >=10 for alignment) 
##Also can check out general chipseq tutorial:
##https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html

##Galaxy submit command:
##qsub -cwd -pe threads 16 -j yes -o ~/path/to/logs/folder /~/path/to/scripts/folder/4c_window_method_pipeline.bash

################################ Set initial variables and perform alignment ###############################

##Initial vairables:
input_dir="~/path/to/trimmed/fastq/folder"
output_dir="~/path/to/desired/outupt/folder"
bowtie2_ref="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/GRCh38"

##Perform alignment with bowtie2:

cd $input_dir

for i in *.trimmed.fastq
do

base_name=$(basename $i ".trimmed.fastq")

/slipstream/home/mmariani/programs/miniconda3/bin/bowtie2 \
-p 16 \
-x $ref \
-U $i \
| samtools view -@ 15 -bS -F 4 -q 10 - \
| samtools sort -@ 15 -o $output_dir"/"$base_name".sorted.mapped.bam" -

##Index sorted, mapped bam and collect alignment stats: 
samtools index $output_dir"/"$base_name".sorted.mapped.bam"
samtools flagstat $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.flagstat"
samtools stats $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.stats"
samtools idxstats $output_dir"/"$base_name".sorted.mapped.bam" > $output_dir"/"$base_name".sorted.mapped.bam.idxstats"

done

#################################### perform peak calling ###############################

##Make sure that you have the necessary R libraries installed
##for the R script and set the path to the R libraries folder:
export R_LIBS="~/path/to/R/libraries"
##e.g.,
##export R_LIBS="/slipstream/home/mmariani/R/x86_64-pc-linux-gnu-library/3.6"

##Now use my Rscript that calculates the p-value based
##of Kim et al. windowed method:

##Seth required variables:
suffix=_"L002_R1_001.sorted.mapped.bam$" ##Can change the suffix of your files as needed
read_length=60 ##read length after trimming viewpoint
wm_output_dir="~/path/to/desired/output/folder/for/window/method" ##set output folder

##Run the newest version of script (version 2)
Rscript "~/path/to/scripts/folder/read_poisson_pvalues_mm_v2.R" \
"~/path/to/aligned/bams/from/above/folder" \
$suffix \
$read_length \
$wm_output_dir

############################ Consolidate peaks called from above ###############################

##Here we use a function that is part of macs2
##that will combine smaller peaks that are 
##within a certain distance of larger peaks

cd $wm_output_dir

for i in *.sorted.mapped.bdg
do

~/path/to/macs2 bdgpeakcall \
-i $i \
-c 5 \ ##Kim et al. default
-l 20000 \ ##Kim et al. default
-g 10000 \ ##Kim et al. default
-o $(basename $i ".bdg")".narrowPeak"

##Understand what the above parameters are doing.
##For viral integration detection in HHV-6A I used
##larger peak size parameters than I use for other
##herpesviruses such as EBV, VZV, HSV-1 - in which
##case I would start with the default parameters from 
##the Kim et al. paper. 

done
