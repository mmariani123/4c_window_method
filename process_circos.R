#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(Rsamtools)
library(parallel)
library(pbmcapply)
##library(karyoploteR) ##Error loading library, may need reinstall
library(rlist)
library(magick)

setDTthreads(threads=16)

set.seed(seed=1)

args = commandArgs(trailingOnly=TRUE)

##source("/slipstream/home/mmariani/scripts/4c_window_method/output_circos_mm.R")

source(args[4])

np.bed.files <- list.files(args[1],
				        ##"/slipstream/home/mmariani/scripts/4c_window_method/output",
                ##pattern=".narrowPeak.bed",
                pattern=args[2],
                full.names = TRUE)

##HSV1_E_1_N704_S13_L002_R1_001.sorted.mapped.narrowPeak.bed is empty for round 1 
##but not round 2
##np.bed.files <- np.bed.files[-7]

np.bed.frames <- lapply(np.bed.files, FUN=function(x){if(file.info(x)$size!=0){read.table(x,header=FALSE,sep="\t",stringsAsFactors=FALSE)}})

count=1
for(i in seq_along(np.bed.files)){
  if(file.info(np.bed.files[i])$size==0){next}
  np.bed.frames[[count]]$sample <- gsub(args[3],"",basename(np.bed.files[count]))
  np.bed.frames[[count]]$sample <- gsub(args[2],"",np.bed.frames[[count]]$sample)
  count=count+1
}

np.bed.frame <- do.call(rbind, np.bed.frames)
colnames(np.bed.frame) <- c("chrom","chromStart","chromEnd","sample")

##NC_001806.2

viewpoint.1.list.hsv1 <- list("chrHSV1", 1, 152222)

cytoband.hg38.hsv1 <- args[5]

chromosome.index.hg38.hsv1=c(
  "chr1",
  "chr2",
  "chr3",
  "chr4",
  "chr5",
  "chr6",
  "chr7",
  "chr8",
  "chr9",
  "chr10",
  "chr11",
  "chr12",
  "chr13",
  "chr14",
  "chr15",
  "chr16",
  "chr17",
  "chr18",
  "chr19",
  "chr20",
  "chr21",
  "chr22",
  "chrM",
  "chrX",
  "chrY",
  "chrHSV1")

for(i in seq_along(unique(np.bed.frame$sample))){
  np.bed.subset = subset(np.bed.frame, sample==unique(np.bed.frame$sample)[i])
  output_circos_mm(bed = np.bed.subset,
                   samples = np.bed.subset$sample,
                   viewpoint = viewpoint.1.list.hsv1,
                   filename.out = paste0(args[3],"/",basename(unique(np.bed.frame$sample)[i]),".pdf"),
                   ##filename.out = paste0("/slipstream/home/mmariani/scripts/4c_window_method/output","/",unique(np.bed.frame$sample)[i],".pdf"),
                   cytoband=cytoband.hg38.hsv1,
                   chromosome.index = chromosome.index.hg38.hsv1,
                   title = unique(np.bed.frame$sample)[i],
                   reduce_to_mid_line = FALSE,
                   lwd=2,
                   colors=c("black"))
                   ##colors = c("red","blue"))
}

##output a combined plot:
plots.files <- list.files(args[1],
                         pattern=".pdf",
                         full.names=TRUE)

plots.in <- lapply(plots.files, FUN=function(x){
  plot.now = cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(x, density=600))
  return(plot.now)
})

plots.out <- cowplot::plot_grid(plotlist = plots.in,
                                align="vh",
                                axis="lrtb",
                                ncol=2)

ggsave(plot=plots.out,
       filename=paste0(args[3],"/all.circos.together.pdf"),
       device="pdf",
       height=24,
       width=12)

##Combine narrow peaks files

##args.1 <- "/slipstream/home/mmariani/projects/hsv1_4c/output/round_1_and_round_2_combined"

np.peak.files <- list.files(path = args[1],
                            pattern=".narrowPeak$",
                            full.names = TRUE)

np.peak.list <- list()
for(j in np.peak.files){
  print(j)
  num.lines=as.numeric(unlist(strsplit(system(paste0("wc -l ",j),intern=TRUE), split=" "))[[1]])
  print(num.lines)
  if(num.lines<=1){ next }else{
    frame.now <- read.table(file=j,
                            stringsAsFactors = FALSE,
                            sep="\t",
                            skip=1,
                            header=FALSE)
    np.peak.list = list.append(np.peak.list, frame.now)
  }
}

##args.2   <- ".combined.sorted.mapped.narrowPeak.bed$"
##check.dt <- setDT(do.call(rbind, np.peak.list))[,V4:=gsub(paste0(gsub(".bed\\$","",args.2),"_narrowPeak\\d+"),"",V4,perl=TRUE),]

fwrite(file=paste0(args[3],"/hsv1.4c.combined.narrow.peaks.narrowPeak"),
       x=setDT(do.call(rbind, np.peak.list))[,V4:=gsub(paste0(gsub(".bed\\$","",args[2]),"_narrowPeak\\d+"),"",V4,perl=TRUE),],
       sep="\t",
       col.names = FALSE,
       row.names = FALSE,
       quote=FALSE)
