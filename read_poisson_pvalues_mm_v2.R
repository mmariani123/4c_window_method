#!/usr/bin/env Rscript

##Mike Mariani Frietze Lab UVM 2021

##Updated version 2, 01/22/2021
##This should do a better job with iterating, 
##and also handles the ends of the chromosomes which
##will be binned at under 10kbp

##Calculate poisson probabilites from hg38-aligned
##4C data (10kb windows)

##Good discussion:
##https://math.stackexchange.com/questions/1315565/finding-a-p-value-based-on-the-poisson-distribution#:~:text=P(X%3Dx)%3D,e%E2%88%92%CE%BB%CE%BBxx!
##Also, https://rstudio-pubs-static.s3.amazonaws.com/456645_107fa2aa82de4b1da6c78c418bab9fe9.html#:~:text=R%20function%20dpois(x%2C%20lambda,TRUE%20for%20left%20tail%2C%20lower.

##Note that while a bam file is 0-based and a sam file is 1-based
##Rsamtools is 1-based

##Note the below only outputs calculations for a 
##window if at least one read is in that region, 
##thus it will not output a 10kb regions of all zeros 
## to the output. 

library(data.table)
library(dplyr)
library(Rsamtools)
library(parallel)
library(pbmcapply)
##library(karyoploteR) ##Error loading karyoploteR, may need re-install

setDTthreads(threads=16)
set.seed(seed=1)
args = commandArgs(trailingOnly=TRUE)
##args[1] = input.dir
##args[2] = filename pattern
##args[3] = read length
##args[4] = output.dir

calc_prob_mm <- function(counts.now,mean){
  prob.out=1-sum(dpois(0:counts.now,lambda=mean))
  ##Note soometimes the above sum is double.precision
  ##and if the sum is very close to 1 1L-sum 
  ##can lead to a negative value; therefore,
  ##if negative, just set 1-sum= smallest value. 
  ##This makes sense as probabilities should not be negative.
  if(prob.out<=0){prob.out=.Machine$double.xmin}
  return(prob.out)
}

poisson_prob_mm_v2 <- function(bam.file, read.length, out.dir){

print(paste0("processing bam file ",bam.file," ... "))
  
##bam.file = bam.files[[15]]
##read.length=60
##out.dir="/slipstream/home/mmariani/projects/vzv_interactions_summary/mm1_poisson_bedgraphs"

##http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
##ucsc hg38 chrom widths:
hg38.chrom.sizes <- data.frame(
  chrom=c("chr1",	
          "chr2",	
          "chr3",	
          "chr4",
          "chr5",
          "chr6",
          "chr7",
          "chr8",
          "chr9",
          "chr11",	
          "chr10",
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
          "chrX",
          "chrY"),
          size=c(248956422,
                 242193529,
                 198295559,
                 190214555,
                 181538259,
                 170805979,
                 159345973,
                 145138636,
                 138394717,
                 135086622,
                 133797422,
                 133275309,
                 114364328,
                 107043718,
                 101991189,
                 90338345,
                 83257441,
                 80373285,
                 58617616,
                 64444167,
                 46709983,
                 50818468,
                 156040895,
                 57227415),
  stringsAsFactors = TRUE
)

bam.check <- scanBam(bam.file)
bam.now   <- data.frame(chrom=bam.check[[1]][[3]],
                        pos=bam.check[[1]][[5]],
                        stringsAsFactors = FALSE)
bam.now$chrom <- as.character(bam.now$chrom)
chrom.pos     <- setDT(bam.now)[, .N, by=.(pos,chrom)]
sub.chrom.pos <- chrom.pos[chrom.pos$chrom %in% hg38.chrom.sizes$chrom]

##Version2, let's speed it up even more, does not do single bp position,
##cosiders each mapped read as a count of 1
poisson.outs <- list()
count.now <- 1
for(i in 1:nrow(hg38.chrom.sizes))
{
  chrom.bam = sub.chrom.pos[sub.chrom.pos$chrom %in% hg38.chrom.sizes$chrom[i]]
  ##for(j in seq(from=1,by=1e4,to=hg38.chrom.sizes[i,"size"])+1e4)
  for(j in seq(from=1,by=1e4,to=hg38.chrom.sizes[i,"size"]))
  {
    ##check if at the end of chromosome:
    remainder=hg38.chrom.sizes$size[i]%%1e4
    ##print(remainder)
    if(j+ceiling(remainder)-1==hg38.chrom.sizes$size[i]){
      ##stop("error")
      ##report status:
     if(hg38.chrom.sizes$size[i]==(j+ceiling(remainder)-1)){
       print(paste0("chrom size: ",
                    hg38.chrom.sizes$size[i], 
                    " and iterator + remainder: ",
                    (j+ceiling(remainder)-1),
                    " should be the same number. Is this true? ... ",
                    (hg38.chrom.sizes$size[i]==(j+ceiling(remainder)-1))
                   )
       )
     }else{stop("serious eror, window iteration is not matching chroms sizes")}
      
        empty.scores <- data.table(pos=seq(from=j,by=1,to=hg38.chrom.sizes$size[i]),
                                   chrom=as.character(hg38.chrom.sizes[i,"chrom"]),
                                   score=0,
                                   counts=rep(0,times=hg38.chrom.sizes$size[i]-j+1)
      )
      sub.bam = chrom.bam[pos >= j & pos < j+1e4]
      ##sub.bam.joined = dplyr::left_join(empty.scores,sub.bam,by=c("pos"))[,c(1,2,5)]
      ##sub.bam.joined[is.na(sub.bam.joined)] <- 0
      colnames(sub.bam) <- c("pos", "chrom", "counts")
      ##sub.bam.joined$counts <- as.numeric(sub.bam.joined$counts)
      ##The "score" is the probability 
      sub.bam[,score:=ifelse(counts==0,.Machine$double.xmin,-log10(calc_prob_mm(counts.now=counts,mean=sum(sub.bam$counts)/(hg38.chrom.sizes$size[i]-j+1)))),by="counts"]
      ##sub.bam[,p.value:=ifelse(counts==0,1,sum(dpois(0:counts,lambda=mean(sub.bam$counts)))),by=c("pos")]
      if(nrow(sub.bam)!=0)
      {
        if((any(!is.finite(sub.bam$score)) | any(is.na(sub.bam$score)) | any(sub.bam$score < 0))==TRUE){
          stop()
        }else{
          sub.bam <- dplyr::full_join(empty.scores,sub.bam[,c(2,1,4,3)],by=c("chrom","pos"))[,c(1,2,5,6)]
          colnames(sub.bam) <- c("pos","chrom","scores","counts")
          sub.bam[is.na(sub.bam)] <- 0
          poisson.outs[[count.now]] <- sub.bam
          count.now<-count.now+1
        }
      }
    }else{  ##If not yet reached end portion of chromosome:
      ##remainder=hg38.chrom.sizes$size[1]%%10e4
      ##print(paste0("remainder = ",remainder))
      ##print(hg38.chrom.sizes$size[1]-remainder)
      ##Report staus every 10 million bp:
      if((j-1)%%1e7==0){print(paste0(hg38.chrom.sizes[i,"chrom"],":",j))}
      
      empty.scores <- data.table(pos=seq(j,(j+1e4)-1),
                                 chrom=as.character(hg38.chrom.sizes[i,"chrom"]),
                                 score=0,
                                 counts=rep(0,times=1e4)
      )
      sub.bam = chrom.bam[pos >= j & pos < j+1e4]
      ##sub.bam.joined = dplyr::left_join(empty.scores,sub.bam,by=c("pos"))[,c(1,2,5)]
      ##sub.bam.joined[is.na(sub.bam.joined)] <- 0
      colnames(sub.bam) <- c("pos", "chrom", "counts")
      ##sub.bam.joined$counts <- as.numeric(sub.bam.joined$counts)
      ##The "score" is the probability 
      sub.bam[,score:=ifelse(counts==0,.Machine$double.xmin,-log10(calc_prob_mm(counts.now=counts,mean=sum(sub.bam$counts)/1e4))),by="counts"]
      ##sub.bam[,p.value:=ifelse(counts==0,1,sum(dpois(0:counts,lambda=mean(sub.bam$counts)))),by=c("pos")]
      if(nrow(sub.bam)!=0)
      {
        if((any(!is.finite(sub.bam$score)) | any(is.na(sub.bam$score)) | any(sub.bam$score < 0))==TRUE){
          stop()
        }else{
          sub.bam <- dplyr::full_join(empty.scores,sub.bam[,c(2,1,4,3)],by=c("chrom","pos"))[,c(1,2,5,6)]
          colnames(sub.bam) <- c("pos","chrom","scores","counts")
          sub.bam[is.na(sub.bam)] <- 0
          poisson.outs[[count.now]] <- sub.bam
          count.now<-count.now+1
        }
      }
    }
  }
}

poisson.final <- do.call(rbind,poisson.outs)
##convert from 1-based to 0-based below:
poisson.final$chromStart <- poisson.final$pos-1
poisson.final$chromEnd <- poisson.final$pos
poisson.final <- poisson.final[,c(2,5,6,4,3)]

##One post had fwrite running almost
##20x faster than write.table
data.table::fwrite(x=poisson.final[,c(1,2,3,5)], 
                   file = paste0(out.dir,
                                 "/",
                                 gsub(".bam","",basename(bam.file)),
                                 ".bdg"), 
                   append = FALSE, 
                   quote = "auto",
                   col.names = TRUE,
                   sep = "\t",
                   nThread = 16,
                   showProgress=FALSE)

##return(poisson.final)

}

##bam.files <- list.files(path="/slipstream/home/mmariani/projects/hhv6_detection/hhv6a_illumina_mm1_hg38_aligned_bowtie2_for_nature_method",
##                        pattern="_L002_R1_001.sorted.mapped.bam$",
##                        full.names = TRUE)
##out.dir="/slipstream/home/mmariani/projects/hhv6_detection/mm1_illumina_nature_method_peaks/poisson_v2_bedgraphs"
##lapply(bam.files,
##       poisson_prob_mm_v2,
##       read.length=60,
##       out.dir=out.dir)

##bam.check <- scanBam(bam.files[[1]])
##Look at the fields:
##names(bam.check[[1]])
##bam.check[[1]][[5]] ##left-most mapping position
##bam.check[[1]][[7]] ##mapq score

##For batching:
bam.files <- list.files(path=args[1],
                        pattern=args[2],
                        full.names = TRUE)

print(bam.files)
lapply(bam.files,
	   poisson_prob_mm_v2,
	   read.length=args[3], ##minus viewpoint+enzyme *usually minus 20bp
	   out.dir=args[4]
       )
				  
##saveRDS(outputs,paste0(args[4],"/","poisson.bdg.object.RDS"))
