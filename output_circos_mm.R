#!/usr/bin/env Rscript

##Mike Mariani UVM 2019

##My function for outputting circos plots:
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(magrittr)

output_circos_mm <- function(bed,
                             samples,
                             viewpoint,
                             filename.out,
                             cytoband,
                             chromosome.index,
                             title,
                             colors=NULL,
                             height,
                             width,
                             reduce_to_mid_line=FALSE,
                             lwd)
{
  circos.clear()
  ##bed <- sbw_13_peaks_bed
  ##viewpoint=vzv.vp.1
  ##chromosome.index <- chromosome_index
  ##cytoband <- cytoband.hg38.vzv
  ##viewpoint=viewpoint.1.list
  ##bed <- sub.bed.frame[,c(1,2,3)]
  ##chromosome.index <- chromosome.index.hg38.hhv6a
  ##cytoband <- cytoband.hg38.hhv6a
  ##colors=TRUE
  ##title="hello"
  ##samples=sub.bed.frame[,4]
  ##filename.out=filename.out
  ##png(file=filename.out,
  ##    height=480,
  ##    width=480)
  pdf(file=filename.out,
      height=height,
      width=width)
  ##HHV6(newest)
  circos.initializeWithIdeogram(
    cytoband=cytoband,
    chromosome.index = chromosome.index
  )
  if(!is.null(colors) & is.null(samples))
  {
    values <- log10(bed[,"score"])
    print("samples are null")
    if(length(values)==1)
    {
      values_range <- seq(from=min(values),to=max(values+2),length.out = 3)
    }else{
      values_range <- seq(from=min(values),to=max(values),length.out = 3)
    }
  ##values_range <- seq(from=min(values),to=max(values),length.out = 3)
    col_fun = colorRamp2(values_range, c("blue", "purple", "red"))
  }else if(!is.null(colors) & !is.null(samples)){
    ##print(samples[1])
    ##print(colors[1])
    ##print(length(colors))
    ##print(bed[bed$sample==samples[1],])
    if(length(unique(samples))!=length(colors)){stop("Number of colors must be equal to number of samples")}
    bed$color <- "black"
    for(color.count in 1:length(colors)){
      ##print(unique(bed$sample))
      ##print(color.count)
      ##print(bed[bed$sample==samples[color.count],])
      ##print(colors[color.count])
      ##print(samples[color.count])
      bed[bed$sample %in% unique(samples)[color.count],"color"] <- colors[color.count]
      bed[bed$sample %in% unique(samples)[color.count],"color"] %>% head(n=10) %>% print()
    }
    ##bed$color %>% head(n=10) %>% print()
    values <- unique(samples)
    ##col_fun = colorRamp2(c(1:length(values)),colors=sample(colors(),length(values)))
  }
  bed.2.mat <- matrix(rep(unlist(viewpoint),each=nrow(bed)),ncol=3,nrow=nrow(bed))
  bed.2 <- as.data.frame(bed.2.mat, stringsAsFactors=FALSE)
  ##print(bed.2)
  colnames(bed.2) <- c("chrom","chromStart","chromEnd")
  bed.2$chromStart <- as.integer(x=bed.2$chromStart)
  bed.2$chromEnd <- as.integer(x=bed.2$chromEnd)
  bed.2$color <- bed$color
  if(!is.null(colors) & is.null(samples)){
    bed.2$score <- bed[,"score"]
    bed.2$color <- bed[,"color"]
    circos.genomicLink(
      bed,
      bed.2,
      col = col_fun(values))
  }else if(!is.null(colors) & !is.null(samples)){
      ##print(bed$color)
      ##colors <- I(brewer.pal(length(unique(values)), name = 'Accent'))
      ##print(nrow(bed))
      ##print(nrow(bed.2))
      ##print(bed$color)
      circos.genomicLink(
      bed,
      bed.2,
      col=bed$color,
      reduce_to_mid_line=reduce_to_mid_line,
      lwd=lwd)
      ##col = factor(samples, levels=unique(samples)))
  }else{
      print(nrow(bed))
      print(nrow(bed.2))
      circos.genomicLink(
        bed,
        bed.2)
  }
  ##title(paste0(virus," interactions FDR < ",as.character(fdr)))
  title(title)
  if(!is.null(colors) & is.null(samples)){
    circle_size = unit(1, "snpc")
    lgd_links = ComplexHeatmap::Legend(at=c(round(min(values),digits=3),round((min(values)+(max(values)-min(values)/2)),digits=3),round(max(values),digits=3)), col_fun = col_fun, title_position = "topleft", title = "score")
    lgd_list_vertical = packLegend(lgd_links)
    ##pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical), 
    ##                      height = grobHeight(lgd_list_vertical), just = c("left", "center")))
    draw(lgd_list_vertical, x = unit(0.95, "npc"), y = unit(0.95, "npc"), just = c("right", "top"))
    ##upViewport()
  }else if(!is.null(colors) & !is.null(samples)){
    circle_size = unit(1, "snpc")
    lgd_links = ComplexHeatmap::Legend(labels = values, 
                                       title = "sample", 
                                       ##legend_gp = gpar(fill = 1:length(values)),
                                       legend_gp = gpar(fill = colors)
                                       ##col=colorRamp2(breaks=c(samples), colors=colors, transparency = 0, space = "LAB")
                                       )
    lgd_list_vertical = packLegend(lgd_links)
    ##pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical), 
    ##                      height = grobHeight(lgd_list_vertical), just = c("left", "center")))
    draw(lgd_list_vertical, x = unit(0.92, "npc"), y = unit(0.92, "npc"), just = c("right", "top"))
    ##upViewport()
  }
  dev.off()
  circos.clear()
}

##viewpoint.chosen <- list("chrVZV",79298,79894)
##cytoband.hg38.vzv <- "/slipstream/home/mmariani/references/4c_references/cytoband.ucsc.hg38.and.vzv.txt"
##chromosome.index.hg38.vzv=c(
##  "chr1",
##  "chr2",
##  "chr3",
##  "chr4",
##  "chr5",
##  "chr6",
##  "chr7",
##  "chr8",
##  "chr9",
##  "chr10",
##  "chr11",
##  "chr12",
##  "chr13",
##  "chr14",
##  "chr15",
##  "chr16",
##  "chr17",
##  "chr18",
##  "chr19",
##  "chr20",
##  "chr21",
##  "chr22",
##  "chrX",
##  "chrY",
##  "chrVZV")
