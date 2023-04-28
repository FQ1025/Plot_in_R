#library(Gviz) #require for IGrange/GenomicRanges/GenomeInforDb/grid
library(GenomicRanges)
library(dplyr)
library(patchwork)
library(rtracklayer)
library(tidyverse)

library(ggbio)
library(ggplotify)
library(biovizBase)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationFilter)

library(cowplot)
library(scales)

###----functions ------
ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
plot.track = function(region = GRanges(),
                      bigwigfile = character(),
                      downsample = T,
                      log.scale = T,
                      max.downsample = 5000,
                      downsample.rate = 0.1,
                      smooth = 200,
                      track.col = "red",
                      plot.type = "bar",
                      heat.bin = 20,
                      y.axis = c(0,100),
                      y.label = character(),
                      seed = 42
){
  chr = rtracklayer::import(con = bigwigfile, as = "NumericList", which = region )[[1]]
  chr = RcppRoll::roll_mean(x = chr, n = smooth, fill = 0L) # smooth line
  tmp = data.frame(
    position = start(region):end(region),
    score = chr,
    stringsAsFactors = FALSE)
  if( downsample ){
    window.size = width( region)
    sampling = min(max.downsample, window.size * downsample.rate)
    set.seed(seed)
    tmp = slice_sample(.data = tmp, n = sampling)
  }
  if( log.scale ){
    tmp$score = log10(tmp$score+1)
  }
  
  a = ggplot(tmp,aes(x = position, y = score))
  
  if(plot.type=="heatmap"){
    tmp$bin = floor(x = tmp$position / heat.bin )
    tmp = tmp %>% group_by(bin) %>% mutate( score = mean(score))
    tmp = unique(tmp[, c("bin", "score")])
    a = ggplot( tmp, aes(x = bin, y = 1, fill = score)) + 
      geom_tile() + scale_fill_viridis_c(option = "plasma") + # tokyo 
      theme_cowplot() + 
      ylab(y.label)+
      scale_y_continuous( expand = c(0,0) ) +
      scale_x_continuous( expand = c(0,0) ) +
      theme(axis.text = element_blank(),axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0,hjust = 1,vjust = 0.5),
            axis.ticks = element_blank(),axis.line = element_blank(),
            text = element_text(family = "sans",size = 20))
    return(a)
  }else if(plot.type=="coverage"){
    a = a + geom_area(fill=track.col,color=track.col)
  }else if(plot.type=="line"){
    a = a + geom_line()
  }else if(plot.type=="bar"){
    a = a + geom_bar(stat="identity",fill=track.col,color=track.col)
  }
  a = a + theme_cowplot() + 
    ylab(y.label)+
    scale_y_continuous( limits = y.axis,expand = c(0,0),breaks = y.axis[2],
                        labels = paste0("[",y.axis[1],"-",y.axis[2],"]") ) +
    scale_x_continuous( limits = c(start(region),end(region)),expand = c(0,0) ) +
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          axis.title.y = element_text(angle = 0,hjust = 1,vjust = 0.5),
          axis.ticks = element_blank(),axis.line.x = element_blank(),
          axis.text.y = element_text(size = 14),
          text = element_text(family = "sans",size = 20))
  a
  return(a)
}
plot.gene = function(database = NULL,
                     names.expr = "gene_name",
                     genome.range = GRanges(),
                     expand.for.symbol = c(0,2.5),
                     symbol.size = 4,
                     gene.col = "#2d266c",
                     gene.list = character(),
                     label = F){
  if(!is_empty(gene.list)){
    gtf = read.delim( "/lustre/user/liclab/fangq/ref/UCSC/mm10/Annotation/Genes/mm10.gene.bed",
                      header = F,sep = "\t",col.names = c("chr","start","end","ncbi","strand","symbol") )
    tmp = GRanges(gtf %>% dplyr::filter(chr ==paste0("chr",seqnames(genome.range)) & start>=start(genome.range) & end <= end(genome.range)) %>%
                    dplyr::filter(symbol %in% gene.list ))
    wh = range(tmp,ignore.strand = T)
    b = autoplot(database,GeneNameFilter(gene.list),
                 stat = "reduce",
                  names.expr = names.expr,label = label,
                  axis.text.x = NULL,layout = "linear",coord = "genome",#geom = "full",
                  label.color = "black",color = gene.col,fill = gene.col,#,#(collapse)
                  label.size = symbol.size,indel.color = "red") +
    scale_x_continuous( limits = c(start(genome.range),end(genome.range)),expand = c(0,0) ) +
    scale_y_continuous( expand = expand.for.symbol ) +
    theme_cowplot() +
    xlab(paste0("Chr",seqnames(genome.range)," : ",start(genome.range )," - ",end(genome.range))) +
    ylab("Gene") + 
    theme(axis.title.x = element_text(size = 22),
          axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
          axis.title.y = element_text(size = 20,angle = 0,vjust = 0.5,hjust = 1))
  b
  return(b)
  }
  b = autoplot(database, 
               which = genome.range, 
               axis.text.x = NULL,
               label.color = "black",
               color = gene.col,
               fill = gene.col,label = label,
               names.expr = names.expr,
               label.size = symbol.size,
               indel.color = "red",
               label.repel = T) +
    scale_x_continuous( limits = c(start(genome.range),end(genome.range)),expand = c(0,0) ) +
    scale_y_continuous( expand = expand.for.symbol ) +
    theme_cowplot() +
    xlab(paste0("Chr",seqnames(genome.range)," : ",start(genome.range )," - ",end(genome.range))) +
    ylab("Gene") + 
    theme(axis.title.x = element_text(size = 22),
          axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
          axis.title.y = element_text(size = 20,angle = 0,vjust = 0.5,hjust = 1))
  b
  return(b)
}

plot_isulation = function( GR = GRanges(),
                           sample.name = character()){
  require(dplyr)
  require(GenomicRanges)
  str = GR@ranges@start+1
  ed = GR@ranges@start+GR@ranges@width
  
  mat = read.delim(list.files(paste0("/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/dense/40k/cworld_TAD/rawIS/",sample.name,"/",
                          as.character(GR@seqnames),"/"),pattern = "insulation.bedGraph",full.names = T),skip = 1,header = F)
  boundary = read.delim(list.files(paste0("/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/dense/40k/cworld_TAD/rawIS/",sample.name,"/",
                                                as.character(GR@seqnames),"/"),pattern = "insulation.boundaries$",full.names = T),skip = 28,header = T)
  mat = mat %>% filter( V2>=str & V3<= ed)
  # boundary = boundary %>% filter( start >=str & start <= ed ) %>% dplyr::select(boundaryHeader)
  # boundary = gsub("bin\\w+\\|\\w+\\|\\w+:","",boundary$boundaryHeader) %>% as.data.frame()
  # boundary$start = apply(boundary,1,function(x){x = strsplit(x,"-")[[1]][1]})
  # boundary$end = apply(boundary,1,function(x){x = strsplit(x,"-")[[1]][2]})
  boundary = boundary %>% filter( start >=str & start <= ed )

  rects <- data.frame(start=as.numeric(boundary$start), end=as.numeric(boundary$end), group=seq_along(start))
  p = ggplot(mat ) + 
    geom_line( aes(x = V2,y = V4),size = 1 ) + 
    geom_rect(data=rects, 
              inherit.aes=FALSE, 
              aes(xmin=start, xmax=end, 
                  ymin=(-1), ymax=1, group=group), 
              color="transparent", fill="orange", alpha=0.3) +
    scale_y_continuous(limits = c(-1,1),breaks = c(-1,1),
                       expand = c(0,0)) +
    scale_x_continuous( breaks = c(str,ed),
                        labels = paste0(round(c((str-1)/1000000,(ed-1)/1000000),digits = 1),"Kb"),
                        expand = c(0,0))+
    geom_abline(slope = 0,intercept = 0,
                color = alpha("#e64047",alpha = 0.8),
                show.legend = T ,linetype="dashed",size = 1) + 
    theme_cowplot() +ylab("IS") +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(size = 1),
          axis.title.y = element_text(angle = 0,hjust = 0.5,vjust = 0.5,size = 20),
          axis.text.y = element_text(family = "sans",size = 20))
  p
  return(p)
}
  

### track for gene using Gviz----

# gtf = read.delim( "/lustre/user/liclab/fangq/ref/UCSC/mm10/Annotation/Genes/mm10.gene.bed",
#                   header = F,sep = "\t",col.names = c("chr","start","end","ncbi","strand","symbol") )
# tmp = GRanges(gtf %>% filter(chr =="chr5" & start>=147287555 & end <= 147340000))
# 
# tmp
# grtrack = GeneRegionTrack( tmp, 
#                            genome = "mm10",
#                            chromosome = "chr5",showId = T,geneSymbol = T,
#                            name = "gene model",start = 147290000,end = 147320000 )
# grtrack
# plotTracks(grtrack,
#                transcriptAnnotation = "symbol",
#                #collapseTranscripts = T,
#                col = "grey",fill = "grey",
#                fontfamily = "sans",fontfamily.group = "sans",
#                fontcolor.group = "black",
#                fontsize  =10,fontsize.group = 25,title.width = 0)


### track for AB compartment-----
# a$type = ifelse(a$V4>0,"A",ifelse(a$V4<0,"B",NA))
# ggplot(a,aes(x = V2,y = V4,fill = type,color = type)) + geom_bar(stat="identity")+ 
#   scale_x_continuous( expand = c(0,0) ) + scale_y_continuous( expand = c(0,0) ) + 
#   scale_color_manual(values = c("darkred","navy")) + scale_fill_manual(values = c("darkred","navy"))+
#   theme_null() +  guides(fill="none",color = "none") +
#   theme() 
