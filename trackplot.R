#library(Gviz) #require for IGrange/GenomicRanges/GenomeInforDb/grid
require(GenomicRanges)
require(dplyr)
require(patchwork)
require(rtracklayer)
require(tidyverse)

# library(ggbio)
# library(biovizBase)
require(ggplotify)

#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(AnnotationFilter)

require(cowplot)
require(scales)
require(ensembldb)
require(stringr)
require(ggplot2);require(cowplot);require(ggarchery)
require(dplyr)
ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79

###----functions ------
plot_track = function(region = GRanges(),
                      bigwigfile = character(),
                      downsample = F,
                      log.scale = F,
                      smooth = F,
                      max.downsample = 5000,
                      downsample.rate = 0.1,
                      smooth.rate = 200,
                      track.col = "red",
                      plot.type = "coverage",
                      heat.bin = 20,
                      y.axis = c(0,100),
                      y.label = character(),
                      seed = 1234){
  if( downsample ){
    window.size = width( region)
    sampling = min(max.downsample, floor(window.size * downsample.rate))
    set.seed(seed)
    tmp = slice_sample(.data = tmp, n = sampling)
  }
  chr = rtracklayer::import(con = bigwigfile,format = "bw", as = "NumericList",which = region )[[1]]
  # chr = RcppRoll::roll_mean(x = chr, n = smooth, fill = 0L) # smooth line
  
  tmp = data.frame(
    position = start(region):end(region),
    score = chr,
    stringsAsFactors = FALSE)
  if( downsample ){
    window.size = width( region)
    sampling = min(max.downsample, floor(window.size * downsample.rate))
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
    a = a + geom_ribbon(aes(ymax=score,ymin=0),fill=track.col,color=track.col,position = 'identity')
    # a = a + geom_area(fill=track.col,color=track.col)
  }else if(plot.type=="line"){
    a = a + geom_line()
  }else if(plot.type=="bar"){
    a = a + geom_bar(stat="identity",fill=track.col,color=track.col)
  }
  a = a + theme_cowplot() + 
    ylab(y.label)+
    scale_y_continuous( limits = y.axis,expand = c(0,0),breaks = y.axis[2],
                        labels = paste0("[",y.axis[1],"-",y.axis[2],"]"),position = 'left' ) +
    scale_x_continuous( limits = c(start(region),end(region)),expand = c(0,0) ) +
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          axis.title.y = element_text(angle = 0,hjust = 1,vjust = 0.5),
          axis.ticks = element_blank(),axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_text(size = 10,hjust = 2),
          text = element_text(family = "sans",size = 20))
  a
  return(a)
}
####----**** This is the true trackplot ****----####
####----**** NOTE: only support gencode GTF file up to now ****----####
# gtf = rtracklayer::import('/lustre/user/liclab/fangq/ref/UCSC/mm10/Annotation/gencode.vM10.annotation.gtf')
plot_gene = function(a = Granges(),
                     gene_widths_proportion = 0.1,
                     gtf = data.frame()){
  if(purrr::is_empty(gtf)){
    stop("No GTF file supported!!!")
  }
  aa = gtf[gtf@seqnames %in% levels(a@seqnames)  & end(gtf@ranges)>=start(a) & start(gtf@ranges)<=end(a) ,]
  tmp = as.data.frame(aa)
  tmp$type=factor(tmp$type,
                  levels = rev(c('start_codon',"stop_codon",'exon',"CDS","UTR",'transcript','gene','Selenocysteine')),
                  ordered = T)
  
  tmp$start = ifelse(tmp$start<start(a),start(a),tmp$start)
  tmp$end = ifelse(tmp$end>end(a),end(a),tmp$end)
  tmp = tmp %>% dplyr::group_by(gene_name) %>% 
    dplyr::mutate(tx_num = length(unique(transcript_id))) %>%
    as.data.frame()
  ### rank the y position in grouped gene list
  xx = tmp[tmp$type=='gene',c('gene_name','start','end','tx_num')] %>% unique()
  if(nrow(xx)==1){
    group = xx
    group$rank=0
  }else{
    group = as.list(1:nrow(xx))
    for(g in 1:length(group)){
      group[[g]]=xx[1,]
      for(i in 2:nrow(xx)){
        j1=group[[g]][nrow(group[[g]]),]
        j=xx[i,]
        if(j$start>=j1$end){
          group[[g]]=rbind(group[[g]],j)
        }
      }
      group[[g]]$rank=g-1
      xx=xx[!(xx$gene_name %in% group[[g]]$gene_name),]
      if(nrow(xx)<=1){
        if(nrow(xx)==1){
          group[[g+1]]=xx
          group[[g+1]]$rank=g
          for(l in g+2:length(group)){
            group[[l]]=NA
          }
          break
        }
        else{
          for(l in g+1:length(group)){
            group[[l]]=NA 
          }
          break
        }
      }
    } 
  }
  group = do.call(rbind,group)%>% drop_na() 
  group = group %>% dplyr::group_by(rank) %>% mutate(avg = quantile(tx_num,0.7))
  
  for(i in unique(group$rank)){
    now_rank=i
    if(now_rank!=0){
      above_rank=i-1
      value = unique(group$avg[group$rank==above_rank])
      group$rank[group$rank==now_rank]=value*0.5+now_rank
    }
    
  }
  tt = tmp[tmp$type=='transcript',c('gene_name','start','end','width','transcript_id','transcript_name')]
  tt = merge(tt,group[,c('gene_name','rank')],by='gene_name')
  tt = tt %>% arrange(start,desc(width))
  tt  =tt %>% group_by(gene_name) %>% mutate(tx_rank=(-1)*(rank+0.5*0:(length(transcript_id)-1)))
  
  k = merge(tmp,tt[,c('transcript_id','rank','tx_rank')],by='transcript_id')
  k = k %>% group_by(gene_name) %>% mutate(tx_num = length(unique(transcript_id)))
  k = k %>% group_by(rank) %>% mutate(avg_tx = quantile(unique(tx_num),0.75) )
  k$y_mid=k$tx_rank
  k$y_mid[k$rank==1] = k$y_mid[k$rank==1] -unique(k$avg_tx[k$rank==0])
  k = k %>% dplyr::group_by(gene_name) %>% mutate(y_min = y_mid-gene_widths_proportion*1,y_max = y_mid+gene_widths_proportion*1,
                                                  y_min2=y_mid-gene_widths_proportion/2*1, y_max2=y_mid+gene_widths_proportion/2*1)
  k$tx_x = ifelse(k$strand=='-',k$end,k$start)
  k$tx_xend = ifelse(k$strand=='-',k$start,k$end)
  
  
  cols = c('grey70','black','black','grey40','red','red')
  names(cols) = c('transcript','exon',"CDS","UTR",'start_codon',"stop_codon")
  p = ggplot() +
    geom_rect(data=k[k$type %in%c("CDS"), ],
              aes(xmin=start,xmax=end,ymin=y_min,ymax=y_max,fill=type,group=gene_name)) +
    geom_rect(data=k[k$type %in%c("UTR",'exon'),],
              aes(xmin=start,xmax=end,ymin=y_min2,ymax=y_max2,fill=type,group=gene_name)) +
    geom_segment(data=k[k$type %in%c('start_codon',"stop_codon"),],
                 aes(x=(start+end)/2,xend=(start+end)/2,y=y_min,yend=y_max,group=gene_name),
                 color = 'red',linewidth=unit(0.5,'cm')) +
    # geom_segment(data=k[k$type %in%c('transcript'),],
    #             aes(x=tx_x,xend=tx_x,y=y_min,yend=y_max,group=gene_name),
    #             color = 'black',
    #             linewidth=unit(0.7,'cm')) +
    scale_fill_manual(values = cols) +
    ggarchery::geom_arrowsegment(data=k[k$type %in%c('transcript'),],
                                 aes(x=tx_x, y=y_mid, xend=tx_xend, yend=y_mid,group=gene_name),
                                 arrow_positions = seq(0.1,1,0.2),
                                 arrows = arrow(length=unit(0.1, "cm"),angle = 30)) +
    scale_x_continuous(limits = c(start(a),end(a)),expand = c(0,0),
                       breaks = c(start(a),end(a)),label=NULL) +
    scale_y_continuous(limits = c(min(k$y_mid)-1,1),label=NULL) +
    geom_text(data=k[k$type %in%c('transcript'),],
              aes( x=(start+end)/2,y = y_mid-gene_widths_proportion*2.3,label=gene_name),
              size=3.5,nudge_x = 0.5) +
    theme_cowplot() + xlab(paste0(levels(a@seqnames),":",start(a),"-",end(a))) + ylab("Genes") +
    theme(
      axis.line = element_blank(),axis.ticks=element_blank(),legend.position="NULL",
      axis.title.y = element_text(angle = 0,hjust = 1,vjust = 0.5,size=20)
    )
  return(p)
}
####----**** This is the true trackplot ****----####

plot_insulation = function( GR = GRanges(),
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

# plot.gene = function(database = NULL,
#                      names.expr = "gene_name",
#                      genome.range = GRanges(),
#                      expand.for.symbol = c(0,2.5),
#                      symbol.size = 4,
#                      gene.col = "#2d266c",
#                      gene.list = character(),
#                      label = F){
#   if(!is_empty(gene.list)){
#     gtf = read.delim( "/lustre/user/liclab/fangq/ref/UCSC/mm10/Annotation/Genes/mm10.gene.bed",
#                       header = F,sep = "\t",col.names = c("chr","start","end","ncbi","strand","symbol") )
#     tmp = GRanges(gtf %>% dplyr::filter(chr ==paste0("chr",seqnames(genome.range)) & start>=start(genome.range) & end <= end(genome.range)) %>%
#                     dplyr::filter(symbol %in% gene.list ))
#     wh = range(tmp,ignore.strand = T)
#     b = autoplot(database,GeneNameFilter(gene.list),
#                  stat = "reduce",
#                   names.expr = names.expr,label = label,
#                   axis.text.x = NULL,layout = "linear",coord = "genome",#geom = "full",
#                   label.color = "black",color = gene.col,fill = gene.col,#,#(collapse)
#                   label.size = symbol.size,indel.color = "red") +
#     scale_x_continuous( limits = c(start(genome.range),end(genome.range)),expand = c(0,0) ) +
#     scale_y_continuous( expand = expand.for.symbol ) +
#     theme_cowplot() +
#     xlab(paste0("Chr",seqnames(genome.range)," : ",start(genome.range )," - ",end(genome.range))) +
#     ylab("Gene") + 
#     theme(axis.title.x = element_text(size = 22),
#           axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
#           axis.title.y = element_text(size = 20,angle = 0,vjust = 0.5,hjust = 1))
#   b
#   return(b)
#   }
#   b = autoplot(database, 
#                which = genome.range, 
#                axis.text.x = NULL,
#                label.color = "black",
#                color = gene.col,
#                fill = gene.col,label = label,
#                names.expr = names.expr,
#                label.size = symbol.size,
#                indel.color = "red",
#                label.repel = T) +
#     scale_x_continuous( limits = c(start(genome.range),end(genome.range)),expand = c(0,0) ) +
#     scale_y_continuous( expand = expand.for.symbol ) +
#     theme_cowplot() +
#     xlab(paste0("Chr",seqnames(genome.range)," : ",start(genome.range )," - ",end(genome.range))) +
#     ylab("Gene") + 
#     theme(axis.title.x = element_text(size = 22),
#           axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
#           axis.title.y = element_text(size = 20,angle = 0,vjust = 0.5,hjust = 1))
#   b
#   return(b)
# }


  
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
