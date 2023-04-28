library(HiTC)


### pearson correlation-----
dot.corrplot = function(input.mat = matrix(),
                        is.tri = logical(),upper.lower = character(),
                        point.range = c(0,11),
                        cor.ranges = c(-1,1),cor.method = "pearson",
                        col.palette = colorRampPalette(c("#b55b41","white","#4f7e90"))(100),
                        breaks = numeric()){
  require(Matrix)
  require(ggplot2)
  require( data.table )
  tmp = as.matrix(cor(input.mat, method=cor.method))
  breaks = seq(floor(round(min(tmp),digits = 2)/0.05)*0.05,1,by = 0.05)
  if( is.tri ){
    if(upper.lower == "upper"){ tmp[upper.tri(tmp)]=NA
    }else if(upper.lower=="lower"){tmp[lower.tri(tmp)]=NA}
  }
  tmp = melt( tmp,value.name = "Cor")
  ggplot(tmp, aes(x = Var1, y = Var2)) + 
    geom_point(aes(color = Cor,size = abs(Cor),)) + 
    scale_y_discrete( limits = rev(levels(tmp$Var2)),position = "left") + 
    scale_x_discrete( position = "top" ) +
    scale_color_gradientn(  colours = col.palette,
                            values = cor.ranges,
                            breaks = breaks,
                            name = gsub(substr(cor.method,1,1),toupper(substr(cor.method,1,1)),cor.method))+
    theme_cowplot() + 
    scale_radius( range = point.range,breaks = breaks,name = "Correlation" )+
    theme( axis.ticks = element_blank(),
           axis.title = element_blank(),
           axis.line = element_blank(),
           axis.text = element_text( size = 24 ),
           text = element_text(size = 24),
           axis.text.x = element_text(angle = 30,hjust = 0))
}
heatmap.corrplot = function(input_data,
                              index,condition_col,
                              cellsize = 35,
                            cor.method = "pearson",
                            breaks = numeric(),
                              show_label = T){
  num = length(levels(index[,condition_col]))
  cor = as.matrix(cor(input_data, method=cor.method))
  color = colorRampPalette( brewer.pal(name = "Paired",n = 10) )(num)
  names(color) = levels(index[,which(colnames(index)==condition_col)])
  annotation_color = list(cell.type = color)
  breaks = seq(floor(round(min(cor),digits = 2)/0.05)*0.05,1,by = 0.05)
  bk = seq( min(cor),max(cor),length.out=100)
  bkcolor = c(colorRampPalette( c("#ffffff","#00607d"))(100))
  
  cor.heat = pheatmap::pheatmap(cor,clustering_method = "ward.D",
                                 annotation_col = column_to_rownames(index,"sample"),
                                 display_numbers = show_label,
                                 annotation_colors = annotation_color,
                                 fill = "white",border_color = NA,
                                 color = bkcolor,breaks = bk,legend_breaks = breaks,
                                 fontsize = 20,fontsize_number = cellsize/2-2,annotation_names_col = F,
                                 angle_col = 315,cellheight = cellsize,cellwidth = cellsize,
                                 treeheight_row = 18,treeheight_col = 18)
  print(cor.heat)
  return(cor.heat)
}
# heatmap.corrplot( PC1.COMB,index = index[,-3],condition_col = "cell.type",
#                   cor.method = "pearson",cellsize = 30,show_label = T )

## plot hic.matrix--------
plot.hic.matrix = function( mat,
                            matrix.type = "HiTC",
                            chr = c("chr1","chr1"),region = c(1,10),mat.resolution = 100000,ab.resolution = 500000,
                            wgchr.start = character(),
                            wgchr.end = character(),
                            AB.files = character(),
                            compartment.col = c("#e38d2e","#153e67"),
                            colorlist = c( "#e6424a","white","#1b378c" ),
                            quantile = 0.8,
                            name = character(),
                            wd = 8,hg = 7){
  if(matrix.type == "homer"){
    nr = mat[,1]
    tmp = mat[,-c(1:2)] %>% as.data.frame()
    rownames(tmp) = NULL
    rownames(tmp ) = nr$V1
    colnames(tmp) = nr$V1
    filter = nr$V1[-grep("chr[Y|M]-",nr$V1)] 
    tmp = tmp[filter ,c(filter)]
    if( length(chr) == 2 ){ 
      #tmp = as.matrix( mat[[paste0(chr[1],chr[2])]]@intdata)
      tmp = tmp[grep(paste0(chr[2],"-"),filter,value = T), grep(paste0(chr[1],"-"),filter,value = T)]
    }else if( chr == "WholeGenome" ){ 
      for (i in c(1:19,"X")) { tmp[ filter[grep(paste0("chr",i,"-"),filter)], filter[grep(paste0("chr",i,"-"),filter)] ] = NA}
    }else if( !is.null(wgchr.start) & !is.null(wgchr.end) ){
      if(wgchr.end=="X"){for (i in c(1:19,"X")) { tmp[ filter[grep(paste0("chr",i,"-"),filter)], filter[grep(paste0("chr",i,"-"),filter)] ] = NA}}
      else{
        idx = paste0("chr",as.integer(wgchr.start):as.integer(wgchr.end))
        for (i in idx) { tmp[ filter[grep(paste0(i,"-"),filter)], filter[grep(paste0(i,"-"),filter)] ] = NA}
        l = character()
        for(i in idx){ l = c(l,grep(paste0(i,"-"),filter,value = T,fixed = T)) }
        tmp = tmp[l,c(l)]
        
      }
    }
    
    ###--------AB compartment plot----------
    if(length(AB.files)!=0){
      ab = read.delim(AB.files,skip = 1,header=F)
      chr.bottom = ab[ab$V1==chr[1],c(1,2,4)]
      loc = grep("^\\d",unlist(strsplit(colnames(tmp),"-")),value = T)
      options( scipen = 200 )
      x = data.frame(V1=factor(chr[1],levels = levels(ab$V1)),
                     V2 = as.integer(setdiff(loc,chr.bottom$V2)),
                     V4 = 0 )
      chr.bottom = rbind(x,chr.bottom)
      
      
      chr.left =ab[ab$V1==chr[2],c(1,2,4)]
      loc = grep("^\\d",unlist(strsplit(rownames(tmp),"-")),value = T)
      options( scipen = 200 )
      x = data.frame(V1=factor(chr[2],levels = levels(ab$V1)),
                     V2 = as.integer(setdiff(loc,chr.left$V2)),
                     V4 = 0 )
      chr.left = rbind(x,chr.left)
      
      chr.bottom = chr.bottom %>% arrange( V2 )
      chr.left = chr.left %>% arrange( V2 )
      top = HeatmapAnnotation(foo = anno_barplot(chr.bottom$V4, height = unit(1, "cm"),
                                                 bar_width = 1,baseline = 0,
                                                 border = F,axis = F,extend = 0,
                                                 gp = gpar(fill = ifelse(chr.bottom$V4>0,compartment.col[1],
                                                                         ifelse(chr.bottom$V4<0,compartment.col[2],NA)),
                                                           col = NA )),
                              annotation_label = chr[1],annotation_name_side = "left"
      )
      left = rowAnnotation(foo = anno_barplot(chr.left$V4, width = unit(1, "cm"),
                                              bar_width = 1,baseline = 0,
                                              border = F,axis = F,extend = 0,
                                              gp = gpar(fill = ifelse(chr.left$V4>0,compartment.col[1],
                                                                      ifelse(chr.left$V4<0,compartment.col[2],NA)),
                                                        col = NA )),
                           annotation_label = chr[2],annotation_name_rot = 0,annotation_name_side = "top")
    }else { top=NULL ; left = NULL}
    
  }else if(matrix.type == "HiTC"){
    tmp = extractRegion( mat[[paste0(chr[1],chr[2])]],MARGIN = c(1,2),chr = chr[1],from = region[1],to = region[2] )
    tmp =  normPerExpected(tmp, method="loess")
    tmp = forceSymmetric(tmp)
    tmp = as.matrix(tmp@intdata)
    ###--------AB compartment plot----------
    if(length(AB.files)!=0){
      ab = read.delim(AB.files,skip = 1,header=F)
      chr.bottom = ab[ab$V1==chr[1],c(1,2,4)]
      loc.ab = seq(region[1],region[2]-1,by = ab.resolution)
      chr.bottom = chr.bottom[chr.bottom$V2 %in% loc.ab, ]
      loc.mat = seq(region[1],region[2]-1,by = mat.resolution)
      ab.resolution/mat.resolution
      options( scipen = 200 )
      chr.bottom = data.frame( chr = chr[1],
                               start = loc.mat,
                               pc1 = rep(chr.bottom$V4,each = ab.resolution/mat.resolution))
      
      
      ab = read.delim(AB.files,skip = 1,header=F)
      chr.left = ab[ab$V1==chr[1],c(1,2,4)]
      loc.ab = seq(region[1],region[2]-1,by = ab.resolution)
      chr.left = chr.left[chr.left$V2 %in% loc.ab, ]
      loc.mat = seq(region[1],region[2]-1,by = mat.resolution)
      ab.resolution/mat.resolution
      options( scipen = 200 )
      chr.left = data.frame( chr = chr[1],
                             start = loc.mat,
                             pc1 = rep(chr.left$V4,each = ab.resolution/mat.resolution))
      
      chr.bottom = chr.bottom %>% arrange( start )
      chr.left = chr.left %>% arrange( start )
      top = HeatmapAnnotation(foo = anno_barplot(chr.bottom$pc1, height = unit(1, "cm"),
                                                 bar_width = 1,baseline = 0,
                                                 border = F,axis = F,extend = 0,
                                                 gp = gpar(fill = ifelse(chr.bottom$pc1>0,compartment.col[1],
                                                                         ifelse(chr.bottom$pc1<0,compartment.col[2],NA)),
                                                           col = NA )),
                              annotation_label = chr[1],annotation_name_side = "left")
      left = rowAnnotation(foo = anno_barplot(chr.left$pc1, width = unit(1, "cm"),
                                              bar_width = 1,baseline = 0,
                                              border = F,axis = F,extend = 0,
                                              gp = gpar(fill = ifelse(chr.left$pc1>0,compartment.col[1],
                                                                      ifelse(chr.left$pc1<0,compartment.col[2],NA)),
                                                        col = NA )),
                           annotation_label = chr[2],annotation_name_rot = 0,annotation_name_side = "top")
    }else { top=NULL ; left = NULL}
  }
  
  tmp = as.matrix(tmp)
  limit.up = quantile(tmp,quantile,na.rm = T)
  tmp[tmp>limit.up] = limit.up
  
  col = colorRampPalette(rev(colorlist))(100)
  p = ComplexHeatmap::Heatmap( tmp,col = col,
                               na_col = "white",
                               heatmap_legend_param = list( title = NULL,
                                                            at = c(min(tmp,na.rm = T),0.5*max(tmp,na.rm = T),max(tmp,na.rm = T)),
                                                            labels = c("min","mean","90%max")),
                               cluster_columns = F,
                               cluster_rows = F,
                               show_heatmap_legend = T,
                               show_column_names = F,
                               show_row_names = F,
                               bottom_annotation = top,
                               left_annotation = left,
                               width = 8,height = 8)
  #print(p)
  return(p)
}

Slope.cutoff = function(value = numeric(),slope.cutoff = 1,plot = logical()){
  ma = max(value,na.rm = T)
  mi = min(value,na.rm = T)
  n = length(value)
  value = apply(data.frame(value),1,function(x){(x-mi)/(ma-mi)})
  value = data.frame(y = sort(value),x = c(1:n)/n)  ## sort and rank
  value$intercept = value$y-slope.cutoff*value$x
  cutoff.rank = rownames(value[which(value$intercept == min(value$intercept)),])
  if(plot){
    value$type = ifelse(as.numeric(rownames(value)) >= as.numeric(cutoff.rank) ,"Significant","Non-signif")
    p = ggplot( value, aes(x = x,y = y,color = type) ) + geom_point(size = 1) + 
      geom_abline( slope = slope.cutoff,intercept = min(value$intercept),
                   color = alpha("#e64047",alpha = 0.6),
                   show.legend = T ,linetype="dashed",size = 1) + 
      geom_hline( yintercept = value[rownames(value)==cutoff.rank,"y"],
                  color = alpha("#42776a",alpha = 0.8),size = 1,
                  linetype="dashed" ) +
      scale_color_manual(values = alpha(c("grey","darkred"),alpha = 0.5),name = "Gene Type") + theme_cowplot()
    plot(p)
  }
  return(as.numeric(cutoff.rank))
}

library(ggsci) ### p value plot

# stats = diff %>% 
#   #mutate(log.expr = log10(`ST-HSC`+1)) %>% 
#   t_test( log2FoldChange ~ GAD.type,
#           paired = F,
#           conf.level = 0.05,
#           p.adjust.method = "BH",
#           ref.group = "stable" ) %>% add_xy_position(x = "GAD.type")
# 
# ggplot(diff) + geom_boxplot( aes( x = GAD.type, y = log2FoldChange,fill = GAD.type),outlier.size = 1 ) +  
#   #scale_y_continuous( limits = c(0,3) )+
#   stat_pvalue_manual(stats,
#                      label = "p.adj.signif",
#                      y = 12.5,
#                      remove.bracket = T,
#                      #bracket.size = 0.5,
#                      x = "group2",
#                      linetype = 1,hide.ns = T,label.size = 10,
#                      bracket.nudge.y = 0.1,
#                      tip.length = 0,
#                      vjust = 1 ) + 
#   geom_hline(yintercept = 0,color = alpha("darkred",alpha = 0.8),
#              show.legend = T ,linetype="dashed",size = 1 )+
#   scale_fill_manual(values = alpha(c("grey",pal_npg(palette = c("nrc"))(4)[c(1,4)]),alpha = 1),name = "GAD score") +
#   theme_cowplot() +
#   xlab("GAD score(GR/ST)") + ylab( "log2FoldChange(GR/ST)" ) + 
#   theme( text = element_text( size = 20 ),axis.text = element_text(size = 18) )

GO_Plot = function( GO_dataset, padj_threshold=0.05, show_num = 15, 
                    fill_color = "#e63f46",wid = 10, hgt = 8,
                    output_prefix,
                    GO_term = character(),
                    keywords = character()){
  tmp = GO_dataset@result
  if(length(keywords)!=0){
    seek = numeric()
    for( i in keywords){
      j = grep(i,tmp$Description)
      seek = c(seek,j)
    }
    if(length(seek)==0){return("No aim GO terms was found!")}
    tmp = tmp[seek,]
  }
  if(!is_empty( GO_term )){
    tmp = tmp %>% dplyr::filter( ONTOLOGY %in% GO_term )
    tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
    tmp$logpval = (-log10( tmp$pvalue ))
    tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
    tmp = tmp[1:min(show_num,nrow(tmp)),] 
    tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
    
    go = ggplot()+
      geom_bar( data=tmp,
                aes_string(x = "Description" ,y = "logpval"),
                fill = fill_color,
                stat = "identity",width = 0.5) + ylab( "-log(p.value)" ) +
      coord_flip( )+
      facet_grid( drop = T,rows = "ONTOLOGY",scales = "free",space = "free_y")+
      theme_cowplot()+
      theme( text = element_text( size = 28),
             axis.title = element_text(size = 28),
             axis.text= element_text(size = 15,family = "sans"),
             line = element_line(size = 1),
      )
  }else{
    tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
    tmp$logpval = (-log10( tmp$pvalue ))
    tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
    tmp = tmp[1:min(show_num,nrow(tmp)),] 
    tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
    
    go = ggplot()+
      geom_bar( data=tmp,
                aes_string(x = "Description" ,y = "logpval"),
                fill = fill_color,
                stat = "identity",width = 0.5) + ylab( "-log(p.value)" ) +
      coord_flip( )+
      theme_cowplot()+xlab(NULL)+
      theme( text = element_text( size = 28),
             axis.title = element_text(size = 28),
             axis.text= element_text(size = 15,family = "sans"),
             line = element_line(size = 1),
      )
    go
    #ggsave( go,filename = paste0( output_prefix,".pdf" ),width = wid, height = hgt, units = "in" )
  }
}

#### mannual exchange p value to label----------
p.label = function(p.value){
  for(i in 1:length(p.value)){
    if(is.na(p.value[i])){p.value[i] = NA
    }else if(p.value[i]>0.05){p.value[i] = NA
    }else if(p.value[i]<=0.05&p.value[i]>0.01){p.value[i] = c("*")
    }else if(p.value[i]<=0.01&p.value[i]>0.001){p.value[i] = c("**")
    }else if(p.value[i]<=0.001&p.value[i]>0.0001){p.value[i] = c("***")
    }else{p.value[i] = c("****")
    }
  }
  return(p.value)  
}


### triangle plot------
chr10 = extractRegion( GMP$chr17chr17,c(1,2),chr = "chr17",from = 15000000,to = 25000000 )
# library(InteractionSet)
# tmp = InteractionSet::GInteractions(chr10@xgi,chr10@ygi,as.matrix(chr10@intdata))
# tmp = GInteractions(chr10@xgi, chr10@ygi,regions = )
tmp = as.matrix(chr10@intdata)
tmp[lower.tri(tmp)] = NA
tmp = melt( tmp,value.name = "Contact")
tmp$Contact[tmp$Contact>quantile(tmp$Contact,0.9,na.rm = T)] = quantile(tmp$Contact,0.9,na.rm = T)
max(tmp$Contact,na.rm = T)
min(tmp$Contact,na.rm = T)
summary(tmp$Contact)
min(tmp$Var1)
min(tmp$Var2)
theta = cos(45)
tmp$x = (tmp$Var1-min(tmp$Var1))*theta+(tmp$Var2-min(tmp$Var2))*theta
tmp$y = (tmp$Var2-min(tmp$Var2))*theta-(tmp$Var1-min(tmp$Var1))*theta


p = ggplot(tmp, aes(x = x, y = y)) +
  # geom_tile(aes(fill = Contact)) +
  geom_point( aes(color = Contact),shape = 18,size = 2.7 ) +
  # scale_fill_continuous(type = "viridis") +
  scale_y_discrete( limits = rev(levels(sqrt(2)/2*tmp$Var2)),position = "left",expand = c(0,0)) +
  scale_x_discrete( position = "top",expand = c(0,0) ) +
  # geom_abline( slope = 1,intercept = -1,
  #              color = "white",
  #              show.legend = T ,linetype="solid",size = 3) +
  scale_color_gradientn(  colours = colorRampPalette(c("white","red","black"))(1000),
                          na.value = "white",
                          # values = c(0,floor(max(tmp$Contact))),
                          breaks = seq(0,round(max(tmp$Contact,na.rm =T),digits = 0),length.out = 5))+
  theme_cowplot() +
  theme( axis.ticks = element_blank(),
         axis.title = element_blank(),
         # axis.line = element_blank(),
         axis.text = element_text( size = 24 ),
         text = element_text(size = 24),
         axis.text.x = element_text(angle = 30,hjust = 0))
p

####genova plot------------
source("~/ffqq/proj/R_pack/trackplot.R")
# library(GENOVA,lib.loc = "~/ffqq/proj/R_pack/GENOVA/")
library( HiTC )
library(RColorBrewer)
library(tidyverse)
library(GENOVA)
ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79


genova_ko1 <- load_contacts(
  signal_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-KO-1/iced/10000/E95-KO-1.ds_10000_iced.matrix',
  indices_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-KO-1/raw/10000/E95-KO-1.ds_10000_abs.bed',
  sample_name = "E95-KO-1", #centromeres = centromeres,
  colour = "blue"
)

genova_ko2 <- load_contacts(
  signal_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-KO-2/iced/10000/E95-KO-2.ds_10000_iced.matrix',
  indices_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-KO-2/raw/10000/E95-KO-2.ds_10000_abs.bed',
  sample_name = "E95-KO-2", #centromeres = centromeres,
  colour = "green"
)
genova_wt1 <- load_contacts(
  signal_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-WT-1/iced/10000/E95-WT-1.ds_10000_iced.matrix',
  indices_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-WT-1/raw/10000/E95-WT-1.ds_10000_abs.bed',
  sample_name = "E95-WT-1", #centromeres = centromeres,
  colour = "black"
)

genova_wt2 <- load_contacts(
  signal_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-WT-2/iced/10000/E95-WT-2.ds_10000_iced.matrix',
  indices_path = '/lustre/user/liclab/fangq/proj/embryo/hic/matrixs/E95-WT-2/raw/10000/E95-WT-2.ds_10000_abs.bed',
  sample_name = "E95-WT-2", #centromeres = centromeres,
  colour = "grey"
)

save( genova_ko1,genova_ko2,genova_wt1,genova_wt2,
      file = "/lustre/user/liclab/fangq/proj/embryo/hic/analysis/TAD/genova.RData",compress = T )
load("/lustre/user/liclab/fangq/proj/embryo/hic/analysis/TAD/genova.RData",verbose = T)
ppp = c("#7C65A9","#806AA7","#846FA5","#8874A2","#8D78A0","#917D9E",
        "#95829C","#99879A","#9D8C97","#A19195","#A59693","#AA9A91",
        "#AE9F8F","#B2A48C","#B6A98A","#BAAE88","#BEB386","#C2B883",
        "#C6BD81","#CBC17F","#CFC67D","#D3CB7B","#D7D078","#DBD576",
        "#DFDA74","#E3DF72","#E8E370","#ECE86D","#F0ED6B","#F4F269")
ppp = colorRampPalette(c("#432371","#F4F269"))(50)
options("GENOVA.colour.palette" = ppp)
library(GENOVA)


GR = GRanges("chr5:146000000-149000000")
gr = GRanges("5:146000000-149000000")
s = 1460e5
e = 1490e5
chr = "chr5"
scale_limits = 30
hgts = 10e5

ko1.pyr <- pyramid(
  genova_ko1,
  chr = chr, start = s, end = e,
  crop_y = c(0, hgts),
  colour = c(0,scale_limits)
) + ggplot2::ggtitle("KO1") +
  theme(text = element_text(size = 22),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.direction = "horizontal",
        legend.position = "none")

ko2.pyr <- pyramid(
  genova_ko2,
  chr = chr,start = s, end = e,
  crop_y = c(0, hgts),
  colour = c(0,scale_limits)
) + ggplot2::ggtitle("KO2") +
  theme(text = element_text(size = 22),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.direction = "horizontal",
        legend.position = "none")
wt1.pyr <- pyramid(
  genova_wt1,
  chr = chr, start = s, end = e,
  crop_y = c(0, hgts),
  colour = c(0,scale_limits)
) + ggplot2::ggtitle("WT1") +
  theme(text = element_text(size = 22),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.direction = "horizontal",
        legend.position = "none")
wt2.pyr <- pyramid(
  genova_wt2,
  chr = chr, start = s, end = e,
  crop_y = c(0, hgts),
  colour = c(0,scale_limits)
) + ggplot2::ggtitle("WT2") +
  theme(text = element_text(size = 22),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.direction = "horizontal",
        legend.position = "none")

library(patchwork)
 
# + plot_layout(widths = unit(c(20,20,20,20), c('cm','cm','cm','cm')),
#                                                         heights = unit(c(12,2,2,2), c('cm','cm','cm','cm')),
#                                                         guides = "collect" )
# 

# gmp.k4me1 = "/lustre/user/liclab/fangq/proj/tagHiC_rep/ChIP-seq/GSE60103/mm10fq/bw/H3K4me1/GSM1441289_H3K4me1GMP.ucsc.10.bw"
# g.k4me1 = "/lustre/user/liclab/fangq/proj/tagHiC_rep/ChIP-seq/GSE60103/mm10fq/bw/H3K4me1/GSM1441293_H3K4me1GN.ucsc.10.bw"
# gmp.k27ac = "/lustre/user/liclab/fangq/proj/tagHiC_rep/ChIP-seq/GSE60103/mm10fq/bw/H3K27ac/GSM1441273_K27Ac_GMP.ucsc.10.bw"
# g.k27ac = "/lustre/user/liclab/fangq/proj/tagHiC_rep/ChIP-seq/GSE60103/mm10fq/bw/H3K27ac/GSM1441277_K27Ac_GN.ucsc.10.bw"
# 
# gmp.rna = "/lustre/user/liclab/fangq/proj/tagHiC_rep/RNA-seq/STAR_results/bam/gmp.bw"
# g.rna = "/lustre/user/liclab/fangq/proj/tagHiC_rep/RNA-seq/STAR_results/bam/gr.bw"
# 
gene.track = plot.gene( database = ensdb,genome.range = gr,gene.col = "#005f77",gene.list = c("Cdx2"))
ko1.is = plot_isulation( GR,sample.name = "E95-KO-1" )
ko2.is = plot_isulation( GR,sample.name = "E95-KO-2" )
wt1.is = plot_isulation( GR,sample.name = "E95-WT-1" )
wt2.is = plot_isulation( GR,sample.name = "E95-WT-2" )
pdf("TAD/chr5-1465-1490.pdf",width = 20,height = 10)
((ko1.pyr/ ko1.is/ gene.track@ggplot )|( ko2.pyr / ko2.is/ gene.track@ggplot)) / 
  ( (wt1.pyr /wt1.is/ gene.track@ggplot) | (wt2.pyr/wt2.is/ gene.track@ggplot)) 
dev.off()
# gmp.t1 = plot.track( GR,bigwigfile = gmp.k4me1,y.axis = c(0,50),log.scale = F,downsample = T,
#                      smooth = 1000,plot.type = "coverage",track.col = "#065983",
#                      y.label = "H3K4me1")
# gmp.t2 = plot.track( GR,bigwigfile = gmp.k27ac,y.axis = c(0,50),log.scale = F,downsample = T,
#                      smooth = 1000,plot.type = "coverage",track.col = "#c298be",
#                      y.label = "H3K27ac")
# gmp.t3 = plot.track( GR,bigwigfile = gmp.rna,y.axis = c(0,1000),log.scale = F,downsample = T,
#                      smooth = 1000,plot.type = "coverage",track.col = "#ef476f",
#                      y.label = "RNA-seq")
# 
# 
# gr.t1 = plot.track( GR,bigwigfile = g.k4me1,y.axis = c(0,50),log.scale = F,downsample = T,
#                     smooth = 1000,plot.type = "coverage",track.col = "#065983",
#                     y.label = "H3K4me1")
# gr.t2 = plot.track( GR,bigwigfile = g.k27ac,y.axis = c(0,50),log.scale = F,downsample = T,
#                     smooth = 1000,plot.type = "coverage",track.col = "#c298be",
#                     y.label = "H3K27ac")
# gr.t3 = plot.track( GR,bigwigfile = g.rna,y.axis = c(0,1000),log.scale = F,downsample = T,
#                     smooth = 1000,plot.type = "coverage",track.col = "#ef476f",
#                     y.label = "RNA-seq")
# 
# pdf("~/ffqq/proj/tagHiC_rep/r-code/tagHi-C_AB/test.pdf",width = 30,height = 15)
# (gmp.pyr /gmp.t1/gmp.t2 /gmp.t3 /gene.track@ggplot) + plot_layout(widths = unit(c(20,20,20,20,20), c('cm','cm','cm','cm','cm')),
#                                                                   heights = unit(c(12,2,2,2,4), c('cm','cm','cm','cm','cm')),
#                                                                   guides = "collect" )| 
#   ( gr.pyr /gr.t1/gr.t2/ gr.t3/gene.track@ggplot ) + plot_layout(widths = unit(c(20,20,20,20,20), c('cm','cm','cm','cm','cm')),
#                                                                  heights = unit(c(12,2,2,2,4), c('cm','cm','cm','cm','cm')),
#                                                                  guides = "collect" )
# dev.off()

###### saddle plot------
library(GENOVA)
library(tidyverse)
# st.test = load_contacts(
#   signal_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/iced/ST/ST_500000.iced.matrix',
#   indices_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/raw/ST/raw/500000/ST_500000_abs.bed',
#   sample_name = "ST-HSC", #centromeres = centromeres,
#   colour = "blue"
# )
# mpp.test = load_contacts(
#   signal_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/iced/MPP/MPP_500000.iced.matrix',
#   indices_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/raw/MPP/raw/500000/MPP_500000_abs.bed',
#   sample_name = "MPP", #centromeres = centromeres,
#   colour = "orange"
# )
# gmp.test = load_contacts(
#   signal_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/iced/GMP/GMP_500000.iced.matrix',
#   indices_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/raw/GMP/raw/500000/GMP_500000_abs.bed',
#   sample_name = "GMP", #centromeres = centromeres,
#   colour = "black"
# )
# gr.test = load_contacts(
#   signal_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/iced/G/G_500000.iced.matrix',
#   indices_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/raw/G/raw/500000/G_500000_abs.bed',
#   sample_name = "GR", #centromeres = centromeres,
#   colour = "red"
# )
# mk.test = load_contacts(
#   signal_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/iced/MK/MK_500000.iced.matrix',
#   indices_path = '/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/hicpro/merged.matrix/raw/MK/raw/500000/MK_500000_abs.bed',
#   sample_name = "MK", #centromeres = centromeres,
#   colour = "green"
# )
# library(data.table)
# bed = fread( "/lustre/user/liclab/fangq/proj/tagHiC_rep/HiC/homer_merged_from_aVP/PC1.500k/GMP.500k.PC1.txt",header = F,nThread = 6,skip = 1)
# bed = bed[bed$V6>0]
# bed = bed[,c(2,3,4)]

library(cowplot)
library(tidyverse)

# cs =  compartment_score( list( st.test,mpp.test,
#                                gmp.test,mk.test,gr.test ),
#                          bed = bed )
# cs$compart_scores$GMP
# visualise(cs,chr = "chrX")
# saddle.out = saddle( explist = list( st.test,mpp.test,
#                                      gmp.test,mk.test,gr.test ),
#                      CS_discovery = cs,
                     # bins = 50)
saddle.fq <- function( saddle.vactor ){
  require(GENOVA)
  require(ggplot2)
  p = visualise(saddle.out,raw = T) + 
    theme_cowplot() + xlab("Quantile") + ylab("Quantile") + 
    scale_altfill_continuous(low = "#2c6e49",high = "#941b0c",
                             name = "Difference",) + 
    scale_fill_gradientn(colours = c( colorRampPalette(rev(c("#ffffff","#407ba7","#004e89","#002962","#00043a")))(100),
                                      colorRampPalette(rev(c("#800016","#a0001c","#c00021","#ff002b","#ffffff")))(500)),
                         name = "Observed/Expected") + 
    scale_x_continuous(expand = c(0,0),breaks = c(0,max(p$data$q1)),labels = c("B","A")) +
    scale_y_continuous(expand = c(0,0),breaks = c(0,max(p$data$q2)),labels = c("B","A")) +
    coord_fixed(ratio = 1)+
    theme(
      # axis.title = element_blank(),
      axis.line = element_blank(),
      text = element_text(size = 18),
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      strip.background = element_rect(fill = "white",colour = "white"),
      strip.placement = "outside",legend.position = "top",
      panel.border = element_rect(color = "black"),
    )
  return(p)
}
# p = saddle.fq( saddle.out )
# p
