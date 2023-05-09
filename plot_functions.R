drop_Nas = function (data,dims = 1,partitial = F,valid_proportions = as.numeric()){
  if (dims == 1){
    if(partitial){
      keep_rows = apply(data,1,function(x){count(x,value = NA)/length(x)<=(1-valid_proportions)}) ### keep rows that contains at least one non-NA
    }else{
      keep_rows = apply(data,1,function(x){any(!is.na(x))}) ### keep rows that contains at least one non-NA
    }
    data = data[keep_rows,]
  }
  else if (dims == 2){
    if(partitial){
      keep_cols = apply(data,2,function(x){count(x,value = NA)/length(x)<=(1-valid_proportions)}) ### keep rows that contains at least one non-NA
    }else{
      keep_cols = apply(data,2,function(x){any(!is.na(x))})
    }
    data = data[,keep_cols]
  }
  else{
    warning("dims can only equal to 1 and 2, 1=row, 2=column ")
  }
  return(data)
}

### plot functions ###
###***--- For geneal plot themes ---***###
Themes = function(pt_size = 8,
                  fontsize = 28,
                  axis_x_fontsize = 28,
                  axis_y_fontsize = 28,
                  titlefontsize = 28,
                  linesize = 2,
                  rotate = T){
    require(cowplot)
    require(ggplot2)
    a = theme_cowplot(font_family = "sans")+
            theme(  text = element_text( size = fontsize),
                    axis.ticks = element_line(linewidth = linesize),
                    axis.text.x = element_text( size = axis_x_fontsize),
                    axis.text.y = element_text( size = axis_y_fontsize),
                    axis.title = element_text(size = fontsize,family = "sans"),
                    axis.line = element_line(linewidth = linesize),
                    strip.background = element_blank(),
                    strip.text = element_text(size = titlefontsize,face = "bold"),
                    plot.title = element_text(hjust = 0.5,size = titlefontsize),
                    legend.position = "right")
    if(rotate){
      a = a+ theme(axis.text.x = element_text(angle = 30,hjust = 1))
    }
    return(a)
}


###***--- For statistics test plotting ---***###
data_summary <- function(x) { 
  # m <- mean(x)
  m = median(x)
  ymin <- m-sd(x) 
  ymax <- m+sd(x) 
  return(c(y=m,ymin=ymin,ymax=ymax)) 
}
# stat_summary(fun.data=data_summary,size = 1.2)

# stats_p2 = na.omit(ungroup(mat2)) %>% 
#           rstatix::anova_test(MeanSCR~Carrier.x) %>% 
#           # rstatix::adjust_pvalue() %>% #rstatix::add_xy_position(dodge = 0.8) %>%  
#           rstatix::add_significance(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

# stat_pvalue_manual(data = as.data.frame(ttest), label="p.signif",
#                         size = 0.2,bracket.size = 0.8,step.increase = 0,
#                        tip.length=0.02,label.size = 6,family="sans") +
#     stat_summary(aes( x = Var1, y = Freq, group = label ),fatten =0,
#                  fun.data=data_summary,size = 1, position = position_dodge(width = 0.8)) +


# nikolai1 = ggplot( tab,aes(x = Filtration, y = PSM_Number )) + 
#   geom_violin(aes(fill = Filtration,color = Filtration),
#               na.rm = T,size = 1,width = 0.6) +#, outlier.size = 0) +
#   scale_fill_manual( values = alpha(colorRampPalette(c("navy","#8acef7"))(4),0.85 )) + 
#   scale_color_manual( values = alpha(colorRampPalette(c("navy","#8acef7"))(4),0.8) ) + 
#   geom_text( aes( y = 0.95*median, x = Filtration,group = Filtration,label = median ),
#              check_overlap = F,size = 7,hjust = 0.5,vjust = -6,
#   ) +
#   stat_summary(fun = mean,size = 4,color = "white",
#                geom = "point" ) + 
#   stat_summary(fun.data = data_summary,size = 1,color = "white",
#                geom = "errorbar",width = 0.1 ) + 
#   Themes(pt_size = 8,alpha = 1,fontsize = 28,linesize = 2) + 
#   scale_x_discrete( labels=c("Raw","Contaminants","PIF!=NA","PIF>0.5") ) + 
#   # scale_y_continuous( expand = c(0.05,-0.2)) + 
#   xlab(NULL) + ylab( "PSM Number" ) +
#   theme( axis.text.x = element_text(angle = 30,hjust = 1),
#          panel.background = element_rect( color = NA ),
#          strip.background = element_rect(fill = NA),
#          strip.text = element_text(size = 28,family = "sans",face = "bold.italic"),
#          legend.position = "none") 

###***--- For PCA plotting ---***###
PCA_point = function(data = data.frame(),
                     Group_shape = character(),
                     Group_color = character(),
                     deviation = numeric(),
                     color = character(),
                     shape = c(16:21),
                     pt_size = 8,alpha = 1,
                     fontsize = 28,linesize = 2
                     ){
    require(rlang)
    require(ggplot2)
    require(cowplot)
    if(!rlang::is_empty(Group_shape)){
        p = ggplot() +
            geom_point( data = data,aes_string(x = "PC1", y = "PC2", color = Group_color, shape = Group_shape),
                        size = pt_size ) +
            # xlab(paste0("PC1 (",round(deviation[1],1),"% )")) +
            # ylab(paste0("PC2 (",round(deviation[2],1),"% )")) +
            scale_color_manual( values = alpha(color ,alpha = alpha) )+ 
            scale_shape_manual(values = c(16,17,18,19,20,21) ) + theme_cowplot()+
            theme( text = element_text( size = fontsize),
                    axis.ticks = element_line(size = linesize),
                    axis.text = element_text( size = fontsize),
                    axis.title = element_text(size = fontsize,family = "sans"),
                    axis.line = element_line(size = linesize))
    }else {
       p = ggplot(data) +
            geom_point( data = data,aes_string(x = "PC1", y = "PC2", color = Group_color),
                        size = pt_size  ) + 
            # xlab(paste0("PC1 (",round(deviation[1],1),"% )")) +
            # ylab(paste0("PC2 (",round(deviation[2],1),"% )")) +
            scale_color_manual( values = alpha(color ,alpha = alpha) )+ 
            theme_cowplot()+
            theme(  text = element_text( size = fontsize),
                    axis.ticks = element_line(size = linesize),
                    axis.text = element_text( size = fontsize),
                    axis.title = element_text(size = fontsize,family = "sans"),
                    axis.line = element_line(size = linesize))
    }
    if(!rlang::is_empty(deviation)){
        p = p + xlab(paste0("PC1 (",round(deviation[1],1),"% )")) +
                ylab(paste0("PC2 (",round(deviation[2],1),"% )"))
    }
    return(p)
}

require(ggplot2)

###***--- For single cell UMAP / tSNE left-bottom small axis plotting ---***###
# NoAxes() + patchwork::inset_element(axis, left = 0, bottom = 0, right = 0.15, top = 0.17, align_to = 'full',on_top = FALSE) 
axis = ggplot( data = data.frame(0:1) ) + 
  geom_point(aes(x = 0:1,y = 0:1),color = alpha("black",0)) + 
  cowplot::theme_cowplot() + xlab("UMAP1") + ylab("UMAP2") + 
  scale_x_continuous( labels = NULL,limits = c(0,1) ) + 
  scale_y_continuous( labels = NULL,limits = c(0,1) ) + 
  theme( axis.title = element_text( size = 10, hjust = 0.5, vjust = 0),
         axis.ticks = element_blank(),
         axis.line = element_line(linewidth = 0.8,arrow = arrow(type='closed', length = unit(5,'pt') ) )) 

###***--- For single cell boxplotting ---***###
library(reshape2)
BoxPlot = function( seurat.data,
                    features = character(),
                    ncol = 3,
                    cell.type = character() ){
  mat = t(as.matrix(ho.hsc[["SCT.regressed"]] @data[ gene.list,]))
  cell.cluster = data.frame(Cluster = paste0("t",cell.type,as.integer(Idents(ho.hsc))),
                            row.names = names(Idents(ho.hsc)))
  mat = cbind(cell.cluster,mat) %>% rownames_to_column("Cells")
  mat = melt( mat,id.vars = c("Cluster","Cells"),variable.name = "Genes")
  p = ggplot( mat ) +
    geom_boxplot( aes( x = Cluster, y = value, fill =  Cluster) ) + 
    #scale_fill_manual(values =brewer.pal( length(levels(Idents(ho.hsc))),"Paired" ),) +
    facet_wrap( .~Genes, ncol = ncol) +
    theme_cowplot() + ylab( "Normalized Expression" )+ xlab(NULL)+
    theme( text = element_text( family = "sans",size = 20 ),
           axis.text =  element_text( family = "sans",size = 16 ),
           strip.placement = "inside",
           strip.background = element_rect( color = NA,fill = NA ),
           legend.position = "none",
           panel.border = element_rect( colour = "black",size = 1 ))
  p
}

###***--- For GO plotting with text and filtration ---***###
GO_Plot = function( GO_dataset, padj_threshold=0.05, show_num = 15, 
                    fill_color = "#e63f46",wid = 10, hgt = 8,
                    fontsize = 2,fontcolor = "black",
                    barwidth = 0.8,
                    GO_term = character(),
                    keywords = character(),
                    discard = character()){
  tmp = GO_dataset@result
  tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
  tmp$logpval = (-log10( tmp$pvalue ))
  tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
  
  if(!is_empty( GO_term )){
    tmp = tmp %>% filter( ONTOLOGY %in% GO_term )
  }
  if(!is_empty( keywords )){
    seek = numeric()
    for( i in keywords){
      j = grep(i,tmp$Description,value = T)
      seek = c(seek,j)
    }
    if(length(seek)==0){return("No aim GO terms was found!")}
  }
  if(show_num==0){
    all_term = unique(seek)
  }else {
    all_term = unique(c(tmp$Description[1:min(length(tmp$Description),show_num)],seek))
  }
  if(length(discard)!=0){
    remove = numeric()
    for( i in discard){
      j = grep(i,tmp$Description,value = T)
      remove = c(remove,j)
    }
    if(length(remove)==0){return("No discarded GO terms was found!")}
    tmp = tmp[tmp$Description %in% setdiff(tmp$Description,remove),]
  }
  
  tmp = tmp[tmp$Description %in% all_term,] %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
  tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
  go = ggplot()+
    geom_bar( data=tmp,
              aes_string(x = "Description" ,y = "logpval"),
              fill = fill_color,stat = "identity",width = barwidth) + ylab( "-log(p value)" ) +
    geom_text( data=tmp,aes(x = Description ,y = 0.1, label = Description),
               hjust = 0,size = fontsize,color = fontcolor,family = "sans" )+
    xlab(NULL)+scale_y_continuous(expand = c(0,0)) + 
    coord_flip( )+
    # facet_grid( drop = T,rows = "ONTOLOGY",scales = "free",space = "free_y")+
    theme_cowplot()+
    theme( text = element_text( size = 30),
           axis.title = element_text(size = 30),
           axis.text.y = element_blank(),
           axis.text.x= element_text(size = 30,family = "sans"),
           axis.line = element_line(size = 1.5),
           axis.ticks = element_line(size = 1.5),
           axis.ticks.y = element_blank()
           
    )
  return(go)
}




# load( "../exvivo/parameter.rds/cellcycle.genes.RData" )

# for (i in c("Ho","Tx")) {
#   load( file = list.files( "data",pattern = i,full.names = T ) )
#   seurat.tmp = CreateSeuratObject( counts = cell.ho, 
#                                    project = i,
#                                    meta.data = info.ho,
#                                    min.features = 1000 )
#   seurat.tmp[["percent.MT"]] = PercentageFeatureSet( seurat.tmp,pattern = "^mt-" )
#   assign( paste0(i,".seurat"),seurat.tmp )
#   rm(seurat.tmp)
# }

# Ho.TX.seurat = merge( Ho.seurat,Tx.seurat,add.cell.id = c("Ho", "Tx") )
# table(Ho.TX.seurat$Ho_Tx)   
# save(Ho.seurat,Tx.seurat,Ho.TX.seurat,file = "data/seurat.RData")
# Ho   Tx 
# 1270 1058
# new.seurat = lapply(list(Ho.seurat,Tx.seurat), function(x) {
#   seurat.tmp = x
#   seurat.tmp = SCTransform(seurat.tmp,new.assay.name = "SCT",
#                            vars.to.regress = "percent.MT",do.scale = F,do.center = F,
#                            return.only.var.genes = T,
#                            variable.features.n = nrow( seurat.tmp ))
#   seurat.tmp = CellCycleScoring( seurat.tmp,
#                                  s.features = cc.s$MGI.symbol,
#                                  g2m.features = cc.g2m$MGI.symbol,
#                                  set.ident = T  )
#   seurat.tmp = FindVariableFeatures( seurat.tmp,
#                                      nfeatures = 2000,
#                                      selection.method = "vst",
#                                      verbose = FALSE )
#   seurat.tmp = SCTransform( seurat.tmp,
#                             new.assay.name = "SCT.regressed",
#                              assay = "RNA",
#                              do.scale = T,
#                              vars.to.regress = c("percent.MT","S.Score","G2M.Score"),
#                              return.only.var.genes = T,
#                              variable.features.n = nrow(seurat.tmp)) # default parameters
#   return(seurat.tmp)
# })
# Ho.seurat = new.seurat[[1]]
# Tx.seurat = new.seurat[[2]]
# save(Ho.seurat,Tx.seurat,Ho.TX.seurat,new.seurat,file = "data/seurat.RData",compress = T)
# save(Ho.seurat,Tx.seurat,Ho.TX.seurat,file = "data/seurat.small.RData",compress = T)
