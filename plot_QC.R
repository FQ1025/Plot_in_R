### plot functions ###

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

data_summary <- function(x) { 
  # m <- mean(x)
  m = median(x)
  ymin <- m-sd(x) 
  ymax <- m+sd(x) 
  return(c(y=m,ymin=ymin,ymax=ymax)) 
}
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

