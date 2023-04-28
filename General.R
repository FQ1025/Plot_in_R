###***--- For geneal plot themes ---***###
Themes = function(pt_size = 8,
                  alpha = 1,
                  fontsize = 28,
                  linesize = 2){
    require(cowplot)
    require(ggplot2)
    a = theme_cowplot()+
            theme(  text = element_text( size = fontsize),
                    axis.ticks = element_line(size = linesize),
                    axis.text = element_text( size = fontsize),
                    axis.title = element_text(size = fontsize,family = "sans"),
                    axis.line = element_line(size = linesize))
    return(a)
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



###***--- For statistics test plotting ---***###
data_summary <- function(x) { 
       m <- mean(x) 
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