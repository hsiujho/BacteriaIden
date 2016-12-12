#
#
#
#
#
#

phylo_dendrogram=function(phylo,group_var,toRA=T,ranklv="OTU",method="complete",type="rectangle",autolegend=T){
  if(any(rank_names(phylo)==ranklv)){
    phylo%<>%tax_glom(ranklv,NArm=F)%>>%(prune_taxa(taxa_sums(.)>0,.))
  }
  if(toRA){
    phylo%<>%transform_sample_counts(function(x)x/sum(x)*100)
  }

  hc <- otu_table(phylo)@.Data %>% t() %>% dist() %>% hclust(method=method)
  hcdata <- dendro_data(hc, type=type)
  p1=ggplot() +
    geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_dendro()+theme(axis.title=element_blank())
  p1=p1+
    scale_x_continuous(breaks=1:nrow(label(hcdata)),labels=label(hcdata)$label)

  if(!missing(group_var)){
    if(any(group_var%in%sample_variables(phylo))){
      if(is.factor(get_variable(phylo,group_var))){
        a0=get_variable(phylo,group_var) %>%
          data.frame(stringsAsFactors=F) %>%
          setNames(group_var)
        rownames(a0)=sample_names(phylo)
        a1=rownames_to_column(a0,"label")
        hcdata$label%<>%left_join(a1,by="label")

        p2=p1+
          theme(axis.text.x = element_text(angle = 90, hjust = 0
                                           , vjust=0.5,color=hcdata$label[[group_var]]))

        if(autolegend){
          a2=p2$theme$axis.text.x$colour %>>%
            (data.frame(levels(.),1:nlevels(.)
                        ,x=ggplot_build(p2)$panel$ranges[[1]]$x.range[2]
                        ,y=ggplot_build(p2)$panel$ranges[[1]]$y.range[2]) %>%
               setNames(c(group_var,"i","x","y")))
          #ggplot_build(p2)$panel$ranges[[1]]$y.range
          #ggplot_build(p2)$panel$ranges[[1]]$x.range
          plot(hc)
          a3=wordlayout(x=a2$x, y=a2$y, words=a2$group
                        ,xlim=ggplot_build(p2)$panel$ranges[[1]]$x.range
                        ,ylim=ggplot_build(p2)$panel$ranges[[1]]$y.range) %>%
            data.frame() %>% rownames_to_column(group_var)
          dev.off()
          a2$yj=a2$y-a2$i*a3$ht-(a2$i-1)*a3$ht/3
          a2$xj=a2$x-max(a3$width)

          p3=p2+annotate("text", x = a2$xj, y = a2$yj, label = a2[[group_var]],color=1:3,hjust=0)
        } else {
          p3=p2
        }
        return(list(hc=hc,hcdata=hcdata,fig=p3))
      }
    }
  }
  p3=p1+theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))
  return(list(hc=hc,hcdata=hcdata,fig=p3))
}
