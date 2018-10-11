#
#
#
#
#
#

phylo_dendro_bar_2=function(phylo,ranklv="Phylum",topn=6
                            ,dendro_height=30
                            ,dendro_ypos=100
                            ,bar_width=0.4 #0~0.5
                            ,group_var
                            ,group_height=-5
                            ,group_grid=-1
                            ,group_ypos=-2
                            ,hc_method="ward.D"
                            ,dist_method="bray"
                            ,cell_border_colour="black"
                            ,cell_border_size=0.1
                            ,using_path=T
                            ,annotation_colors
){
  #  ranklv="Class"
  if(!"SampleID"%in%sample_variables(phylo)){
    phylo@sam_data$SampleID=sample_names(phylo)
  }

  b1=phylo_dendrogram(phylo, ranklv=ranklv,method=hc_method,dist_method=dist_method,using_path=using_path)

  c1=phylo_bar_plot(phylo,fill = ranklv,x = "SampleID",topn = topn)

  c3=c1$barplot$data %>>%
    (split(.,.$SampleID)) %>>%
    lapply(function(x){
      arrange_at(x,ranklv) %>>%
        mutate(cumsum=cumsum(Abundance)
               ,ymin=c(0,head(cumsum,-1))
               ,ymax=ymin+Abundance)
    }) %>>%
    bind_rows() %>>%
    merge(b1$hcdata$labels,by.x="SampleID",by.y="label") %>>%
    mutate(xmin=x-bar_width,xmax=x+bar_width)

  # dendro_height=30
  # dendro_ypos=100

  dendro_origin_y_range=b1$hcdata$segments[,c("y","yend")] %>>%
    sapply(range) %>>%
    unlist() %>>%
    range()
  shrink_y=function(val){
    up=val-dendro_origin_y_range[1]
    down=diff(dendro_origin_y_range)
    up/down*dendro_height+dendro_ypos
  }
  p0=ggplot(c3,aes_string(xmin="xmin",xmax="xmax",ymin="ymin",ymax="ymax",fill=ranklv))+
    geom_rect(colour=cell_border_colour,size=cell_border_size)
  if(using_path){
    p1=p0+
      geom_path(aes(x=x,y=shrink_y(y),group=k)
                   ,data=b1$fig$data
                   ,inherit.aes = F)

  } else {
    p1=p0+
      geom_segment(aes(x=x,y=shrink_y(y)
                       ,xend=xend,yend=shrink_y(yend))
                   ,data=b1$hcdata$segments
                   ,inherit.aes = F)
  }
  p1=p1+
    theme_minimal()+
    scale_x_continuous(breaks = b1$hcdata$labels$x
                       ,labels = b1$hcdata$labels$label,expand = c(0.01,0))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5)
          ,panel.grid.major.y = element_line(colour="gray", size=0.5)
          ,panel.grid.major.x = element_blank()
          ,panel.grid.minor = element_blank())+
    labs(y="Relative abundnace (%)",x="Sample")

  # p1+theme(axis.text.x = element_blank())
  # p1 + coord_fixed(ratio = 0.2)

  # 增加分組資訊
  # group_var=c("group","seq")
  # group_height=5
  # group_grid=1
  # group_ypos=-5

  if(missing(group_var)){
    p2=p1+scale_y_continuous(breaks = seq(0,100,20)
                             ,labels=seq(0,100,20)
                             ,expand = c(0,0))
    return(p2)
  } else {
    p2=p1
    group_plot_list=list()
    for(k in 1:length(group_var)){
      df=my_get_variable(phylo,group_var[k]) %>>%
        merge(b1$hcdata$labels,by.x=0,by.y="label") %>>%
        mutate(x1=x-bar_width,x2=x+bar_width
               ,y1=group_ypos+group_height*k+group_grid*(k-1)
               ,y2=group_ypos+group_height*(k-1)+group_grid*(k-1))

      group_var_item=unique(df[[group_var[k]]])
      df_col=data.frame(item=sort(group_var_item),col_id=scales::hue_pal()(length(group_var_item)),stringsAsFactors = F)
      if(!missing(annotation_colors)){
        acol=annotation_colors[[group_var[k]]]
        if(!is.null(acol)){
          df_col2=data.frame(item=names(acol)
                             ,col_id2=acol,stringsAsFactors = F)
          df_col=merge(df_col,df_col2,by="item",all.x=T) %>>%
            mutate(col_id=ifelse(is.na(col_id2),col_id,col_id2)) %>>%
            select(-col_id2)
        }
      }
      df_p1=ggplot(df,mapping=aes_string(xmin="x1",xmax="x2",ymin="y1",ymax="y2",fill=group_var[k]))+
        geom_rect(colour=cell_border_colour,size=cell_border_size)+
        scale_fill_manual(group_var[k],values = setNames(df_col$col_id,df_col$item))
      group_plot_list[[group_var[k]]]=df_p1
      for(l in 1:NROW(df_col)){
        p2=p2+geom_rect(data=filter_at(df,group_var[k],all_vars(.==df_col$item[l]))
                        ,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2)
                        ,colour=cell_border_colour
                        ,size=cell_border_size
                        ,inherit.aes = F,fill=df_col$col_id[l])
      }
    }
    group_y_posi=(group_ypos+(0:length(group_var))*group_height) %>>%
      ((head(.,-1)+tail(.,-1))/2)+(0:(length(group_var)-1))*group_grid
    p2=p2+scale_y_continuous(breaks = c(group_y_posi,seq(0,100,20))
                             ,labels=c(group_var,seq(0,100,20))
                             ,expand = c(0,0))+
      theme(legend.position = "none")

    #取出ggplot的legend
    group_legend=lapply(group_plot_list, g_legend)
    bar_legend=g_legend(p0)
    dev.off()
    empty_rect=rectGrob(gp = gpar(lwd = 2, col = NA, fill = NA))

    return(list(main=p2,group_legend=group_legend
                ,bar_legend=bar_legend,empty_rect=empty_rect))
  }
}


function(){
  require(BacteriaIden)
  rm(list=ls())
  phylo=subset_samples(MS16036_gg_phylo
                       ,group%in%c("N","NE","N_FNE")&ST=="After") %>>%
    (prune_taxa(taxa_sums(.)>0,.))
  p1=phylo_dendro_bar_2(phylo,ranklv="Phylum",cell_border_colour = "gray")
  p1=phylo_dendro_bar_2(phylo,ranklv="Phylum",cell_border_colour = "gray",using_path = F)
  p1=phylo_dendro_bar_2(phylo,ranklv="Class",topn=10)
  p1=phylo_dendro_bar_2(phylo,ranklv="Genus",topn=15)

  p2=phylo_dendro_bar_2(phylo,ranklv="Class",topn=10,group_var = "group")
  grid.arrange(p2$main+theme(axis.text.x = element_blank())
               ,arrangeGrob(p2$empty_rect
                            ,p2$bar_legend$grobs[[1]]
                            ,p2$group_legend$group$grobs[[1]]
                            ,p2$empty_rect,ncol=1,heights = c(1,3,2,1))
               ,ncol=2,widths=c(2,1))


  p2=phylo_dendro_bar_2(phylo,ranklv="Genus",topn=15
                        ,group_var = c("group","ST"))

  grid.arrange(p2$main+theme(axis.text.x = element_blank())
               ,arrangeGrob(p2$empty_rect
                            ,p2$bar_legend$grobs[[1]]
                            ,p2$empty_rect,ncol=1,heights = c(0.5,3,0.5))
               ,arrangeGrob(p2$empty_rect
                            ,p2$group_legend$group$grobs[[1]]
                            ,p2$group_legend$ST$grobs[[1]]
                            ,p2$empty_rect,ncol=1,heights = c(0.5,2,1,0.5))
               ,ncol=3,widths=c(2,1,.5))

  # 移除 Others
  p2$main$data %<>% filter(Genus!="Others")
  grid.arrange(p2$main+theme(axis.text.x = element_blank())
               ,arrangeGrob(p2$empty_rect
                            ,p2$bar_legend$grobs[[1]]
                            ,p2$empty_rect,ncol=1,heights = c(0.5,3,0.5))
               ,arrangeGrob(p2$empty_rect
                            ,p2$group_legend$group$grobs[[1]]
                            ,p2$group_legend$ST$grobs[[1]]
                            ,p2$empty_rect,ncol=1,heights = c(0.5,2,1,0.5))
               ,ncol=3,widths=c(2,1,.5))

  #
  p3=phylo_dendro_bar_2(phylo,ranklv="Genus",topn=15
                        ,group_var = c("group","ST")
                        ,dendro_height=30
                        ,dendro_ypos=115
                        ,group_height=5
                        ,group_grid=1
                        ,group_ypos=102
  )
  p3$main
  # 自訂variable 的顏色
  p4=phylo_dendro_bar_2(phylo,ranklv="Class",topn=10,group_var = c("group","ST")
                        ,annotation_colors = list(group=c(N="#e41a1c",NE="#377eb8",N_FNE="#4daf4a")))
  p4$main
}
