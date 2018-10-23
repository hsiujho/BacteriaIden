#
#
#
#
#
#

phylo_alpha_div_boxplot=function(phylo,ranklv,group_var,test_for="mean",add_points=T){
  if(missing(ranklv)) ranklv="OTU"
  if(any(rank_names(phylo)==ranklv)) {
    a0=my_tax_glom(phylo,ranklv,NArm=F)
  } else {
    a0=phylo
  }
  index=c("Observed","Chao1","Shannon","Simpson","InvSimpson")
# sample_names 由數字開頭會在名稱前端增加"X"字元, 特殊符號會被轉換成"."
  a1=estimate_richness(a0,measures=index) %>%
    rownames_to_column("Subjects") %>%
    mutate(SampleID=gsub("\\.","-",Subjects) %>>%
             (ifelse(grepl("^\\d.*",sample_names(a0)),gsub("^X","",.),.))) %>%
    melt("SampleID",index,"Index")

  ylabtext=ifelse(any(rank_names(phylo)==ranklv),sprintf("Expression at %s level",ranklv),"Expression at OTU level")
  if(missing(group_var)) group_var=""
  if(any(sample_variables(phylo)==group_var)){
    b1=my_get_variable(phylo,group_var)
    a1=merge(a1,b1,by.x="SampleID",by.y=0)
    l_item=length(unique(a1[[group_var]]))
    test.formula=as.formula(sprintf("value~%s",group_var))
    if(l_item==2) {
      t1=group_by(a1,Index) %>>%
        do(k1=wilcox.test(test.formula,data=.,exact=F)
           ,k2=t.test(test.formula,data=.)) %>>%
        summarise(Index,median.test.pvalue=k1$p.value
                  ,mean.test.pvalue=k2$p.value)
    } else if (l_item>2) {
      t1=group_by(a1,Index) %>>%
        do(k1=kruskal.test(test.formula,data=.)
           ,k2=anova(lm(test.formula,data=.))) %>>%
        summarise(Index,median.test.pvalue=k1$p.value
                  ,mean.pvalue=k2$`Pr(>F)`[1])
    }
    t2=group_by_at(a1,c("Index",group_var)) %>>%
      summarise(N=length(value),Mean=mean(value),Median=median(value))
    p1=ggplot(a1,aes_string(x=group_var,y="value")) +geom_boxplot()+facet_wrap(~Index, scales = "free_y",nrow=1)+ylab(ylabtext)
    p11=ggplot_build(p1)
    p12=lapply(p11$layout$panel_params,function(x){
      return(c(x$x.range,x$y.range))
    }) %>>%
      (do.call(rbind,.)) %>>%
      data.frame() %>>%
      setNames(c("x_min","x_max","y_min","y_max")) %>>%
      mutate(Index=p11$layout$layout$Index) %>>%
      merge(t1,by="Index") %>>%
      mutate(pvalue=ifelse(rep(test_for=="mean",5),mean.test.pvalue,median.test.pvalue) %>>%
               (ifelse(.<0.001,"P < 0.001",ifelse(.>0.999,"P > 0.999",sprintf("P = %.3f",.)))))
    p2=p1+geom_text(
      data = p12,
      mapping = aes(x = x_min, y = y_max, label = pvalue)
      ,vjust=-0.05
      ,hjust=-0.05
    )+
      theme(legend.background=element_blank(),
            plot.background = element_blank())
    if(add_points){
      p2=p2+
        geom_jitter(width = 0.25,colour="black", alpha=0.3)
    }
    a2=sprintf("SampleID+%s~Index",group_var) %>% formula() %>% dcast(a1,.)
    return(list(div_tbl=a2,summar_tbl=t2,test_tbl=t1,boxplot=p2))
  } else {
    p1=ggplot(a1,aes_string(x=factor(""),y="value")) +geom_boxplot()+facet_wrap(~Index, scales = "free_y",nrow=1)+xlab("")+ylab(ylabtext)
    if(add_points){
      p1=p1+
        geom_jitter(width = 0.25,colour="black", alpha=0.3)
    }
    p2=p1+
                  theme(legend.background=element_blank()
                        ,plot.background = element_blank()
                        ,axis.ticks.x.bottom = element_blank()
                        ,axis.title.x = element_blank())
    t2=group_by_at(a1,"Index") %>>%
      summarise(N=length(value),Mean=mean(value),Median=median(value))
    a2=dcast(a1,SampleID~Index)
    return(list(div_tbl=a2,summar_tbl=t2,boxplot=p2))
  }
}
