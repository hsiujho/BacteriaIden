#
#
#
#
#
#

phylo_bar_plot=function(
  phylo
  ,fill="Phylum" #生物層次
  ,x #分組變數, 至少一個變數名稱, 若要進一步套用facet_*, 則類別變數由大至小給定, 再將此函數的結果自行套用facet_*函數.
  ,topn=6 #只標註表現量最多的taxonomy, 其餘用Others彙總.
  ,NArm=F
){
  if(!fill%in%rank_names(phylo)) stop("The parameter of fill must be one of rank_names(phylo).")
  p1=my_tax_glom(phylo,fill,NArm=NArm) %>>% {
    if(NArm) {
      .
    } else {
      add_label_to_unknown_taxon(.)
    }
  } %>>%
    transform_sample_counts(function(x)x/sum(x)*100) %>>%
    psmelt()
  fillv=rank_names(phylo) %>>% ("["(.,1:which(.==fill)))
  p2=group_by_at(p1,.vars=fillv) %>>%
    summarise(Abundance=sum(Abundance)) %>>%
    arrange(desc(Abundance)) %>>%
    top_n(topn,Abundance) %>>%
    "[["(fill) %>>%
    as.character()

#  p3=subset_taxa(p1,lazyeval::interp("x %in% p2", x = fill)) %>>% psmelt() %>% group_by_(.dots=c(x,fill)) %>>% summarise(Abundance=mean(Abundance)) %>>% ungroup() %>>% mutate_each_(funs(as.character),c(x,fill))

  p3=p1[p1[[fill]]%in%p2,] %>%
    group_by_at(.vars=c(x,fillv)) %>>%
    summarise(Abundance=mean(Abundance)) %>>%
    ungroup() %>>%
    mutate_at(.vars=c(x,fillv),funs(as.character))

#  p3=dplyr::filter(p1,Phylum%in%p2)
#lazyeval::interp("~f(x,y)", x = fill,y=p2,f=as.name("%in%"))

  #為每一個分群增加others的類別
  #library(lazyeval)

  p4=group_by_at(p3,.vars=x) %>>%
    summarise(Abundance=100-sum(Abundance)) %>>%
    ungroup() #%>>%
    #mutate_(.dots=lazyeval::interp("x", x = "Others")%>>%setNames(fill))
  if(all(p4$Abundance<=0)){
    p5=p3
    #設定fill菌叢的順序
  #  p6=group_by_(p3,.dots=fill) %>>% summarise(Abundance=sum(Abundance)) %>>% arrange(desc(Abundance)) %>>% "[["(fill) %>>% c("Others")
    p6=unique(p2)
  } else {
    p4[[fill]]="Others"
    p5=bind_rows(p3,p4)
    p6=c(unique(p2),"Others")
  }
  p5[[fill]]%<>%factor(levels=p6)

  p = ggplot(arrange_at(p5,.vars=fill), aes_string(x = tail(x,1), y = "Abundance", fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", color = "black")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))+
    theme(legend.background=element_blank(),
          plot.background = element_blank())

  blank_theme =theme(axis.title.y=element_blank()
                     ,axis.ticks.y=element_blank()
                     ,axis.text.y=element_blank()
                     ,legend.background=element_blank()
                     ,plot.background = element_blank())
  pie=ggplot(arrange_(p5,fill),aes_string(x="factor(1)",y="Abundance",fill=fill))+
    geom_bar(stat = "identity", color = "black",width = 1)+
    coord_polar(theta="y")+ylab("Relative abundance")+
    blank_theme+facet_wrap(formula(sprintf("~%s",tail(x,1))),ncol=5)

  return(list(barplot=p,piechart=pie))
}

#data(GlobalPatterns)
#phylo_bar_plot(phylo=GlobalPatterns,fill="Phylum",x="SampleType")
#phylo_bar_plot(phylo=GlobalPatterns,fill="Genus",x="SampleType",topn=15,NArm=F)
#
