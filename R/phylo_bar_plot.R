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
){

  p1=tax_glom(phylo,fill) %>>% transform_sample_counts(function(x)x/sum(x)*100) %>>% psmelt()
  p2=group_by_(p1,fill) %>>% summarise(Abundance=sum(Abundance)) %>>% arrange(desc(Abundance)) %>>% top_n(topn,Abundance) %>>% "[["(fill) %>>% as.character()
  #print(p1)

#  p2=taxa_sums(p1) %>>% (names(.)[tail(order(.),topn)])
#  p2=taxa_sums(p1)
#  p2=names(p2)[tail(order(p2),topn)]
#  p31=taxa_names(p1) %in% p2

#  print(head(p31))
#  p3=subset_taxa(p1,lazyeval::interp("x %in% p2", x = fill)) %>>% psmelt() %>% group_by_(.dots=c(x,fill)) %>>% summarise(Abundance=mean(Abundance)) %>>% ungroup() %>>% mutate_each_(funs(as.character),c(x,fill))
#  p21=p1[[fill]]%in%p2
  p3=p1[p1[[fill]]%in%p2,] %>% group_by_(.dots=c(x,fill)) %>>% summarise(Abundance=mean(Abundance)) %>>% ungroup() %>>% mutate_each_(funs(as.character),c(x,fill))

#  p3=dplyr::filter(p1,Phylum%in%p2)
#lazyeval::interp("~f(x,y)", x = fill,y=p2,f=as.name("%in%"))

  #為每一個分群增加others的類別
  #library(lazyeval)

  p4=group_by_(p3,.dots=x) %>>% summarise(Abundance=100-sum(Abundance)) %>>% ungroup() %>>% mutate_(.dots=lazyeval::interp("x", x = "Others")%>>%setNames(fill))

  p5=bind_rows(p3,p4)
  #設定fill菌叢的順序
#  p6=group_by_(p3,.dots=fill) %>>% summarise(Abundance=sum(Abundance)) %>>% arrange(desc(Abundance)) %>>% "[["(fill) %>>% c("Others")
  p6=c(p2,"Others")
  p5[[fill]]%<>%factor(levels=p6)

  p = ggplot(arrange_(p5,fill), aes_string(x = tail(x,1), y = "Abundance", fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", color = "black")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))

  blank_theme =theme(axis.title.y=element_blank()
                     ,axis.ticks.y=element_blank()
                     ,axis.text.y=element_blank())
  pie=ggplot(arrange_(p5,fill),aes_string(x="factor(1)",y="Abundance",fill=fill))+
    geom_bar(stat = "identity", color = "black",width = 1)+
    coord_polar(theta="y")+ylab("Relative abundance")+
    blank_theme+facet_wrap(formula(sprintf("~%s",tail(x,1))),ncol=5)

  return(list(barplot=p,piechart=pie))
}

#data(GlobalPatterns)
#phylo_bar_plot(phylo=GlobalPatterns,fill="Phylum",x="SampleType")
#
#
