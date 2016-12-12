#
# 使用LEfSe的輸出, 並指定欄位名稱c("Statistic","meandiff","log2meandiff"),
# 依生物層畫出有顯著的taxanomy
# 輸出為list, 每個元素為ggplot2物件
#
#
#

LEfSe_bar=function(LEfSe_outp,colv){
  v0=LEfSe_outp
  rank_lv=c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  v1=left_join(v0$LEfSe,v0$meandiff,by="Taxonomy") %>%
    dplyr::filter(V5!="-"&Group!="") %>%
    mutate(Rank_order=ordered(Rank,rank_lv))
  v3=rank_lv[-1] %>% "["(.%in%unique(v1$Rank_order))
  if(length(v3)==0) stop("There have not significant data.")
  v4=lapply(v3,function(i)
  {

    v2=dplyr::filter(v1,Rank==i) %>%
      group_by(Rank,Group) %>%
      arrange_(sprintf("desc(abs(%s))",colv)) %>%
      mutate(taxo_label=strsplit(Taxonomy,"\\.") %>%
               sapply(function(x){rev(x)[1]})) %>%
      select_("Taxonomy","taxo_label","Rank","Group",colv) %>%
      mutate_(abs_stat=sprintf("abs(%s)",colv))%>% top_n(10,abs_stat)#%>%arrange_(colv)
    v2$tax_lab_ord=  with(v2,factor(x=taxo_label,levels=taxo_label))

    eval(parse(text=sprintf("p4=ggplot(v2,aes(x=tax_lab_ord,y=%s))",colv)))

    p4=p4 + geom_bar(stat = "identity",position="dodge") +theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5))+xlab("")
    #  p4+ facet_wrap(~Rank_order+Group,scales = "free_x",ncol=2)
    #  p4=p4+ facet_grid(.~Rank+Group,scales = "free_x", space ="free")
    p4=p4+ facet_grid(.~Group,scales = "free_x", space ="free")+ggtitle(i)
    return(p4)
  })
  #ml <- marrangeGrob(v4, nrow=3, ncol=2,top="")
  #ml
  return(v4)
}
