#
#
#
#
#
#

phylo_div_boxplot=function(phylo,group.varn=""){

  if(all(sample_variables(phylo)!=group.varn)|is.na(group.varn)){
    sample_data(phylo)$dummy_var=""
    group.varn="dummy_var"
  } else {
    sample_data(phylo)$dummy_var=get_variable(phylo,group.varn)
  }
  sample_data(phylo)$SampleID=sample_names(phylo)

  b1=lapply(c("Phylum","Genus"),function(x){
    tax_glom(phylo,x,NArm=F) %>%
      estimate_richness() %>% add_rownames("Subjects") %>%
      melt("Subjects",colnames(.)[-1],"Index") %>% mutate(Rank=x)
  }) %>% bind_rows() %>% mutate_each(funs(as.character),Index) %>% filter(!Index%in%c("se.chao1","ACE","se.ACE")) %>%
    mutate(Rank=factor(Rank,levels=c("Phylum","Genus"))
           ,Index=factor(Index,levels=c("Shannon","Simpson","InvSimpson","Observed","Chao1","Fisher"))
           ,SampleID=gsub("\\.","-",Subjects)) %>%
    left_join(get_variable(phylo,c("SampleID","dummy_var")),by="SampleID")

  p1=filter(b1,Rank=="Phylum") %>%
    ggplot(aes(x=dummy_var,y=value)) +geom_boxplot()+facet_grid(Index~Rank,scales = "free",switch="y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,v=0.5))+ylab("")+xlab(ifelse(group.varn=="dummy_var","",group.varn))
  p2=filter(b1,Rank=="Genus") %>%
    ggplot(aes(x=dummy_var,y=value)) +geom_boxplot()+facet_grid(Index~Rank,scales = "free")+ theme(axis.text.x = element_text(angle = 90, hjust = 1,v=0.5))+ylab("")+xlab(ifelse(group.varn=="dummy_var","",group.varn))

  grid.arrange(p1,p2,ncol=2)
}
