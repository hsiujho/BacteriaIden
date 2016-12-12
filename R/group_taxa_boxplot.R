#
#
#
#
#
#
#

group_taxa_boxplot=function(phylo,group.varn=NA,ranklv="Phylum",topn=5,facet.ncol=3){

  if(all(sample_variables(phylo)!=group.varn)|is.na(group.varn)){
    sample_data(phylo)$dummy_var=""
    group.varn="dummy_var"
  }

  a2=phylo %>%
    tax_glom(ranklv,NArm=T) %>%
    transform_sample_counts(function(x){x/sum(x)*100})

  a3=taxa_sums(a2) %>% order(decreasing=T) %>% "["(-(1:topn)) %>%
    merge_taxa(a2,.,2) %>% psmelt() %>%
    mutate_(rank_var=ranklv,x_var=group.varn) %>%
    mutate_each(funs(as.character),rank_var)
  a3$rank_var%<>%ifelse(is.na(.),"Others",.)

  a3$rank_var%<>%factor(levels=group_by(a3,rank_var)%>%
                          summarise(RA=sum(Abundance))%>%arrange(desc(RA))%$%rank_var)

  ggplot(a3,aes(x=x_var,y=Abundance))+geom_boxplot()+ylab("Relative abundance")+
    xlab(ifelse(group.varn=="dummy_var","",group.varn))+
    facet_wrap(~rank_var,scales = "free_y",ncol=facet.ncol)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,v=0.5))
}
