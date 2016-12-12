#
#a0: object of phyloseq
#ranklv: one of rank_names
#group.varn: name of group variable
#topn: number of leading taxanomy
#
#
#

group_taxa_mean=function(a0,ranklv,group.varn,topn=5){
  a2=a0 %>%
    tax_glom(ranklv,NArm=T) %>%
    transform_sample_counts(function(x){x/sum(x)*100})
  sample_data(a2)$SubjID=sample_names(a2)

  a3=taxa_sums(a2) %>% order(decreasing=T) %>% "["(1:topn) %>%
    "["(taxa_names(a2),.) %>% "%in%"(taxa_names(a2),.) %>%
    prune_taxa(a2) %>% psmelt() %>%
    acast(formula(sprintf("%s~%s",ranklv,"SubjID")),value.var ="Abundance") %>%
    rbind(Others=100-colSums(.)) %>%
    melt(varnames=c(ranklv,"SubjID"),value.name="Abundance") %>%
    mutate_each_(funs(as.character),"SubjID") %>%
    left_join(
      sample_data(a2) %>% data.frame() %>% select_("SubjID",group.varn)
      ,by="SubjID") %>%
    group_by_(group.varn,ranklv) %>%
    dplyr::summarise(RA=mean(Abundance)) %>% ungroup() %>%
    mutate_each_(funs(as.character),c(group.varn,ranklv))
  return(a3)
}
