#
#
#
#
#
#

phylo_alpha_div_boxplot=function(phylo,ranklv,group_var){
  if(missing(ranklv)) ranklv="OTU"
  if(any(rank_names(phylo)==ranklv)) {
    a0=tax_glom(phylo,ranklv,NArm=F)
  } else {
    a0=phylo
  }
  index=c("Observed","Chao1","Shannon","Simpson","InvSimpson")
  a1=a0 %>% estimate_richness(measures=index) %>% rownames_to_column("Subjects") %>% mutate(SampleID=gsub("\\.","-",Subjects)) %>% melt("SampleID",index,"Index")

  ylabtext=ifelse(any(rank_names(phylo)==ranklv),sprintf("Expression at %s level",ranklv),"Expression at OTU level")
  if(missing(group_var)) group_var=""
  if(any(sample_variables(phylo)==group_var)){
    b1=data.frame(sample_names(phylo),get_variable(phylo,group_var),stringsAsFactors = F) %>% setNames(c("SampleID",group_var))
    a1%<>%left_join(b1,by="SampleID")
    p1=ggplot(a1,aes_string(x=group_var,y="value")) +geom_boxplot()+facet_wrap(~Index, scales = "free_y",nrow=1)+ylab(ylabtext)
    a2=sprintf("SampleID+%s~Index",group_var) %>% formula() %>% dcast(a1,.)
  } else {
    p1=ggplot(a1,aes_string(x=factor(""),y="value")) +geom_boxplot()+facet_wrap(~Index, scales = "free_y",nrow=1)+xlab("")+ylab(ylabtext)
    a2=dcast(a1,SampleID~Index)
  }
  return(list(table=a2,boxplot=p1))
}
