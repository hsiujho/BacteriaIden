#
#
#
#
#
#

phylo_rarecurve=function(phylo,group.varn,ranklv="Phylum",facet.row=".",facet.col="."
                         , rarecurve.step = 20, rarecurve.sample = NA){

  if(all(sample_variables(phylo)!=group.varn)|is.na(group.varn)){
    sample_data(phylo)$dummy_var=""
    group.varn="dummy_var"
  }

  if(all(sample_variables(phylo)!=facet.row)){
    facet.row="."
  }

  if(all(sample_variables(phylo)!=facet.col)){
    facet.col="."
  }

  if(any(rank_names(phylo)==ranklv)) {
    a2=phylo %>% tax_glom(ranklv,NArm=F)
  } else {
    a2=phylo
    ranklv="OTU"
  }

  sample_data(a2)$SampleID=sample_names(a2)
  if(group.varn=="dummy_var"){
    sample_data(a2)$colour.var=""
  } else {
    sample_data(a2)$colour.var=get_variable(a2,group.varn)
  }

  b1=otu_table(a2)@.Data%>%t()
  S<-specnumber(b1)
  if(is.na(rarecurve.sample))
    rarecurve.sample <- min(rowSums(b1))
  b2=rarecurve(b1, step = rarecurve.step, sample = rarecurve.sample, col = "blue", cex = 0.6)
  dev.off()
#  b5=rarefy(b1[1,], seq(10000,30000,5000)) #輸出是期望值, 所以不用重覆
#  b6=sapply(1:1000,function(i){
#    rep(colnames(b1),b1[1,]) %>% sample(30000) %>% unique() %>% length()
#  }) %>% mean()

  b3=lapply(1:length(b2),function(i){
    b2[[i]]%>%data.frame(y=.)%>%rownames_to_column()%>%
      mutate(x=sub("N","",rowname)%>%as.numeric,SampleID=rownames(b1)[i],rowname=NULL)
  }) %>% bind_rows()

  varn=c("SampleID","colour.var")
  if(facet.row!=".") varn=c(varn,facet.row)
  if(facet.col!=".") varn=c(varn,facet.col)
  b3%<>%left_join(get_variable(a2,varn),by="SampleID")

  if(group.varn=="dummy_var"){
    p1=ggplot(b3,aes(x,y,group=SampleID))
  } else {
    p1=ggplot(b3,aes(x,y,group=SampleID,colour=colour.var))
  }
  p1=p1+geom_line()+
    ggtitle(sprintf("Rarefaction curves"))+
    ylab(sprintf("Rarefied No. of %s",ranklv))

  if(group.varn!="dummy_var")
    p1=p1+scale_color_discrete(guide = guide_legend(title = group.varn))

  if(facet.row!="."|facet.col!=".")
    p1=p1+facet_grid(formula(sprintf("%s~%s",facet.row,facet.col)))

  return(p1)
}
