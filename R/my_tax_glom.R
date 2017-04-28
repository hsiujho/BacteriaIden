#
#
#
#利用data.table加速tax_glom
#
#
#


my_tax_glom=function(phylo,ranklv,NArm=F){
  num_ranklv=which(rank_names(phylo)==ranklv)
  keep_rank=rank_names(phylo)[1:num_ranklv]
  # u0=phylo
  # u0@sam_data<-NULL
  # u1=psmelt(u0) %>>%
  #   select_(.dots=c("OTU","Sample","Abundance",keep_rank)) %>>%
  #   data.table()
  v0=otu_table(phylo) %>>%
    melt(varnames=c("OTU","Sample"),value.name="Abundance") %>>%
    mutate_each(funs(as.character),OTU,Sample)
  v1=tax_table(phylo) %>>% data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("OTU")
  u1=left_join(v0,v1,by="OTU") %>>%
    select_(.dots=c("OTU","Sample","Abundance",keep_rank)) %>>%
    data.table()
  if(NArm){
    u1=u1[u1[[ranklv]]!=""]
  }
  u2=u1[,.(Abundance=sum(Abundance),OTU=min(OTU)),by=c("Sample",keep_rank)]
  u3=dcast(u2,OTU~Sample,value.var = "Abundance",fill=0) %>>%
    column_to_rownames("OTU") %>>%
    data.matrix() %>>%
    otu_table(taxa_are_rows=T)
  u4=unique(u2,by=c("OTU",keep_rank)) %>>%
    (.[,c("Sample","Abundance"):=NULL]) %>>%
    data.frame(stringsAsFactors=F) %>>%
    mutate_each_(funs(as.character),vars=keep_rank) %>>%
    column_to_rownames("OTU") %>>%
    as.matrix() %>>%
    tax_table()
  u5=phyloseq(u3,sample_data(phylo),u4)
  return(u5)
}



# k0=subset_samples(MS16049_gg_phylo,group=="ST") %>>%
#   (prune_taxa(taxa_sums(.)>0,.))
#
# system.time(
#   u6<-my_tax_glom(k0,"Genus")
# )
# system.time(
#   k2<-tax_glom(k0,"Genus",NArm=F)
# )
# system.time(
#   u7<-my_tax_glom(k0,"Genus",NArm=T)
# )
# system.time(
#   k3<-tax_glom(k0,"Genus",NArm=T)
# )
