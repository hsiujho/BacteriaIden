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
    filter(Abundance>0) %>>%
    mutate_each(funs(as.character),OTU,Sample)
  v1=tax_table(phylo) %>>% data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("OTU")
  u1=left_join(v0,v1,by="OTU") %>>%
    select_(.dots=c("OTU","Sample","Abundance",keep_rank)) %>>%
    data.table() %>>%
    (.[,OTU:=min(OTU),by=keep_rank])
  #統一相同taxon的編號, 指定by且使用:=操作, 資料不會化簡
  if(NArm){
    u1=u1[u1[[ranklv]]!=""]
  }
  u2=u1[,.(Abundance=sum(Abundance)),by=c("Sample","OTU",keep_rank)]
  #指定by且使用.()操作, 資料會化簡
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
  if(is.null(phylo@phy_tree)){
    u5=phyloseq(u3,sample_data(phylo),u4)
  } else {
    u5=phyloseq(u3,sample_data(phylo),u4,phy_tree(phylo))
    #因為使用相異的OTUID, phy_tree的部分會不同於原始tax_glom的篩選
  }
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
