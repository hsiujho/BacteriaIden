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
  if(taxa_are_rows(phylo)){
    varnames=c("OTU","Sample")
  } else {
    varnames=c("Sample","OTU")
  }
  # u0=phylo
  # u0@sam_data<-NULL
  # u1=psmelt(u0) %>>%
  #   select_(.dots=c("OTU","Sample","Abundance",keep_rank)) %>>%
  #   data.table()
  v0=otu_table(phylo) %>>%
    melt(varnames=varnames,value.name="Abundance") %>>%
    filter(Abundance>0) %>>%
    mutate_at(funs(as.character),.vars=c("OTU","Sample")) %>>%
    mutate(OTU=factor(OTU,levels=taxa_names(phylo))
           ,Sample=factor(Sample,levels=sample_names(phylo)))
  v1=tax_table(phylo) %>>% data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("OTU") %>>%
    mutate(OTU=factor(OTU,levels=taxa_names(phylo)))
  #每組taxid以加總Sample後reads最多的來取代, 可能有些情況需再釐清, (1)跳層表現量最多變成另一個ID{不會}, (2)若有相同表現量應取哪個ID(第一次出現那個)
  u1=left_join(v0,v1,by="OTU") %>>%
    select_(.dots=c("OTU","Sample","Abundance",keep_rank)) %>>%
    data.table()
  u11=u1[,.(Abundance=sum(Abundance)),by=c("OTU",keep_rank)]
  setkey(u11,OTU)
  u12=u11[u11[,.I[which.max(Abundance)],by=keep_rank]$V1]
  u13=left_join(select(u1,-OTU),select(u12,-Abundance),by=keep_rank)
  u1=as.data.table(u13)
  #統一相同taxon的編號, 指定by且使用:=操作, 資料不會化簡
  if(NArm){
    u1=u1[u1[[ranklv]]!=""]
  }
  u2=u1[,.(Abundance=sum(Abundance)),by=c("Sample","OTU",keep_rank)]
  #指定by且使用.()操作, 資料會化簡
  if(taxa_are_rows(phylo)){
    u3=dcast(u2,OTU~Sample,value.var = "Abundance",fill=0) %>>%
      column_to_rownames("OTU") %>>%
      data.matrix() %>>%
      otu_table(taxa_are_rows=T)
  } else {
    u3=dcast(u2,Sample~OTU,value.var = "Abundance",fill=0) %>>%
      column_to_rownames("Sample") %>>%
      data.matrix() %>>%
      otu_table(taxa_are_rows=F)
  }
  u4=unique(u2,by=c("OTU",keep_rank)) %>>%
    (.[,c("Sample","Abundance"):=NULL]) %>>%
    data.frame(stringsAsFactors=F) %>>%
    mutate_at(funs(as.character),.vars=keep_rank) %>>%
    column_to_rownames("OTU") %>>%
    as.matrix() %>>%
    tax_table()
  u5=phyloseq(u3,u4)
  if(!is.null(phylo@phy_tree)){
    phy_tree(u5)=phy_tree(phylo)
  }
  if(!is.null(phylo@sam_data)){
    sample_data(u5)=sample_data(phylo)
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
