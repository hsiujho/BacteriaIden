#
#
#
#
# my_tax_glom第二版
# 主要是簡化coding
#
#

my_tax_glom_v2=function(phylo,ranklv,NArm=F){
  keep_rank=rank_names(phylo) %>>%
    "=="(ranklv) %>>%
    which() %>>%
    (head(rank_names(phylo),.))
  
  dt1=phylo@otu_table@.Data %>>% {
    if(taxa_are_rows(phylo)){
      .
    } else {
      t(.)
    }
  } %>>%
    (merge(phylo@tax_table@.Data,.,by=0)) %>>%
    data.table()
  dt1[,id:=factor(Row.names,levels=taxa_names(phylo))]
  setkey(dt1,id)
  # setkeyv(dt1,cols=c("id",keep_rank))
  dt1[, Sum := Reduce(`+`, .SD), .SDcols=sample_names(phylo)]
  dt2=dt1[dt1[, .I[which.max(Sum)],by=keep_rank]$V1,c("Row.names",keep_rank),with=F]
  dt2[,OTUID:=factor(Row.names,levels=Row.names)]
  dt3=merge(dt1,dt2,by=keep_rank)
  dt4=dt3[, lapply(.SD,sum),by="OTUID", .SDcols=sample_names(phylo)]
  setkey(dt4,OTUID)
  
  otu=as.matrix(dt4[,sample_names(phylo),with=F])
  rownames(otu)=as.character(dt4$OTUID)
  
  tax=as.matrix(dt2[,keep_rank,with=F])
  rownames(tax)=as.character(dt2$OTUID)
  
  phylo2=phyloseq(otu_table(otu,taxa_are_rows = T),tax_table(tax))
  if(!is.null(phylo@phy_tree)){
    phy_tree(phylo2)=phy_tree(phylo)
  }
  if(!is.null(phylo@sam_data)){
    sample_data(phylo2)=sample_data(phylo)
  }
  return(phylo2)
}

### testing

function(){
  
  require(BacteriaIden)
  require(microbenchmark)
  
  
  data(GlobalPatterns)
  phylo  <- GlobalPatterns
  
  a1=my_tax_glom(phylo,"Class")
  a2=tax_glom(phylo,"Class",NArm=F)
  all.equal(taxa_names(a1),taxa_names(a2))
  
  # microbenchmark(my_tax_glom(GP,"Class"),tax_glom(GP,"Class",NArm=F),  times = 10)
  # 
  # profvis::profvis(
  #   my_tax_glom(GP,"Class")
  # )
  # 
  # profvis::profvis(
  #   tax_glom(GP,"Class",NArm=F)
  # )
  
  
  ####
  
  a3=my_tax_glom_v2(phylo,"Class")
  
  merge(data.frame(V1=sample_sums(a1))
        ,data.frame(V2=sample_sums(a3)),by=0) %$%
    all(V1==V2)
  
  merge(data.frame(V1=taxa_sums(a1))
        ,data.frame(V2=taxa_sums(a3)),by=0) %$%
    all(V1==V2)

  all.equal(a1@otu_table@.Data,a3@otu_table@.Data[taxa_names(a1),])
  all.equal(a1@tax_table@.Data,a3@tax_table@.Data[taxa_names(a1),])
  
  # microbenchmark(my_tax_glom(phylo,"Class"),my_tax_glom_v2(phylo,"Class"),  times = 10)
  # 效率稍減, 但增加code可讀性
  # profvis::profvis(
  #   new_tax_glom(GP,"Class")
  # )
  
}


