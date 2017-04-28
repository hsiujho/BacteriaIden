#
#
#以mothur進行DMM模型的事後分組
#
#
#

# path="E:/Biom/HM16STR"
# savepath=file.path(path,"greengene_20170427")
# phylo=readRDS(file.path(savepath,"HM16STR_gg_phylo.rds"))
#
# a0=phylo_DMM(phylo)

phylo_DMM=function(phylo,ranklv="Genus",NArm = F
                   ,mothur_path="E:/soft/Mothur.win_64_1.39.5/mothur",Num_thread=1){

  c0=tax_glom(phylo,ranklv,NArm = NArm)
  if(!NArm){
    c0=add_label_to_unknown_taxon(c0)
  }
  c1=otu_table(c0)@.Data %>>% t()

  shared.df<-as.data.frame(c1)
  colnames(shared.df)<-paste("gg",colnames(shared.df),sep=".")
  numOtu<-data.frame(numOtu=rep(ncol(shared.df),nrow(shared.df)))
  Group<-data.frame(Group=rownames(shared.df))
  label<-data.frame(label=rep(0.05,nrow(shared.df)))
  out.tsv<-data.frame(label,Group,numOtu,shared.df)

  write.table(out.tsv,file=file.path(mothur_path,"dmm.shared"),sep="\t",row.names=FALSE,quote=FALSE)

  #餵給get.communitytype
  #https://www.mothur.org/wiki/Get.communitytype

  #E:\soft\Mothur.win_64_1.39.5\mothur
  #執行 mothur > get.communitytype(shared=mouse.fmt.genus.shared)

  origin_wd=getwd()
  setwd(mothur_path)
  system(sprintf('mothur "#get.communitytype(shared=dmm.shared, processors=%i)"',Num_thread),intern = T)
  setwd(origin_wd)

  t0=tax_table(c0)@.Data %>>%
    data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("OTUID")

  #DMM結果轉成RDS格式檔

  d0=sort(list.files(mothur_path,pattern="\\.mix\\.posterior",full.names = T))
  d1=sort(list.files(mothur_path,pattern="\\.mix\\.posterior",full.names = F))
  d2=lapply(d0,function(x){
    read.table(x,header = T,stringsAsFactors=F)
  }) %>>% setNames(d1)

  e0=sort(list.files(mothur_path,pattern="\\.mix\\.relabund",full.names = T))
  e1=sort(list.files(mothur_path,pattern="\\.mix\\.relabund",full.names = F))
  e2=lapply(e0,function(x){
    e3=read.table(x,header = T,stringsAsFactors=F)
    e4=grepl("Partition_\\d_\\d+", colnames(e3))
    e5=regmatches(colnames(e3)[e4], regexpr("Partition_\\d", colnames(e3)[e4]))
    colnames(e3)[e4]=e5
    return(e3)
  }) %>>% setNames(e1)

  f0=sort(list.files(mothur_path,pattern="\\.dmm\\.mix\\.",full.names = T))
  f1=sort(list.files(mothur_path,pattern="\\.dmm\\.mix\\.",full.names = F))
  f2=lapply(f0,function(x){
    f3=read.table(x,header = !grepl("\\.design", x),stringsAsFactors=F)
    if(grepl("\\.design", x)) {
      colnames(f3)=c("Sample","CommunityType")
    }
    return(f3)
  }) %>>% setNames(f1)

  ##移除mothur生成的輸出資料

  list.files(mothur_path,pattern="\\.mix\\.posterior",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.mix\\.relabund",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.dmm\\.mix\\.",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.logfile",full.names = T) %>>%
    file.remove()
  file.path(mothur_path,"dmm.shared") %>>% file.remove()

  ##
  return(list(fit=f2,posterior=d2,relabund=e2,taxon=t0))
}
