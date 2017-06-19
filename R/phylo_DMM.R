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

  c0=my_tax_glom(phylo,ranklv,NArm = NArm)
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
    e4=grepl("Partition_\\d+_\\d+", colnames(e3))
    e5=regmatches(colnames(e3)[e4], regexpr("Partition_\\d+", colnames(e3)[e4]))
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

#PAM

  origin_wd=getwd()
  setwd(mothur_path)
  system('mothur "#get.communitytype(shared=dmm.shared,method=pam)"',intern = T)
  setwd(origin_wd)

  g0=sort(list.files(mothur_path,pattern="\\.pam\\.\\d+\\.mix\\.posterior",full.names = T))
  g1=sort(list.files(mothur_path,pattern="\\.pam\\.\\d+\\.mix\\.posterior",full.names = F))
  g2=lapply(g0,function(x){
    read.table(x,header = T,stringsAsFactors=F)
  }) %>>% setNames(g1)

  h0=sort(list.files(mothur_path,pattern="\\.pam\\.mix\\.",full.names = T))
  h1=sort(list.files(mothur_path,pattern="\\.pam\\.mix\\.",full.names = F))
  h2=lapply(h0,function(x){
    h3=read.table(x,header = !grepl("\\.design", x),stringsAsFactors=F)
    if(grepl("\\.design", x)) {
      colnames(h3)=c("Sample","Enterotype")
    }
    return(h3)
  }) %>>% setNames(h1)

  #計算各組菌叢比例平均數, 仿DMM的relabund

  c2=transform_sample_counts(c0,function(x)x/sum(x)*100)
  c3=melt(otu_table(c2),varnames=c("OTUID","Sample"),value.name = "Abundance") %>>%
      mutate_at(funs(as.character),.vars=c("OTUID","Sample"))

  relabund=lapply(g2,function(e1){
    # e1=g2$dmm.0.05.pam.3.mix.posterior
    e2=colnames(e1) %>>%
      (regmatches(., regexpr("Partition_\\d+", .)))
    i0=e1 %>>%
      rownames_to_column("Sample") %>>%
      melt(id.vars="Sample", measure.vars=e2, variable.name="Enterotype",value.name="Y") %>>%
      filter(Y==1) %>>%
      select(-Y)

    i1=c3 %>>%
      left_join(i0,by="Sample") %>>%
      group_by(OTUID,Enterotype) %>>%
      summarise(Mean=mean(Abundance)) %>>%
      dcast(OTUID~Enterotype,value.var="Mean") %>>%
      (left_join(t0,.,by="OTUID"))
    i2=colnames(i1) %>>%
      (regmatches(., regexpr("Partition_\\d+", .)))
    i1$Difference=abs(apply(i1[,i2,drop=F],1,max)-apply(i1[,i2,drop=F],1,min))
    i3=arrange(i1,desc(Difference))
    return(i3)
  }) %>>% setNames(sub("posterior","relabund",names(g2)))

  ##移除mothur生成的輸出資料
  list.files(mothur_path,pattern="\\.pam\\.\\d+\\.mix\\.posterior",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.pam\\.mix\\.",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.logfile",full.names = T) %>>%
    file.remove()
  file.path(mothur_path,"dmm.shared") %>>% file.remove()

  ##
  return(list(DMM=list(fit=f2,posterior=d2,relabund=e2,taxon=t0)
              ,PAM=list(fit=h2,posterior=g2,relabund=relabund)
         ))
}
