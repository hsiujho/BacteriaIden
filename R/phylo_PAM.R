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

phylo_PAM=function(phylo,ranklv="Genus",NArm = F
                   ,mothur_path="E:/soft/Mothur.win_64_1.39.5/mothur",Num_thread=1){

  c0=my_tax_glom(phylo,ranklv,NArm = NArm)
  if(!NArm){
    c0=add_label_to_unknown_taxon(c0)
  }
  c1=otu_table(c0)@.Data %>>% t()

  t0=tax_table(c0)@.Data %>>%
    data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("OTUID")

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

#PAM

  origin_wd=getwd()
  setwd(mothur_path)
  system('mothur "#get.communitytype(shared=dmm.shared,method=pam)"',intern = T)
  setwd(origin_wd)

  g0=sort(list.files(mothur_path,pattern="\\.pam\\.\\d\\.mix\\.posterior",full.names = T))
  g1=sort(list.files(mothur_path,pattern="\\.pam\\.\\d\\.mix\\.posterior",full.names = F))
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
      mutate_each(funs(as.character),OTUID,Sample)

  relabund=lapply(g2,function(e1){
    # e1=g2$dmm.0.05.pam.3.mix.posterior
    e2=colnames(e1) %>>%
      (regmatches(., regexpr("Partition_\\d", .)))
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
      (regmatches(., regexpr("Partition_\\d", .)))
    i1$Difference=abs(apply(i1[,i2,drop=F],1,max)-apply(i1[,i2,drop=F],1,min))
    i3=arrange(i1,desc(Difference))
    return(i3)
  }) %>>% setNames(sub("posterior","relabund",names(g2)))

  ##移除mothur生成的輸出資料
  list.files(mothur_path,pattern="\\.pam\\.\\d\\.mix\\.posterior",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.pam\\.mix\\.",full.names = T) %>>%
    file.remove()
  list.files(mothur_path,pattern="\\.logfile",full.names = T) %>>%
    file.remove()
  file.path(mothur_path,"dmm.shared") %>>% file.remove()

  ##
  return(list(fit=h2,posterior=g2,relabund=relabund))
}
