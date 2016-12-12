#
# 合併並篩選由ClustIden_Blast產生的檔案及註解檔
# 輸出是phyloseq物件
#
#
#
#

blast_pick_anno=function(
  file_path  #由ClustIden_Blast產生的檔案, rds格式限定
  ,anno_path #註解檔的路徑, rds格式限定
){

  if(!file.exists(anno_path)) stop(sprintf("%s does not exist.",anno_path))
  s1=file_path %>% "["(!file.exists(.))
  if(length(s1)>0) warning(sprintf("These files do not exist:\n%s",paste0(s1,collapse = "\n")))
  if(length(s1)==length(file_path)) stop("Nothing can be done")

  #篩選pident+qcovs>=186, 並依pident給定層級深度

  b0=file_path
  on=Sys.time()
  b1=lapply(b0,function(x){
    read_rds(x) %>% filter(pident+qcovs>=186) %>%
      mutate(Deep=cut(pident,breaks=c(85,91,92,95,97,99,Inf)
                      ,labels=c("Class","Order","Family","Genus","Species","strain")
                      ,include.lowest =T,ordered_result=F)
             ,n=qseqid%>%sub(".*;size=","",.)%>%as.integer()) %>%
      group_by(sgi,Deep) %>% summarise(n=sum(n)) %>% arrange(sgi,Deep) %>%
      mutate(SampleID=x%>%sub(".*/","",.)%>%sub(".rds","",.))
  }) %>% bind_rows() %>% dcast(sgi+Deep~SampleID,fill=0,value.var="n") %>%
    mutate_each(funs(as.character),sgi)
  Sys.time()-on

  b2=read_rds(anno_path) %>%
    select(gi,superkingdom,phylum,class,order,family,genus,species,organism,strain) %>%
    dplyr::rename(Kingdom=superkingdom,Phylum=phylum,Class=class,Order=order,Family=family,Genus=genus,Species=species)

  b3=left_join(select(b1,sgi,Deep), b2, by=c("sgi"="gi"))
  k0=c("Class","Order","Family","Genus","Species")
  for(i in 1:4) b3[b3$Deep==k0[i],k0[(i+1):5]]=NA

  b11=as.matrix(b1[,-(1:2)])
  rownames(b11)=paste0(b1[,1],"_",substr(b1[,2],1,1))

  b31=as.matrix(b3[,-(1:2)])
  rownames(b31)=paste0(b3[,1],"_",substr(b3[,2],1,1))

  phylo=phyloseq(otu_table(b11,taxa_are_rows=T),tax_table(b31))
  return(phylo)
}
