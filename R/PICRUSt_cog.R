#
#
#
#
#
#

PICRUSt_cog=function(phylo,path="E:/fastq/picrust/Precalculated"){
  fp1=tax_table(phylo)@.Data %>>% as.data.frame(stringsAsFactors=F) %>>% tibble::rownames_to_column("#OTU_IDs") %>>% left_join(GCN,by="#OTU_IDs")
  for(i in 1:nsamples(phylo)){
    otu_table(phylo)[,i]=otu_table(phylo)[,i]/fp1$`16S_rRNA_Count`
  }
  readfile=function(fn,header=T){
    read.big.matrix(fn,"\t",header=header,type="double") %>>% (.[.[,1]%in%taxa_names(phylo),])
  }

  fn1=list.files(path,"cogf_a",full.names = T)
  b0=rbind(readfile(fn1[1]),do.call(rbind,lapply(fn1[-1],readfile,header=F)))
  b1=b0[,which(regexpr("^COG",colnames(b0))!=-1)]
  rownames(b1)=b0[,1]
  b2=b1[taxa_names(phylo),]

  b3=t(b2) %*% otu_table(phylo)@.Data
  c1=otu_table(b3,T)
  c2=tax_table(cog_anno)
  c3=sample_data(phylo)
  fp_cog=phyloseq(c1,c2,c3)
  return(list(trans=b2,fp=fp_cog))
}
