#
#
#以DirichletMultinomial::dmn進行DMM模型的事後分組
#install
#source("https://bioconductor.org/biocLite.R")
#biocLite("DirichletMultinomial")
#

phylo_DMM_v2=function(phylo,ranklv="Genus",NArm = F, NumComponents=1:10,Num_thread=1){
  a0=my_tax_glom(phylo,ranklv,NArm = NArm)
  if(!NArm){
    a0=add_label_to_unknown_taxon(a0)
  }
  otu=otu_table(a0)@.Data
  if(phylo@otu_table@taxa_are_rows){
    otu=t(otu)
  }

  require(DirichletMultinomial)
  library(doParallel)
  library(foreach)
  registerDoParallel(cores=Num_thread) #設定執行緒數目

  fit <- foreach(k=NumComponents,.packages = "DirichletMultinomial") %dopar% {
    dmn(count=otu,k=k)
  }
  stopImplicitCluster()

  ModelIndex=sapply(fit,function(x){
    c(K=length(x@mixture$Weight),x@goodnessOfFit)
  }) %>>% t()
  (best <- fit[[which.min(ModelIndex[,"Laplace"])]])

  df1=apply(best@group,1,function(x){
    a0=which(x==max(x))
    ifelse(length(a0)==1,a0,sample(a0,1))
  }) %>>%
    (data.frame(partition=.)) %>>%
    rownames_to_column("Sample") %>>%
    mutate(partition=sprintf("pt_%02i",partition))

  p0 <- fitted(fit[[1]], scale=TRUE)     # scale by theta
  p4 <- fitted(best, scale=TRUE)
  colnames(p4) <- sprintf("pt_%02i", 1:ncol(p4))
  diff <- rowSums(abs(p4 - as.vector(p0)))
  o <- order(diff, decreasing=TRUE)
  cdiff <- cumsum(diff[o]) / sum(diff)
  df <- cbind(Mean=p0[o], p4[o,], diff=diff[o], cdiff) %>>%
    data.frame() %>>%
    rownames_to_column("taxon")
  df11=tax_table(a0) %>>%
    data.frame(stringsAsFactors=F) %>>%
    rownames_to_column("taxon") %>>%
    left_join(df,by="taxon") %>>%
    arrange(desc(diff))

  return(list(phylo=a0
       ,model=fit
       ,ModelIndex=ModelIndex
       ,best.fit=df11
       ,best.pt=df1))
}



# ref_path="E:/Biom/Report/000202/CO/20170713"
# data_path=file.path(ref_path,"data")
# phylo1=readRDS(file.path(data_path,"CO_gg_phylo_merge_20170713.rds"))
# phylo=subset_samples(phylo1,!sample_names(phylo1)%in%c("CO-2168-1","CO-2170-1","CO-2238-1","CO-2274-1","CO-2165-1","CO-2237-1")) %>>%
#   (prune_taxa(taxa_sums(.)>0,.))
# sample_names(phylo)=substr(sample_names(phylo),1,7)
#
# a0=phylo_DMM_v2(phylo,ranklv="Genus",NArm = F, NumComponents=1:3,Num_thread=3)
