#
#
#利用
#http://enterotyping.embl.de/enterotypes.html
#的方法
#
#


phylo_PAM_v2=function(phylo,ranklv="Genus",NArm = F, NumComponents=1:10,Num_thread=1){
  c0=my_tax_glom(phylo,ranklv,NArm = NArm)
  if(!NArm){
    c0=add_label_to_unknown_taxon(c0)
  }
  c1=otu_table(c0)@.Data
  if(c0@otu_table@taxa_are_rows){
    c1=t(c1)
  }
  require(distRcpp)
  jsd=JSD(c1,taxa_are_rows = F,numThreads = Num_thread)
  pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(x, k, diss=TRUE)$clustering)
    return(cluster)
  }

  library(doParallel)
  library(foreach)
  registerDoParallel(cores=Num_thread) #設定執行緒數目

  fit <- foreach(k=NumComponents,.packages = "clusterSim") %dopar% {
    if (k==1) {
      return(list(CH=c(K=1,CH=0),WhichGroup=rep(1,nrow(c1))))
    } else {
      WhichGroup=pam.clustering(jsd, k)
      CH=c(K=k,CH=index.G1(c1, WhichGroup,  d = jsd,
                                 centrotypes = "medoids"))
      return(list(CH=CH,WhichGroup=WhichGroup))
    }
  }
  stopImplicitCluster()

  gm=sapply(fit,function(x) x$WhichGroup)
  colnames(gm)=sprintf("k_%02i",NumComponents)

  return(list(phylo=c0
              ,jsd=jsd
              ,CH=t(sapply(fit,function(x) x$CH))
              ,WhichGroup=gm))
}

# ref_path="E:/Biom/Report/000202/CO/20170713"
# data_path=file.path(ref_path,"data")
# phylo1=readRDS(file.path(data_path,"CO_gg_phylo_merge_20170713.rds"))
# phylo=subset_samples(phylo1,!sample_names(phylo1)%in%c("CO-2168-1","CO-2170-1","CO-2238-1","CO-2274-1","CO-2165-1","CO-2237-1")) %>>%
#   (prune_taxa(taxa_sums(.)>0,.))
# sample_names(phylo)=substr(sample_names(phylo),1,7)
#
# a0=phylo_PAM_v2(phylo,ranklv="Genus",NArm = F, NumComponents=1:5,Num_thread=5)
