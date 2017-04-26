#
#
#
#
#
#

core_taxon=function(phylo,prop=0.95,ranklv="Genus",NArm=T
                    ,sampling=F,Num_rep=100,Num_Cores=max(1,detectCores()-1)){
  if(!missing(ranklv)){
    if(any(rank_names(phylo)==ranklv)){
      phylo%<>%tax_glom(ranklv,NArm=NArm)%>>%(prune_taxa(taxa_sums(.)>0,.))
    }
  }
  if(NArm==F){
    phylo=add_label_to_unknown_taxon(phylo)
  }

  n=nsamples(phylo)
  h4=tax_table(phylo)[rowSums(otu_table(phylo)>0)>prop*n,]

  if(sampling){
    h0=phylo@otu_table@.Data
    h10=sample_names(phylo)

    require(foreach)
    require(iterators)
    library(doParallel)
    registerDoParallel(cores=Num_Cores)

    h11=foreach(i=1:n,.combine = rbind) %dopar% {
      threshold=prop*i
      outp=matrix(NA,Num_rep,4)
      for(j in 1:Num_rep){
        h1=rowSums(h0[,sample(h10,i), drop=FALSE]>0)
        outp[j,]=c(i,j,sum(h1!=0),sum(h1>threshold))
      }
      return(outp)
    }
    stopImplicitCluster()
    colnames(h11)=c("Num_Sample","rep","Total","Core")

    h2=data.frame(h11) %>>%
      melt(id.vars=c("Num_Sample","rep"),measure.vars=c("Total","Core"))

    h3=ggplot(h2,aes(x=factor(Num_Sample),y=value))+
      geom_boxplot()+
      xlab("Number of samples")+ylab("Number of genera")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
      # facet_wrap(~variable,nrow=2,scales = "free_y")
      facet_grid(variable~.,scales = "free_y")+
      theme(legend.background=element_blank(),
            plot.background = element_blank())
    return(list(core_taxon=h4,sampling=h3))
  } else {
    return(list(core_taxon=h4))
  }
}
