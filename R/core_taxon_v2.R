#
#
#
#
#
#
# data("GlobalPatterns")
# phylo=my_tax_glom(GlobalPatterns,"Family")
# c0=core_taxon_v2(phylo,sampling=T,Num_samp=seq(1,26,5))
# c0$sampling+ylab("Number of family")

core_taxon_v2=function(phylo,prop=0.95,sampling=F,Num_samp,Num_rep=100
                    ,Num_Cores=max(1,detectCores()-1)){

  n=nsamples(phylo)
  h0=phylo@otu_table@.Data %>>% {
    if (taxa_are_rows(phylo)) {
      .
    } else {
      t(.)
    }
  }
  h4=tax_table(phylo)[rowSums(otu_table(phylo)>0)>prop*n,]

  if(sampling){

    h10=sample_names(phylo)
    if(missing(Num_samp)){
      Num_samp=unique(floor(exp(seq(0,log(n),length=100)))) %>>%
        head(-1) %>>%
        c(n)
    }

    require(foreach)
    require(iterators)
    library(doParallel)
    registerDoParallel(cores=Num_Cores)

    h11=foreach(i=Num_samp,.combine = rbind) %dopar% {
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
      xlab("Number of samples")+ylab("Number of taxon")+
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
