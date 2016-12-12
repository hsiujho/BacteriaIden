#
#
#
#
#
#

phylo_ordination=function(phylo,ranklv,toRA=T
                          ,ord.method="PCoA",ord.distance="bray",group_var
                          ,plot.ellipse=T,ellipse.kind="sd",ellipse.conf=.95
                          ,plot.spider=T,rtest.nrepet=0,...
){
  if(!missing(ranklv)){
    if(any(rank_names(phylo)==ranklv)){
      phylo%<>%tax_glom(ranklv,NArm=F)%>>%(prune_taxa(taxa_sums(.)>0,.))
    }
  }
  if(toRA){
    phylo%<>%transform_sample_counts(function(x)x/sum(x)*100)
  }
  ord1=ordinate(phylo, ord.method, ord.distance,...)
  if(!missing(group_var)){
    if(any(sample_variables(phylo)==group_var)){
      gv=get_variable(phylo,group_var)
      if(!is.factor(gv)) gv%<>%factor()
    } else {
      gv=rep("Single",nsamples(phylo)) %>% factor()
    }
  } else {
    gv=rep("Single",nsamples(phylo)) %>% factor()
  }
  if(nlevels(gv)==1) rtest.nrepet=0
  ord2=plot_ordination(phylo,ord1,justDF=T)[,1:2]
  plot1=ordiplot(ord2,display="sites")
  ord3<-ordiellipse(plot1, gv
                    ,kind=ellipse.kind
                    ,conf=ellipse.conf)
  dev.off()
  #Monte-Carlo test on the between-groups inertia percentage.
  #http://pbil.univ-lyon1.fr/ade4/ade4-html/rtest.between.html
  if(rtest.nrepet>0){
    pvalue=rtest(bca(phyloseq::distance(phylo,ord.distance,...) %>%
                       dudi.pco(scannf=F),gv, scan = FALSE),nrepet=rtest.nrepet) %$%
      pvalue %>>% (ifelse(.<.001,"P<0.001"
                          ,ifelse(.>0.999,"P>0.999"
                                  ,sprintf("P=%.3f",.))))
    #    plot_title=sprintf("%s\np-value=%.5f using rtest",main,pvalue)
  } #else {
  #    plot_title=main
  #  }
  #generate the points of ellipse
  df_ell <- data.frame()
  for(g in names(ord3)){
    df_ell<-rbind(df_ell,data.frame(
      vegan:::veganCovEllipse(cov=ord3[[g]]$cov,center=ord3[[g]]$center
                              ,scale=ord3[[g]]$scale),group=g)
    )
  }
  colnames(df_ell)=c("Axis.1","Axis.2","label")
  if(nlevels(gv)>1){
    p1=plot_ordination(phylo, ord1, type = "samples", color = group_var)
  } else {
    sample_data(phylo)$Single=gv
    p1=plot_ordination(phylo, ord1, type = "samples", color = "Single")+ theme(legend.position="none")
  }
  if(plot.ellipse==T){
    p1=p1+geom_path(data=df_ell, aes(x=Axis.1, y=Axis.2,colour=label), size=1, linetype=1)
  }
  muS=by(ord2,gv,function(x){
    list(mu=colMeans(x),S=cov(x))
  })
  spider=lapply(1:nsamples(phylo),function(x){
    c(ord2[x,1:2],muS[[gv[x]]]$mu)
  }) %>% do.call(rbind,.) %>% data.frame(stringsAsFactors=F)
  for(i in 1:4) spider[[i]]=unlist(spider[[i]])
  spider$label=gv
  if(plot.spider==T){
    p1=p1+geom_segment(mapping=aes(x = Axis.1, y = Axis.2, xend = Axis.1.1, yend = Axis.2.1, colour=label),data=spider)
  }
  if(rtest.nrepet>0){
    xlim=ggplot_build(p1)$panel$ranges[[1]]$x.range
    ylim=ggplot_build(p1)$panel$ranges[[1]]$y.range
    p1=p1+annotate("text",x=xlim[1],y=ylim[2],label=pvalue,hjust=0)
  }
  return(p1)
}
