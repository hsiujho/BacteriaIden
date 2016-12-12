myord=function(
  phylo=a2 ##phyloseq的物件
  ,gv.name="GF.SA.ST" ##在物件sample_data中的欄位名稱
  ,ord.method="PCoA" ##方法, c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")其中一項
  ,ord.distance="bray" ##距離, 可用distance("list") or distanceMethodList查詢
  ,plot.ellipse=T ##是否畫橢圓線
  ,ellipse.kind="sd"
  ,ellipse.conf=.95
  ,plot.spider=T ##是否畫蜘蛛線
  ,main="samples"
  ,coord.fixed=F
  ,legend.name="Group"
  ,rtest.nrepet=0 #number of repetation#
  ,... ##計算unifrac會用到其他參數
){
  gv=get_variable(phylo,gv.name)
  if(!is.factor(gv)) gv%<>%factor()
  if(length(unique(gv))==1) rtest.nrepet=0
  ord1=ordinate(phylo, ord.method, ord.distance,...)
  ord2=plot_ordination(phylo,ord1,justDF=T)[,1:2]
  plot1=ordiplot(ord2,display="sites")
  ord3<-ordiellipse(plot1, gv
                    ,kind=ellipse.kind
                    ,conf=ellipse.conf)
  dev.off()
  #Monte-Carlo test on the between-groups inertia percentage.
  #http://pbil.univ-lyon1.fr/ade4/ade4-html/rtest.between.html
  if(rtest.nrepet>0){
    pvalue=rtest(bca(phyloseq::distance(phylo,ord.distance,...) %>% dudi.pco(scannf=F)
                     ,gv, scan = FALSE),nrepet=rtest.nrepet) %$% pvalue
    plot_title=sprintf("%s\np-value=%.5f using rtest",main,pvalue)
  } else {
    plot_title=main
  }

  #產生橢圓等高線座標
  df_ell <- data.frame()
  for(g in names(ord3)){
    df_ell<-rbind(df_ell,data.frame(
      vegan:::veganCovEllipse(cov=ord3[[g]]$cov,center=ord3[[g]]$center
                              ,scale=ord3[[g]]$scale),group=g)
    )
  }
  colnames(df_ell)=c("Axis.1","Axis.2","label")
  p1=plot_ordination(phylo, ord1, type = "samples", color = gv.name, title = plot_title)
  #+scale_colour_brewer(palette="Dark2",name=legend.name)
  if(coord.fixed){
    p1 = p1+coord_fixed()
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
  return(p1)
}
