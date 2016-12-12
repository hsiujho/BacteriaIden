#
#
#
#
#
#

sunburst_plot=function(phylo,ranklv=rank_names(phylo)[2:6],xshift=2){
  if(nsamples(phylo)!=1) stop("nsamples > 1")
  a5=phylo %>>%
    (prune_taxa(taxa_sums(.)>0,.)) %>>%
    transform_sample_counts(function(x)x/sum(x)*100)
  a5@sam_data=NULL
  a6=lapply(ranklv,function(x){
    tax_glom(a5,x,NArm=F) %>>% psmelt()
  })
  names(a6)=ranklv

  for(k in 1:length(ranklv)){
    rlk=ranklv[k]
    if(k>1){
      for(j in 1:(k-1)){
        rlj=ranklv[j]
        a6[[rlk]][[rlj]]%<>%factor(levels=levels(a6[[rlj]][[rlj]]))
      }
    }
    tx=sprintf("a6[[rlk]]%sarrange(%s,desc(Abundance))","%<>%",paste0(ranklv[1:(k-1)],collapse = ","))
    eval(parse(text=tx))
    a6[[rlk]][[rlk]]%<>%factor(levels=as.character(.[!duplicated(.)]))
    a6[[rlk]]%<>%mutate(ymax=cumsum(Abundance),ymin=c(0,head(ymax,-1)),xmin=k-1,xmax=k)
    a6[[rlk]]$fill_var=as.character(a6[[rlk]][[rlk]])
    #delete the unknown data
    a6[[rlk]]=a6[[rlk]][a6[[rlk]][[rlk]]!="",]
  }
  color=brewer.pal(7, "Dark2")

  a7=bind_rows(a6,.id="Rank") %>>%
    mutate(Rank=factor(Rank,levels=rank_names(a5)[1:7])
           ,fill_var=factor(fill_var,levels=fill_var[!duplicated(fill_var)]))

  a8=arrange(a7,Rank,desc(Abundance)) %>>%
    mutate(rn=as.integer(Rank)) %>>%
    group_by(Rank) %>>%
    mutate(seq=row_number(),l=n())
  a8$col= apply(a8[,c("rn","l","seq")],1,function(v){
    colorRampPalette(color[pmax(1,v[1]-1):pmin(v[1]+1,7)])(v[2])[v[3]]
  })
  maxx=max(a8$xmax)
  p8=ggplot(a8, aes(xmin = xmin+xshift, xmax = xmax+xshift, ymin = ymin, ymax = ymax,colour=alpha(col,1-Rank*.1))) +
    geom_rect(aes(fill = fill_var), colour = "grey50")+
    coord_polar(theta="y")+xlim(0,maxx+xshift)+
    theme_classic()+
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
#    theme(plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "in")) +
    theme(axis.ticks=element_blank())+theme(legend.position="none")
}

#p11=sunburst_plot(phylo=a5,ranklv=rank_names(a5)[1:7],xshift=0)
#p12=sunburst_plot(phylo=a5,ranklv=rank_names(a5)[1:6],xshift=0)
#p13=sunburst_plot(phylo=a5,ranklv=rank_names(a5)[2:6],xshift=2)
#p13+theme(plot.margin = unit(rep(-1,4), "cm"))
#p14=sunburst_plot(phylo=a5)
