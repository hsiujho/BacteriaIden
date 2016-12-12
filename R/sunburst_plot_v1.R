#
#
#
#
#
#

sunburst_plot_v1=function(phylo#,ranklv=rank_names(phylo)[2:6]
                          ,xshift=2,color.threshold=1){
  if(nsamples(phylo)!=1) stop("nsamples > 1")
  a5=phylo %>>%
    (prune_taxa(taxa_sums(.)>0,.)) %>>%
    transform_sample_counts(function(x)x/sum(x)*100)
  a5@sam_data=NULL
  k=6
  b0=tax_glom(a5,"Genus",NArm=F) %>>% transform_sample_counts(function(x)x/sum(x)*100) %>>% psmelt() %>>% arrange(Kingdom,Phylum,Class,Order,Family,Genus) %>>% mutate(ymax=cumsum(Abundance),ymin=c(0,head(ymax,-1)),xmin=k-1,xmax=k,fill_var=as.character(Genus))
  k=5
  b1=group_by(b0,Kingdom,Phylum,Class,Order,Family) %>>% summarise(ymax=max(ymax),ymin=min(ymin),xmin=k-1,xmax=k) %>>% mutate(fill_var=as.character(Family))
  k=4
  b2=group_by(b1,Kingdom,Phylum,Class,Order) %>>% summarise(ymax=max(ymax),ymin=min(ymin),xmin=k-1,xmax=k) %>>% mutate(fill_var=as.character(Order))
  k=3
  b3=group_by(b2,Kingdom,Phylum,Class) %>>% summarise(ymax=max(ymax),ymin=min(ymin),xmin=k-1,xmax=k) %>>% mutate(fill_var=as.character(Class))
  k=2
  b4=group_by(b3,Kingdom,Phylum) %>>% summarise(ymax=max(ymax),ymin=min(ymin),xmin=k-1,xmax=k) %>>% mutate(fill_var=as.character(Phylum))
  k=1
  b5=group_by(b4,Kingdom) %>>% summarise(ymax=max(ymax),ymin=min(ymin),xmin=k-1,xmax=k) %>>% mutate(fill_var=as.character(Kingdom))

  a7=bind_rows(list(
    #Kingdom=b5,
    Phylum=b4,Class=b3,Order=b2,Family=b1,Genus=b0),.id="Rank")

  color=brewer.pal(12, "Paired")

  a8=filter(a7,fill_var!='') %>>%
    mutate(Rank=factor(Rank,levels=rank_names(a5)[1:7])
           ,fill_var=factor(fill_var,levels=fill_var[!duplicated(fill_var)])) %>>%
    mutate(rn=as.integer(Rank)) %>>%
    group_by(Rank) %>>%
    mutate(seq=row_number(),l=n())
#  a8$col= apply(a8[,c("rn","l","seq")],1,function(v){
#    colorRampPalette(color[pmax(1,v[1]-1):pmin(v[1]+1,7)])(v[2])[v[3]]
#  })

  color11=brewer.pal(9, "Pastel1")[1]

  a8$col=ifelse(a8$ymax-a8$ymin>color.threshold,"",color11)
  a8$col[a8$col!=color11]=color

  maxx=max(a8$xmax)

  #  xshift=2,alpha(col,1-Rank*.1)
  p8=ggplot(a8, aes(xmin = xmin+xshift, xmax = xmax+xshift, ymin = ymin, ymax = ymax)) +
    geom_rect(aes(fill = col), colour = "grey50")+
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
