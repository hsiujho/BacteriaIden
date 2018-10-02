#
#
#
#
#
#

phylo_dendro_bar=function(phylo,group_var,ranklv="Phylum",topn=6,xlim,ylim,legend_ncol=2,method="complete"){
  if(any(rank_names(phylo)==ranklv)){
    phylo%<>%my_tax_glom(ranklv,NArm=F)%>>%(prune_taxa(taxa_sums(.)>0,.))
  }
  phylo%<>%transform_sample_counts(function(x)x/sum(x)*100)
  hc <- otu_table(phylo)@.Data %>% t() %>% vegdist() %>% hclust(method=method)

  b0=taxa_sums(phylo) %>>% sort.default(decreasing=T) %>>% names() %>>% "["(1:topn)
  a1=prune_taxa(taxa_names(phylo)%in%b0,phylo)

  rankv=rank_names(phylo)[1:which(rank_names(phylo)==ranklv)]
  a2=psmelt(a1) %>>% mutate_each_(funs(as.character),rankv)

  for(i in 2:length(rankv)){
    j=rankv[i]
    j1=rankv[i-1]
    j2=sprintf(";%s_",tolower(substr(j,1,1)))
    a2[[j]]=ifelse(a2[[j]]=="",paste0(a2[[j1]],j2),a2[[j]])
  }
  a5=a2[,rankv] %>>% unique()
  a2 %<>% dcast(formula(sprintf("%s~Sample",ranklv)),value.var="Abundance")
  a3=as.matrix(a2[,-1])
  rownames(a3)=a2[,1]
  a4=a3[order(rowSums(a3),decreasing = T),]
  colorv=c(brewer.pal(8,"Set1"),brewer.pal(8,"Set2"),brewer.pal(11,"Set3"))

  if(missing(xlim)) xlim=c(0,210)
  if(missing(ylim)) ylim=c(-topn%/%legend_ncol-1,ncol(a4))

  if(!missing(group_var)){
    if(any(group_var==sample_variables(phylo))){
      if(is.factor(get_variable(phylo,group_var))){
        a6=get_variable(phylo,group_var) %>>% data.frame() %>>% setNames(group_var)
        rownames(a6)=sample_names(phylo)
        #        print(a6)
        a7=a6[colnames(a4),]
        c0=ape::as.phylo(hc)
        c0$edge.length=c0$edge.length/max(c0$edge.length)*50

        par(mar=c(1,1,1,1),xpd=NA)
        plot(c0,x.lim=xlim,y.lim=ylim, tip.color =
               as.numeric(a7))
        barplot(a4[,hc$order],col=colorv[1:min(topn,length(colorv))]
                ,horiz=T,yaxt="n",add=T,offset = xlim[2]-100,
                width=0.8,space=c(.75,rep(0.25,ncol(a4)-1))
                ,legend.text = rownames(a4),xaxt="n"
                ,args.legend = list(x = xlim[2],y=0,horiz=F,ncol=legend_ncol,bty="n"))
        legend("topleft",levels(a7),text.col=1:nlevels(a7),bty="n")

        return(list(RA=a4,tax_anno=a5,group=a6))
      }
    }
  }

  par(mar=c(1,1,1,1),xpd=NA)
  plot(ape::as.phylo(hc),x.lim=xlim,y.lim=ylim)
  barplot(a4[,hc$order],col=colorv[1:min(topn,length(colorv))]
          ,horiz=T,yaxt="n",add=T,offset = xlim[2]-100,
          width=0.8,space=c(.75,rep(0.25,ncol(a4)-1))
          ,legend.text = rownames(a4),xaxt="n"
          ,args.legend = list(x = xlim[2],y=0,horiz=F,ncol=legend_ncol,bty="n"))

  return(list(RA=a4,tax_anno=a5))
}
