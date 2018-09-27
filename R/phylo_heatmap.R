#
#
#
#
#
#

phylo_heatmap=function(phylo,ranklv="Genus",topn=20,anno_var,...){
  a0=phylo %>>% transform_sample_counts(function(x)x/sum(x)*100) %>>% my_tax_glom(ranklv,NArm=F)
  b0=taxa_sums(a0) %>>% sort.default(decreasing=T) %>>% names() %>>% "["(1:topn)
  a1=prune_taxa(taxa_names(a0)%in%b0,a0)
  #    subset_taxa(a0,taxa_names(a0)%in%b0)

  rankv=rank_names(phylo)[1:which(rank_names(phylo)==ranklv)]
  a2=psmelt(a1) %>>% mutate_at(.vars=rankv,funs(as.character))

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

  if(!missing(anno_var)){
    ind=anno_var%in%sample_variables(phylo)
    if(any(ind)){
      ind2=anno_var[ind]
      a4=get_variable(a1,ind2) %>% data.frame(stringsAsFactors=F) %>% setNames(ind2)
      rownames(a4)=sample_names(a1)
      pheatmap::pheatmap(a3,annotation_col=a4,...)
      return(list(RA=a3,anno_tax=a5,anno_sample=a4))
    }
  }
  pheatmap::pheatmap(a3,...)
  return(list(RA=a3,anno_tax=a5))
}
