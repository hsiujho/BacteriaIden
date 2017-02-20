
phylo_ISA=function(phylo,group_var
                   ,ranklv="Phylum"
                   ,NArm=F,toRA=T)
{
  if(!missing(group_var)){
    if(any(sample_variables(phylo)==group_var)){
      gv=get_variable(phylo,group_var)
      if(!is.factor(gv)) gv%<>%factor()
    } else {
      stop("Variable not found.")
    }
  } else {
    stop("group_var is missing.")
  }

  if(any(rank_names(phylo)==ranklv)){
    phylo%<>%tax_glom(ranklv,NArm=NArm)%>>%(prune_taxa(taxa_sums(.)>0,.))
  } else {
    stop(sprintf("%s not found",ranklv))
  }

  if(toRA){
    phylo%<>%transform_sample_counts(function(x)x/sum(x)*100)
  }

  bp1=otu_table(phylo)@.Data %>>% t() %>>%data.frame(check.names=F)
  bp2=indicspecies::multipatt(bp1, gv, control = how(nperm=9999))
  bp3=bp2$sign %>>% rownames_to_column("OTU") %>>%
    filter(p.value<0.05) %>>% arrange(index,p.value) %>>%
    left_join(tax_table(phylo)@.Data %>>%
                as.data.frame(stringsAsFactors=F) %>>%
                rownames_to_column("OTU"),by="OTU")

  bp4=psmelt(phylo) %>>%
    mutate_(tmp=ranklv) %>>%
    filter(OTU%in%bp3$OTU&tmp!="") %>>%
    mutate(OTU=factor(OTU,levels=bp3$OTU)) %>>%
    arrange(OTU) %>>%
    mutate(tmp=factor(tmp,levels=tmp[!duplicated(tmp)]))

  bp5=ggplot(bp4,aes_string(x=group_var,y="Abundance"))+
    geom_boxplot()+
    facet_wrap(~tmp,scales = "free_y",ncol=4)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,v=0.5))+
    ylab(ifelse(toRA,"Relative abundance (%)","Abundance"))
  return(list(ISA=bp2,signif=bp3,signif_box=bp5))
}

# cc0=phylo_ISA(phylo,group_var="STGroup",ranklv="Genus")
# cc1=phylo_ISA(phylo,group_var="STGroup",ranklv="Class")
