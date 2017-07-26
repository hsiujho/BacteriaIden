#
#
#
#
#
#
#

phylo_compare_2_groups=function(phylo,group.var){
  if(class(phylo)!="phyloseq") stop("The class of phylo is not phyloseq")
  if(length(group.var)!=1) stop("length of group.var must be 1")
  if(all(sample_variables(phylo)!=group.var)) stop("phylo not contain group variable")
  gv=get_variable(phylo,group.var) %>>% as.character() %>>% factor()
  sample_data(phylo)$gv=gv
  #  print(gv)
  #  print(nsamples(phylo))
  #  cat(sprintf("length(levels(gv))=%i\n",length(levels(gv))))
  if(length(levels(gv))!=2) stop("just for 2 groups")

  e0 = phylo %>>% psmelt()

  e1=e0 %>>% group_by(OTU) %>>% do(k0=anova(lm(Abundance~gv,data=.))$`Pr(>F)`[1]) %>>% summarise(OTU,pvalue=k0[1]) %>>% setNames(c("ID","pvalue_t.test"))

  e2=e0 %>>% group_by(OTU,gv) %>>% summarise(Mean=mean(Abundance)) %>>% dcast(OTU~gv,value.var="Mean") %>>% setNames(c("ID",paste0(colnames(.)[-1],"_mean")))

  e3=e0 %>>% group_by(OTU) %>>% do(k0=wilcox.test(Abundance~gv,data=.,exact = FALSE)) %>>% summarise(OTU,pvalue=k0$p.value) %>>% setNames(c("ID","pvalue_wilcox.test"))

  e4=e0 %>>% group_by(OTU,gv) %>>% summarise(Median=median(Abundance)) %>>% dcast(OTU~gv,value.var="Median") %>>% setNames(c("ID",paste0(colnames(.)[-1],"_median")))

  e5=left_join(e2,e1,by="ID") %>>% left_join(e4,by="ID") %>>% left_join(e3,by="ID")

  dots0=sprintf("%s_mean",levels(gv))
  dots1=paste0(dots0,collapse=">")
  dots=list(Log2Fold=lazyeval::interp(~f(a)-f(b),f=quote(log2),a=as.name(dots0[1]),b=as.name(dots0[2])),Diff=lazyeval::interp(~f(a,b),f=as.name("-"),a=as.name(dots0[1]),b=as.name(dots0[2])),HighGroup=lazyeval::interp(~f1(f2(c,d),a,b),f1=quote(ifelse),f2=as.name(">"),c=as.name(dots0[1]),d=as.name(dots0[2]),a=levels(gv)[1],b=levels(gv)[2]))

  e6=mutate_(e5,.dots=dots)%>>%mutate(Significant=ifelse(pvalue_t.test>=0.05|abs(Log2Fold)<1,"Non",HighGroup)) %>>% arrange(Log2Fold,pvalue_t.test) %>>% left_join(tax_table(phylo)@.Data %>>% data.frame(stringsAsFactors=F) %>>% tibble::rownames_to_column("ID"),by="ID")
  #,Link=sprintf('=HYPERLINK("http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=%s","%s")',ID,ID)
  outp=list(raw=e0,compare=e6)

  if(any(regexpr("Description",colnames(e6))!=-1)){
    ff0=lapply(levels(gv),function(x){
#      x="C02"
      f0=dplyr::filter(e6,Significant==x) %>>% dplyr::select(ID,metadata_COG_Category,metadata_COG_Description)
      if(nrow(f0)==0) return()
      f1=gsub("[[:punct:]]"," ",f0$metadata_COG_Description) %>>% strsplit(" ") %>>% unlist() %>>% (gsub(",","",.)) %>>% table() %>>% data.frame() %>>% arrange(desc(Freq)) %>>% setNames(c("string",x)) %>>% mutate_each(funs(as.character),string) %>>% dplyr::filter(string!="")
      return(f1)
    })
    if(!is.null(ff0[[1]])&is.null(ff0[[2]])) {
      outp$textMining=ff0[[1]]
    }
    if(is.null(ff0[[1]])&!is.null(ff0[[2]])) {
      outp$textMining=ff0[[2]]
    }
    if(!is.null(ff0[[1]])&!is.null(ff0[[2]])) {
      outp$textMining=full_join(ff0[[1]],ff0[[2]],by="string")
    }
  }
  return(outp)
}

function(){
#Example 1: Mcirobiom analysis
  b1=read_rds("E:/Biom/MS16036/greengene 20160531/MS16036_gg_phylo.rds") %>>% subset_samples(ST=="After"&group%in%c("H","NE")) %>>% tax_glom("Class",NArm=F) %>>% {prune_taxa(taxa_sums(.)>0,.)} %>>% transform_sample_counts(function(x) x/sum(x)*100) %>>% phylo_compare_2_groups(group.var="group")
  ggplot(b1$raw,aes(x=group,y=Abundance))+geom_boxplot()+facet_wrap(~Class,scales = "free_y")
  topn_bar(b1$compare,y="Log2Fold",x="Class",g="HighGroup",lab=c("H","NE"),topn=10,colorv=c("firebrick1","royalblue"),ylab="Log2Fold")
#Example 2: Functional prediction
  b0=read_rds("E:/Biom/MS16036/R code/20160616/fp_cog.rds")$fp %>>% subset_samples(ST=="After"&group%in%c("H","NE")) %>>% (prune_taxa(taxa_sums(.)>0,.)) %>>% transform_sample_counts(function(x)x/sum(x)*10000) %>>% phylo_compare_2_groups(group.var="group")
  ggplot(b0$compare,aes(x=Log2Fold,y=-log10(pvalue_t.test),col=Significant,size=abs(Diff)))+geom_point(alpha=0.5)
  b01=dplyr::filter(b0$compare,Significant!="Non"&abs(Diff)>0.001)
  topn_bar(b01,y="Log2Fold",x="ID",g="HighGroup",lab=c("H","NE"),topn=10,colorv=c("firebrick1","royalblue"),ylab="Log2Fold")
  topn_bar(b01,y="Diff",x="ID",g="HighGroup",lab=c("H","NE"),topn=10,colorv=c("firebrick1","royalblue"),ylab="Diff",ylim=c(-10,20))
  b02=dplyr::filter(b0$raw,OTU%in%top_n(b01,12,-Log2Fold)$ID)
  ggplot(b02,aes(x=group,y=Abundance))+geom_boxplot()+facet_wrap(~OTU,scales = "free_y")
  b02=dplyr::filter(b0$raw,OTU%in%top_n(b01,12,-Diff)$ID)
  ggplot(b02,aes(x=group,y=Abundance))+geom_boxplot()+facet_wrap(~OTU,scales = "free_y")
}

