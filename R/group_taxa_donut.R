#
#
#
#
#
#
#

group_taxa_donut=function(phylo,group.varn,ranklv="Phylum",topn=5,facet.row=".",facet.col=".",x_shift=3){

  if(all(sample_variables(phylo)!=group.varn)|is.na(group.varn)){
    sample_data(phylo)$dummy_var=""
    group.varn="dummy_var"
  }

  if(all(sample_variables(phylo)!=facet.row)){
    facet.row="."
  }

  if(all(sample_variables(phylo)!=facet.col)){
    facet.col="."
  }

  a2=phylo %>%
    tax_glom(ranklv,NArm=T) %>%
    transform_sample_counts(function(x){x/sum(x)*100})

  a3=taxa_sums(a2) %>% order(decreasing=T) %>% "["(-(1:topn)) %>%
    merge_taxa(a2,.,2) %>% psmelt() %>%
    mutate_(rank_var=ranklv,x_var=group.varn) %>%
    mutate_each(funs(as.character),rank_var)
  a3$rank_var%<>%ifelse(is.na(.),"Others",.)

  a3$rank_var%<>%factor(levels=group_by(a3,rank_var)%>%
                          summarise(RA=sum(Abundance))%>%arrange(desc(RA))%$%rank_var)

  a4=sprintf("group_by(a3,%s%sx_var,rank_var)"
             ,ifelse(facet.row==".","",paste0(facet.row,","))
             ,ifelse(facet.col==".","",paste0(facet.col,",")))
  df=eval(parse(text=a4))%>%summarise(RA=mean(Abundance))

  d1=df

  d2=split(d1,d1[,c(facet.row,facet.col,"x_var") %>% "["(.!=".")]) %>%
    lapply(function(x){
      mutate(x,ymax=cumsum(RA),ymin=c(0,head(ymax,n=-1)))
    }) %>% bind_rows()

  gv=d2$x_var %>% factor()
  d2$xmax=gv %>% as.numeric()+x_shift+.4
  d2$xmin=gv %>% as.numeric()+x_shift-.4

  p1=ggplot(d2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_rect(aes(fill = rank_var), colour = "grey50")+
    coord_polar(theta="y")+ylab("Relative Abundance")+
    scale_x_continuous(name=ifelse(group.varn=="dummy_var","",group.varn)
                       ,breaks=1:nlevels(gv)+x_shift
                       ,labels=levels(gv),limits=c(0,nlevels(gv)+x_shift+.5))+
    scale_fill_discrete(name=ranklv)+
    theme(legend.background=element_blank(),
          plot.background = element_blank())
  if(facet.row!="."|facet.col!=".")
    p1=p1+facet_grid(formula(sprintf("%s~%s",facet.row,facet.col)))

  return(p1)
}
