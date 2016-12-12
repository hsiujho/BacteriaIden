#
#
#
#
#
#
#

donuts_plot=function(df,y.var,group.var,fill.var,facet.row=".",facet.col=".",x_shift=3){
  d1=df %>% arrange_(group.var,fill.var)
  d1$RA=d1[[y.var]]

  d2=split(d1,d1[,c(facet.row,facet.col,group.var) %>% "["(.!=".")]) %>%
    lapply(function(x){
      mutate(x,ymax=cumsum(RA),ymin=c(0,head(ymax,n=-1)))
    }) %>% bind_rows()

  gv=d2[[group.var]] %>% factor()
  d2$xmax=gv %>% as.numeric()+x_shift+.4
  d2$xmin=gv %>% as.numeric()+x_shift-.4
  d2$fill_var=d2[[fill.var]]

  p1=ggplot(d2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_rect(aes(fill = fill_var), colour = "grey50")+
    coord_polar(theta="y")+ylab("Relative Abundance")+
    scale_x_continuous(name=group.var,breaks=1:nlevels(gv)+x_shift
                       ,labels=levels(gv),limits=c(0,nlevels(gv)+x_shift+.5))+
    scale_fill_discrete(name=fill.var)
  if(facet.row!="."|facet.col!=".")
    p1=p1+facet_grid(formula(sprintf("%s~%s",facet.row,facet.col)))
  return(p1)
}
