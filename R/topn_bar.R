#
#
#
#
#
#
#

topn_bar=function(u0,ranklv,y="fold",x="bn",g="col",topn=10,colorv=c("firebrick1","royalblue"),...){
  if(!missing(ranklv)){
    u1=u0[u0[["Rank"]]==ranklv,]
  } else {
    u1=u0
  }
  u1$val=abs(u1[[y]])%>>%jitter(1e-6)
  u2=u1 %>>% group_by_(g) %>>% top_n(topn,val)
  u3=data.frame(val=rep(0,2*topn-nrow(u2)))
  colnames(u3)=y
  u4=bind_rows(u2,u3) %>>% dplyr::select_(.dots=c(g,x,y)) %>>% arrange_(.dots=sprintf("desc(%s)",y))
  lab=u4[[g]][!duplicated(u4[[g]])]
  ypre=pretty(u4[[y]])
  #  if(head(ypre,1)==0) ypre=c(-max(ypre),ypre)
  #  if(tail(ypre,1)==0) ypre=c(ypre,-min(ypre))
  if(any(is.infinite(u4[[y]]))) u4[[y]] %<>% plyr::mapvalues(from=c(-Inf,Inf),to=range(ypre))

  par(las=1,mar=c(2,5,2,1),xpd=NA,srt=0)
  plot(c(0,topn*2),range(ypre),type='n',bty='n',xaxt='n',yaxt='n',...)
  axis(2,at=ypre,labels=abs(ypre))
  barplot(u4[[y]],width=0.8,space=c(0.75,rep(0.25,topn*2-1)),add = T,yaxt='n',col=rep(colorv,each=topn),border = NA)
  legend("bottom",lab,pch=15,col =colorv,horiz = T)
  par(srt=90)
  text(1:topn,-diff(range(ypre))/20,u4[[x]][1:topn],adj=c(1,0.5))
  text(topn+(1:topn),diff(range(ypre))/20,u4[[x]][topn+(1:topn)],adj=c(0,0.5))
}
