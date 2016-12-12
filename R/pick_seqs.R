#
#
#
#
#
#

pick_seqs=function(ggid,rdspath,pickfrom="cluster"){
  a3=read_rds(rdspath)
  a4=dplyr::filter(a3$annot,X3%in%ggid)
  if(pickfrom=="cluster"){
    a5=a3$clusterfasta[a4$X1]
  } else {
    a5=a4 %>>% (dplyr::filter(a3$map,X10%in%.$X2)) %>>% ("["(a3$nchfasta,.$X9))
  }
  return(a5)
}
