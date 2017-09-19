#
#
#
#
#
#

add_label_to_unknown_taxon=function(phylo,unknown_symbol=c("",NA)){

  # 把未知的分類補上該層級上一層及該層級字首小寫字母加底線, 語法例句如下
  # t0$Genus=ifelse(t0$Genus=="",sprintf("%s;g_",t0$Family),t0$Genus)
  #確認輸入的物件分類層級

  u0=which(apply(tax_table(phylo),2,function(x)sum(x%in%unknown_symbol)) != ntaxa(phylo)) %>>% max()
  if(u0==1) return(phylo)
  v0=sprintf("%%s;%s_",c("p","c","o","f","g","s"))
  for(i in 2:u0){
    phylo@tax_table@.Data[,i]=ifelse(
      phylo@tax_table@.Data[,i]%in%unknown_symbol
      ,sprintf(v0[i-1],phylo@tax_table@.Data[,i-1])
      ,phylo@tax_table@.Data[,i])
  }
  return(phylo)
}
