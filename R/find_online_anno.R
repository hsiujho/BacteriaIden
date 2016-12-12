#
#After Blast, find correponding annotation of gi online
#
#
#
#
#

find_online_anno=function(gi,NumBatch=20){
  b3=gi %>% ntile(length(.)%/%NumBatch) %>% split(gi,.) %>%
    lapply(function(x){
      repeat{
        c0=try(entrez_summary(db="nuccore",x))
        if(!inherits(c0,"try-error")) break
      }
      return(c0)
    })

  k0=c("gi","title","taxid","slen","organism","strain")

  b4=lapply(b3,function(x){
    lapply(x,function(y){
      extract_from_esummary(y,intersect(names(y),k0))
    })
  }) %>% melt() %>% dcast(L2~L3,value.var="value") %>% mutate(L2=NULL)

  taxid=b4[,"taxid"] %>% unique()

  b5=ntile(taxid,length(taxid)%/%NumBatch) %>% split(taxid,.) %>%
    lapply(function(x){
      repeat{
        tax_rec <- try(entrez_fetch(db="taxonomy", id=x, rettype="xml", parsed=TRUE))
        if(!inherits(tax_rec,"try-error")) break
      }
      tax_list <- XML::xmlToList(tax_rec)
      c0=lapply(tax_list,function(x){
        names(x$LineageEx)=paste0(names(x$LineageEx),1:length(x$LineageEx))
        x1=x$LineageEx%>%melt%>%dcast(L1~L2,value.var="value")%>%filter(Rank!="no rank")
        x1$taxid=x$TaxId
        x1$L1=NULL
        return(x1)
      }) %>% do.call(rbind,.) %>% dcast(taxid~Rank,value.var="ScientificName")
    }) %>% bind_rows()

  b6=select(b5,taxid,superkingdom,phylum,class,order,family,genus,species) %>%
    left_join(b4,.,by="taxid") %>%
    select(superkingdom,phylum,class,order,family,genus,species,organism,strain,gi,title,taxid,slen)
  return(b6)
}
