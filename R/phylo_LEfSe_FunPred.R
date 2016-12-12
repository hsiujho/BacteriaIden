#
#
#
#
#
#

phylo_LEfSe_FunPred=function(phylo,gv.name,save_path){
  L1 = phylo %>>% (prune_samples(get_variable(.,gv.name)!="",.)) %>>% (prune_taxa(taxa_sums(.)>0,.))

  L2=data.frame(otu_table(L1)@.Data,check.names = F) %>>% add_rownames("group")

  L6=rbind(c("group",as.character(get_variable(L1,gv.name))),L2)
  lefse_path="E:/soft/python/LEfSe/lefse-a31c10fe09c8/home/ubuntu/lefse_to_export/"
  write.table(L6,file=sprintf("%sz.txt",lefse_path), append = FALSE, quote = F, sep = "\t"
              ,eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F
              , qmethod = c("escape", "double"),fileEncoding = "")

  owd=getwd()
  setwd(lefse_path)
  #以下為LEfSe分析步驟, 參數設定請看各程式碼的檔案
  system("python format_input.py z.txt z.in -c 1 -s -1 -u -1 -o 1000000",intern =T)
  system("python run_lefse.py z.in z.res -y 0",intern =T)
  #LEfSe輸出存放路徑
  outp=read.delim(sprintf("%sz.res",lefse_path), header =F,stringsAsFactors=F) %>% mutate(pvalue=as.numeric(ifelse(V5=="-","",V5))) %>% arrange(desc(V3),pvalue) %>% dplyr::rename(ID=V1,Group=V3,abs_LDA_score=V4)
  setwd(owd)
  saveRDS(outp,file=sprintf("%sLEfSe_%s.rds",save_path,gv.name))
  return()
}
