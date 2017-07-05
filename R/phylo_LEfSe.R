#
#
#
#
#
#

phylo_LEfSe=function(
  phylo=a2  #phyloseq的物件
  ,g1="ph_lv" #分組變數名稱, 物件中的sample_data有的欄位
  ,save_path="E:\\Biom\\20150901 biom_analy\\20151019 GF SA\\" #輸出存檔路徑, LEfSe_"g1".png, LEfSe_"g1".rds
  ,lefse_path="E:\\soft\\python\\LEfSe\\lefse-db64b6287cd6\\home\\ubuntu\\lefse_to_export\\"
#  ,figure_format="png" #png,pdf,svg
  ,labeled_stop_lev=7
  ,abrv_start_lev=5
  ,abrv_stop_lev=7
){
  L1=otu_table(phylo)%>%data.frame(check.names = F,stringsAsFactors=F)%>%tibble::rownames_to_column()
  L2=tax_table(phylo)%>%data.frame(check.names = F,stringsAsFactors=F)%>%tibble::rownames_to_column()
  L2$X16S_rRNA_Count=NULL
  L3=left_join(L1,L2,by="rowname")

  rank_lv=rank_names(phylo)[1:7]
  #x=L3[83,]
  L4=apply(L3,1,function(x){
    x1=lapply(x,rep,each=7)%>%data.frame(check.names = F,stringsAsFactors=F)
    x1$class=sapply(1:7,function(i){
      if(is.na(x1[i,rank_lv[i]])|x1[i,rank_lv[i]]==""){
        return("")
      } else {
        paste(x1[i,rank_lv][1:i],collapse="|",sep="")
      }
    })
    x2=dplyr::filter(x1,class!="")
    x2$rowname=NULL
    for(i in 1:7){x2[[rank_lv[i]]]=NULL}
    return(x2)
  }) %>% bind_rows()

  # L5=mutate_each(L4,funs(as.numeric),-class)%>%group_by(class)%>%summarise_each(funs(sum))
  L5=group_by(L4,class)%>%mutate_all(funs(as.numeric))%>%summarise_all(funs(sum))
  L6=rbind(c("group",as.character(get_variable(phylo,g1))),as.matrix(L5))

  figure_format=c("png","pdf","svg")
  prodClado=function(x6){
    write.table(x6,file=sprintf("%sz.txt",lefse_path), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F, qmethod = c("escape", "double"),fileEncoding = "")

    owd=getwd()
    setwd(lefse_path)
    #以下為LEfSe分析步驟, 參數設定請看各程式碼的檔案
    system("python format_input.py z.txt z.in -c 1 -s -1 -u -1 -o 1000000",intern =T)
    system("python run_lefse.py z.in z.res -y 0",intern =T)
    #-y 1篩選條件嚴格, -y 0不嚴格
    para=sprintf("%s %s %s %s"
                 , sprintf("--labeled_stop_lev %i",labeled_stop_lev) #標籤顯示到第六層級
                 , sprintf("--abrv_start_lev %i --abrv_stop_lev %i",abrv_start_lev,abrv_stop_lev) #以代碼顯示標籤的起迄層級
                 , "--expand_void_lev 2"
                 #在較高層級就中斷的路徑是(1)否(2)延伸到最低層級
                 #, '--title "Cladogram"' #Cladogram的標題
                 , sprintf("--title \"\"")
    )

    for(i in 1:3){
      tx1=sprintf("python plot_cladogram.py z.res z.%s --dpi 300 --format %s",figure_format[i],figure_format[i])
      system(sprintf('%s %s',tx1,para), intern =T )
    }
    # BarChart of LDA score (log10)
    system("python plot_res.py z.res BarChart.png --dpi 300",intern =T)
    system("python plot_res.py z.res BarChart.pdf --format pdf",intern =T)
    system("python plot_res.py z.res BarChart.svg --format svg",intern =T)

    #LEfSe輸出存放路徑
    outp=read.delim(sprintf("%sz.res",lefse_path), header =F,stringsAsFactors=F)
    setwd(owd)
    return(outp)
  }

  L7=prodClado(L6) %>% mutate(pvalue=as.numeric(ifelse(V5=="-","",V5)),lv=nchar(as.vector(V1))-nchar(gsub("\\.","",as.vector(V1))),Rank=rank_lv[lv+1]) %>% arrange(lv,V3,pvalue)%>%dplyr::rename(Taxonomy=V1,Group=V3,`LDA SCORE (log 10)`=V4)

  for(i in 1:3){
    cladofn=sprintf("%sz.%s",lefse_path,figure_format[i])
    if(file.exists(cladofn)) {
  #    file.copy(from=cladofn,to=sprintf("%sLEfSe_%s.png",save_path,g1),overwrite=T)
      file.copy(from=cladofn,to=file.path(save_path,sprintf("LEfSe_%s.%s",g1,figure_format[i])),overwrite=T)
      file.remove(sprintf("%sz.%s",lefse_path,figure_format[i]))
    }
    barchartfn=sprintf("%sBarChart.%s",lefse_path,figure_format[i])
    if(file.exists(barchartfn)) {
      file.copy(from=barchartfn,to=file.path(save_path,sprintf("LEfSe_BarChart_%s.%s",g1,figure_format[i])),overwrite=T)
      file.remove(sprintf("%sBarChart.%s",lefse_path,figure_format[i]))
    }
  }

  file.remove(sprintf("%sz.txt",lefse_path))
  file.remove(sprintf("%sz.in",lefse_path))
  file.remove(sprintf("%sz.res",lefse_path))

  L8=melt(L5,id.vars="class", measure.vars=colnames(L5)[colnames(L5)!="class"],variable.name = "ID",value.name = "value") %>% mutate_at(funs(as.character),.vars="ID") %>% left_join(data.frame(ID=sample_names(phylo),g1=get_variable(phylo,g1),stringsAsFactors = F,check.names = F),by="ID") %>% group_by(class,g1) %>% summarise(mean=mean(value)) %>% dcast(class~g1,value.var ="mean")

  L9=melt(L8, id.vars="class",variable.name = "Group",value.name = "mean")%>%mutate(logmean=ifelse(mean==0,log2(.001),log2(mean)))%>%group_by(class)%>%summarise(meandiff=max(mean)-min(mean),log2meandiff=max(logmean)-min(logmean))%>%right_join(L8,by="class") %>% mutate(Taxonomy=gsub("\\|","\\.",class)%>%gsub("\\[","_",.)%>%gsub("\\]","_",.)%>%gsub("-","_",.)%>%gsub(" ","",.),meanlog10diff=log10(meandiff))

  #LEfSe中class有些特殊字元會被改成底線或其他字元,所以要把L8$class轉換成一樣才能合併,下面是轉換的指令
  #L8$class1=gsub("\\|","\\.",L8$class)%>%gsub("\\[","_",.)%>%gsub("\\]","_",.)%>%gsub("-","_",.)%>%gsub(" ","",.)
  #辨認兩個名稱是否一致
  #all.equal(sort(L8$class1),sort(L7$V1))
  #查詢不一樣的名稱
  #L7$V1[!L7$V1%in%L8$class1]
  #L8$class[!L8$class1%in%L7$V1]

  saveRDS(list(group.means=L8,meandiff=L9,LEfSe=L7)
#          ,file=sprintf("%sLEfSe_%s.rds",save_path,g1)
           ,file=file.path(save_path,sprintf("LEfSe_%s.rds",g1))
          )
  return()
}
