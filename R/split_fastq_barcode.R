#
#
#依barcode將較大FASTQ檔分割成數個小檔
#barcode主要是長度8的字串
#一般是R1(forward), R2(reverse)兩個檔
#因為R1,R2都可能會出現barcode1,barcode2,
#所以都要判斷有無出現barcode
#
#

#parameter

#R1="E:/Biom/MS16022/20160412/test/MS16022-1_S1_L001_R1_001.clean.fastq"
#R2="E:/Biom/MS16022/20160412/test/MS16022-1_S1_L001_R2_001.clean.fastq"

#barcode=read_tsv("E:/Biom/MS16022/barcode.txt") %>%
#  filter(libName=="MS16022-1") %>%
#  select(sampleName,barcode1,barcode2)

#length.barcode=8

#max.mismatches=1

#savepath="E:/Biom/MS16022/20160412/test/a"

#main

#get r1, r2

split_fastq=function(R1,R2,barcode,savepath,max.mismatches=1,length.barcode=8){
  r1.fastq=readFastq(R1)
  r2.fastq=readFastq(R2)

  y2=data.frame(R1=substr(r1.fastq@sread,1,length.barcode)
                ,R2=substr(r2.fastq@sread,1,length.barcode)
                ,stringsAsFactors = F) %>%
    add_rownames()

  string_match=function(x,set){
    x0=substring(x,1:length.barcode,1:length.barcode)
    x1=sapply(set, function(y){
      sum(x0!=y)
    }) %>% "["(order(.)[1])
    return(data.frame(Query=x,Subject=names(x1),mismatch=unname(x1),stringsAsFactors = F))
  }

  a0=barcode %$% c(barcode1,barcode2) %>% unique()
  a1= a0 %>% lapply(substring,1:nchar(a0[1]),1:nchar(a0[1]))
  names(a1)=a0

  y3=unique(y2$R1) %>% lapply(string_match,set=a1) %>% bind_rows() %>%
    filter(mismatch<=max.mismatches)

  y4=unique(y2$R2) %>% lapply(string_match,set=a1) %>% bind_rows() %>%
    filter(mismatch<=max.mismatches)

  y5=y2%>%left_join(y3,by=c("R1"="Query"))%>%
    left_join(y4,by=c("R2"="Query"))%>%
    filter(!is.na(Subject.x)&!is.na(Subject.y))%>%
    left_join(barcode,by=c("Subject.x"="barcode1","Subject.y"="barcode2"))%>%
    left_join(barcode,by=c("Subject.x"="barcode2","Subject.y"="barcode1"))%>%
    mutate(sampleName=ifelse(is.na(sampleName.x),sampleName.y,sampleName.x))

  y7=y5 %$% split(rowname,sampleName) %>% lapply(as.integer)

  #split and export to small fastq

  #i=1
  for(i in 1:length(y7)){
    y9=y7[[i]] %>% "["(r1.fastq,.)
    y10=DNAStringSet(y9@sread,length.barcode+1)
    names(y10)=sprintf("%s___%08i",names(y7)[i],1:length(y10))#substr(y9@id,1,max(width(y9@id)))

    writeXStringSet(y10, sprintf("%s/%s.R1.fastq",savepath,names(y7)[i]), format="fastq",
                    qualities=BStringSet(y9@quality@quality,length.barcode+1))

    y9=y7[[i]] %>% "["(r2.fastq,.)
    y10=DNAStringSet(y9@sread,length.barcode+1)
    names(y10)=sprintf("%s___%08i",names(y7)[i],1:length(y10))#substr(y9@id,1,max(width(y9@id)))

    writeXStringSet(y10, sprintf("%s/%s.R2.fastq",savepath,names(y7)[i]), format="fastq",
                    qualities=BStringSet(y9@quality@quality,length.barcode+1))
  }
  cat("\nNumber of sequences in each fastq file:\n\n")
  sapply(y7,length)%>%print()
  gc()
}




