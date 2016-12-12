
ClustIden=function(ID,R1,R2,gz=T,minlength=450,maxlength=550
                   ,wdpath="E:/fastq",savepath,cores=detectCores()-1){

  setwd(wdpath)

  #解壓縮檔案至工作資料夾

  if(gz==T){
    gunzip(filename=R1,destname="a.R1.fastq",remove=F, overwrite=T)
    gunzip(filename=R2,destname="a.R2.fastq",remove=F, overwrite=T)
  } else {
    file.copy(from=R1,to="a.R1.fastq")
    file.copy(from=R2,to="a.R2.fastq")
  }

  r1=readFastq("a.R1.fastq") %>% length()
  r2=readFastq("a.R2.fastq") %>% length()

  #合併R1和R2

  system("usearch -fastq_mergepairs a.r1.fastq -reverse a.r2.fastq -fastq_minovlen 8 -fastqout a.merged.fastq")

  #分割序列及品質分數

  system('mothur "#fastq.info(fastq=a.merged.fastq)"')

  #merge後的序列摘要

  summary.merge= system('mothur "#summary.seqs(fasta=a.merged.fasta)"',intern = T) %>%
    mothur.summary()

  #依研究及品質門檻篩選序列

  system(sprintf('mothur "#trim.seqs(fasta=a.merged.fasta,qfile=a.merged.qual,qaverage=27,minlength=%i,maxlength=%i,maxambig=0,maxhomop=8,processors=1)"',minlength,maxlength))

  #trim後的序列摘要

  summary.trim= system('mothur "#summary.seqs(fasta=a.merged.trim.fasta)"',intern = T) %>%
    mothur.summary()

  #分割檔案, 避開記憶體的限制

  system('perl fasta-splitter.pl --n-parts 2 A.merged.trim.fasta')

  #檢測chimera

  system(sprintf('usearch -uchime_ref A.merged.trim.part-1.fasta -db rdp_gold.fa -nonchimeras A.merged.trim.part-1.fasta.nch -strand plus -mindiv 3 -threads %i',cores))
  system(sprintf('usearch -uchime_ref A.merged.trim.part-2.fasta -db rdp_gold.fa -nonchimeras A.merged.trim.part-2.fasta.nch -strand plus -mindiv 3 -threads %i',cores))

  #合併檢測後的序列

  shell('type A.merged.trim.part-1.fasta.nch A.merged.trim.part-2.fasta.nch > A.merged.trim.fas.nch')

  #non-chimera後的序列摘要

  summary.nchi= system('mothur "#summary.seqs(fasta=a.merged.trim.fas.nch)"',intern = T) %>%
    mothur.summary()

  #各篩選階段有幾個序列數

  b4=data.frame(
    Step=c("R1","R2","Merged","trim","nchi")
    ,No.Seqs=c(r1,r2,c(read_lines("a.merged.summary") %>% length()
               ,read_lines("a.merged.trim.summary") %>% length()
               ,read_lines("a.merged.trim.fas.summary") %>% length()
    )-1)
    ,sampleID=ID
  )

  #clustering後鑑定
  #  刪除重覆序列

  system('usearch -derep_fulllength a.merged.trim.fas.nch -fastaout A.derep.fas -sizeout')

  #  依序列重覆數及名稱排序

  system('usearch -sortbysize A.derep.fas -fastaout A.sorted.fas')

  #  序列分群並取得代表群集的序列

  system('usearch -cluster_otus A.sorted.fas -otus A.otus.fas')

  #  鑑定群集的序列, 以97_otus當作比對檔

  system(sprintf('usearch -usearch_global A.otus.fas -db gg_13_5_otus/97_otus.udb -userout A.usearch.out -id 0.5 -userfields query+target+id+qcov -strand both -top_hit_only -threads %i',cores))

  #  串聯菌種編號的註解檔

  system('perl pick.usearch.otus.pl -s A')

  #  利用鑑定後的群集序列當作比對檔, 比對non-chimera後的全部序列

  system(sprintf('usearch -usearch_global a.merged.trim.fas.nch -db a.pick.labeled.otus -strand plus -top_hit_only -id 0.97 -uc A.map.uc -threads %i',cores))

  #  匯入鑑定結果

  if(file.info("A.map.uc")[["size"]]!=0)
  {
    a4=read_tsv("A.map.uc",col_names=F) %>%
      group_by(X10) %>%
      summarise(num_seq=length(X10)) %>% ungroup() %>%
      left_join(read_tsv("A.pick.gg.taxonomy",col_names=F),by=c("X10"="X2")) %>%
      dplyr::rename(OTUID=X3) %>% group_by(OTUID) %>%
      summarize(N=sum(num_seq)) %>% mutate(sampleID=ID)
  }

  #保留檔案轉存

  a5=list(
    nchfasta=readAAStringSet("A.merged.trim.fas.nch")
    ,norepfasta=readAAStringSet("A.sorted.fas")
    ,clusterfasta=readAAStringSet("A.otus.fas")
    ,annot=read_tsv("A.pick.gg.taxonomy",col_names=F)
    ,map=read_tsv("A.map.uc",col_names=F)
    ,summary=list(merge=summary.merge
                  ,trim=summary.trim
                  ,nchi=summary.nchi)
  )
  saveRDS(a5,file = file.path(savepath,sprintf("%s.rds",ID)))

  #刪除檔案

  list.files(pattern =".logfile") %>% file.remove()
  list.files(pattern ="^[Aa]\\.") %>% file.remove()

  return(list(pipe=b4,otus=a4))
}
