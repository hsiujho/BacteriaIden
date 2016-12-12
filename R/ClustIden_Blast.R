#BLAST
#
#會使用到ClustIden的輸出
#
#
#

ClustIden_Blast=function(
  FileFrom
  ,FileName
  ,SaveTo
  ,dbpath="E:/Biom/NCBI_DB/refseq_rna/refseq_rna"
  ,NumClust=1
  ,Npara=1
){
  if(NumClust*Npara>detectCores()) stop("NumClust*Npara more than the number of threads.")
  on=Sys.time()
  cl=makeCluster(NumClust)
  parallel::clusterExport(cl,c("FileFrom","SaveTo","Npara","dbpath"),envir=environment())
  parLapply(cl,FileName,function(x){
    require(rBLAST)
    require(magrittr)
    bl <- blast(db=dbpath)
    readRDS(sprintf("%s/%s.rds",FileFrom,x)) %$%
      predict(bl, norepfasta, BLAST_args=sprintf("-max_target_seqs 1 -perc_identity 85 -num_threads %i",Npara)
              , custom_format="qseqid qlen sseqid sgi slen pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore score ppos staxids") %>%
      saveRDS(sprintf("%s/%s.rds",SaveTo,x))
  })
  stopCluster(cl)
  off=Sys.time()-on
  cat(sprintf("\nIt tooks %.2f %s\n\n",off,attr(off,"units")))
}


