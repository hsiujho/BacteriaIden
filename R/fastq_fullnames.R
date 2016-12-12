fastq_fullnames=function(fastq_path){
  #get the fullnames of fastq files
  fn2=list.files(fastq_path,pattern=".fastq",full.names = T)
  data.frame(fn=fn2,stringsAsFactors = F) %>%
    mutate(ID=fn2 %>% sub(".*/","",.) %>% sub("\\..*","",.)
           ,R12=strsplit(fn2,"\\.") %>% sapply(intersect,y=c("R1","R2"))
    ) %>% acast(ID~R12,value.var="fn") %>% return()
}


check_soft=function(){
  return(list(mothur_ver=system('mothur "#quit()"',intern = T) %>% "["(nchar(.)!=0)
              ,usearch_ver=system('usearch',intern = T) %>% "["(nchar(.)!=0)))
}


mothur.summary=function(x){
  x %>%
    "["(regexpr("\t",.)!=-1) %>%
    strsplit("\t") %>% lapply(function(x){
      if(length(x)!=8){
        return(c("",x,rep("",7-length(x))))
      } else {
        return(x)
      }
    }) %>% do.call(rbind,.)
}
