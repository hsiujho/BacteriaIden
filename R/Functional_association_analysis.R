#
#
#
#
#
#

Functional_association_analysis=function(phylo,y,Is_y_discrete=F){
  if(missing(phylo)) stop("phylo is missing")
  if(missing(y)) stop("y is missing")
  if(missing(Is_y_discrete)) stop("Is_y_discrete is missing")
  if(!is.logical(Is_y_discrete)) stop("Is_y_discrete should be logical (TRUE or FALSE)")
  if(Is_y_discrete & length(table(get_variable(phylo,y)))!=2) stop("Just for 2-group")
  sample_data(phylo)$dependent=get_variable(phylo,y)
  a01=phylo %>>% subset_samples(!is.na(dependent)) %>>%
    (prune_taxa(taxa_sums(.)>0,.))
  a0=a01 %>%
    transform_sample_counts(function(x) x/sum(x)*10000) %>%
    psmelt() %>%
    dplyr::rename(COG=OTU)

  lm.summary=function(formula,data,Is_y_discrete){
    if(Is_y_discrete){
      b1=glm(formula,data=data, family = "binomial")
    } else {
      b1=lm(formula,data=data)
    }
    b2=summary(b1)
    b3=data.frame(b2$coefficients)
    return(data.frame(beta0=b3[1,1],p0=b3[1,4]
                      ,beta1=b3[2,1],p1=b3[2,4]
                      ,check.names = F)
    )
  }

  #  lm.summary(lm_formula,filter(a0,Genus=="Prevotella"))

  lm_formula=sprintf("%s~Abundance",y) %>% formula()
  anno=data.frame(tax_table(phylo)@.Data,stringsAsFactors = F) %>% rownames_to_column("COG")
  a1=split(a0,a0$COG) %>>%
    lapply(lm.summary,formula=lm_formula,Is_y_discrete=Is_y_discrete) %>>%
    bind_rows(.id="COG") %>%
    mutate(Significant=ifelse(p1<0.05,"Yes","No") %>% factor(levels=c("Yes","No"))
           ,sign=sign(beta1)) %>%
    left_join(anno,by="COG") %>%
    arrange(Significant,sign,p1)

  ##Volcano plot
  vol=ggplot(a1,aes(x=beta1,y=-log10(p1),colour=Significant))+geom_point(alpha=0.5)+xlab("beta")+ylab("-log10(p-value)")

  ##Scatter plot

  a3=group_by(a1,sign) %>>% top_n(10,-log10(p1)) %>>% arrange(sign,p1)
  a3$lab=sprintf("%s\nbeta=%.3e, p=%.3f",a3$COG,a3$beta1,a3$p1) %>% gsub("e\\+00","",.)

  p0=lapply(c(1,-1),function(x,Is_y_discrete){
    a4=filter(a3,sign==x) %>% ungroup() %>% select_("COG","lab")
    a5=inner_join(a0,a4,by="COG")
    a5$lab %<>% factor(levels=a4$lab)

    fig1=ggplot(data = a5, aes_string(x = "Abundance", y = y, label=y)) + ggtitle(ifelse(x==1,"Positive beta","Negative beta")) +
      geom_point() + facet_wrap(~lab,ncol=5,scales="free_x") +
      xlab("Normalized Relative Frequency") + ylab(y)
    if(Is_y_discrete){
      fig2=fig1+geom_smooth(method="glm", method.args = list(family = "binomial"),se=FALSE)
    } else {
      fig2=fig1+geom_smooth(method="lm",se=FALSE)
    }
    return(fig2)
  },Is_y_discrete=Is_y_discrete)

  sca=marrangeGrob(p0,nrow=2,ncol=1,top="")

  a0$COG %<>% factor(levels=a1$COG)
  a0$Sample %<>% factor(levels=sample_names(a01))
  a11=sprintf("%s~Sample","COG") %>% formula() %>>%
    (dcast(a0,.,value.var = "Abundance"))
  return(list(estimate=a1,scatter=a11,vol=vol,sca=sca))
}

function(){
#example
  phylo_file="E:/Biom/MS16049/greengene 20160714/MS16049_gg_phylo.rds"
  b1=read_rds(phylo_file)
  b0=read_rds("E:/Biom/MS16049/R code/20160718/fp_cog.rds")$fp
  sample_data(b0) <- sample_data(b1)

  phylo=b0 %>>% subset_samples(group%in%c("NP","ST")) %>>%
    (prune_taxa(taxa_sums(.)>0,.))

  a0=Functional_association_analysis(phylo,y="BMI",Is_y_discrete=F)
  a1=Functional_association_analysis(phylo,y="BMIgt25",Is_y_discrete=T)

  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunConVol.pdf",a0$vol)
  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunDisVol.pdf",a1$vol)
  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunConVol_xlim.pdf",a0$vol+xlim(-50,20))
  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunDisVol_xlim.pdf",a1$vol+xlim(-10,25))
  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunConSca.pdf",a0$sca,height=10,width=12)
  ggsave("E:/Biom/MS16036/R code/20160620 DEMO/AssAnaOutp/FunDisSca.pdf",a1$sca,height=10,width=12)

}
