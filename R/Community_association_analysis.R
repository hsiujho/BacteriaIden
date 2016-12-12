#
#
#
#
#
#

Community_association_analysis=function(phylo,ranklv,y,Is_y_discrete=F){
  if(missing(phylo)) stop("phylo is missing")
  if(missing(ranklv)) stop("ranklv is missing")
  if(missing(y)) stop("y is missing")
  if(missing(Is_y_discrete)) stop("Is_y_discrete is missing")
  if(!is.logical(Is_y_discrete)) stop("Is_y_discrete should be logical")
  if(Is_y_discrete & length(table(get_variable(phylo,y)))!=2) stop("Just for 2-group")
  rankv=rank_names(phylo)[1:which(rank_names(phylo)==ranklv)]
  sample_data(phylo)$dependent=get_variable(phylo,y)
  a01=phylo %>>% subset_samples(!is.na(dependent)) %>>%
    (prune_taxa(taxa_sums(.)>0,.))
  a0=a01 %>%
    transform_sample_counts(function(x) x/sum(x)*100) %>>% {
      if(ranklv%in%rank_names(phylo)){
        tax_glom(.,ranklv,NArm=F)
      }
    } %>>%
    psmelt() %>>%
    mutate_each_(funs(as.character),rankv)

  for(i in 2:length(rankv)){
    j=rankv[i]
    j1=rankv[i-1]
    j2=sprintf(";%s_",tolower(substr(j,1,1)))
    a0[[j]]=ifelse(a0[[j]]=="",paste0(a0[[j1]],j2),a0[[j]])
  }

  a2=a0 %>>%
    group_by_(ranklv) %>>%
    summarise(N=n(),mean.RA=mean(Abundance),sd.RA=sd(Abundance))

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

  a1=split(a0,a0[[ranklv]]) %>>%
    lapply(lm.summary,formula=lm_formula,Is_y_discrete=Is_y_discrete) %>>%
    bind_rows(.id=ranklv) %>%
    mutate(Significant=ifelse(p1<0.05,"Yes","No") %>% factor(levels=c("Yes","No"))
           ,sign=sign(beta1)) %>%
    arrange(Significant,sign,p1)

  ##Volcano plot
  vol=ggplot(a1,aes(x=beta1,y=-log10(p1),colour=Significant))+geom_point(alpha=0.5)+xlab("beta")+ylab("-log10(p-value)")

  ##Scatter plot

  a3=group_by(a1,sign) %>>% top_n(10,-log10(p1)) %>>% arrange(sign,p1)
  a3$lab=sprintf("%s\nbeta=%.3e, p=%.3f",a3[[ranklv]],a3$beta1,a3$p1) %>% gsub("e\\+00","",.)

  p0=lapply(c(1,-1),function(x,Is_y_discrete){
    a4=filter(a3,sign==x) %>% ungroup() %>% select_(ranklv,"lab")
    a5=inner_join(a0,a4,by=ranklv)
    a5$lab %<>% factor(levels=a4$lab)

    fig1=ggplot(data = a5, aes_string(x = "Abundance", y = y, label=y)) + ggtitle(ifelse(x==1,"Positive beta","Negative beta")) +
      geom_point() + facet_wrap(~lab,ncol=5,scales="free_x") +
      xlab("Relative Abundance (%)") + ylab(y)
    if(Is_y_discrete){
      fig2=fig1+geom_smooth(method="glm", method.args = list(family = "binomial"),se=FALSE)
    } else {
      fig2=fig1+geom_smooth(method="lm",se=FALSE)
    }
    return(fig2)
  },Is_y_discrete=Is_y_discrete)

  sca=marrangeGrob(p0,nrow=2,ncol=1,top="")

  a0[[ranklv]] %<>% factor(levels=a1[[ranklv]])
  a0$Sample %<>% factor(levels=sample_names(a01))
  a11=sprintf("%s~Sample",ranklv) %>% formula() %>>%
    (dcast(a0,.,value.var = "Abundance"))
  return(list(estimate=a1,scatter=a11,vol=vol,sca=sca))
}

function(){
#example
  phylo_file="E:/Biom/MS16049/greengene 20160714/MS16049_gg_phylo.rds"
  b0=read_rds(phylo_file)

  phylo=b0 %>>% subset_samples(group%in%c("NP","ST")) %>>%
    (prune_taxa(taxa_sums(.)>0,.))

  a0=Community_association_analysis(phylo,ranklv="Genus",y="BMI",Is_y_discrete=F)
  a1=Community_association_analysis(phylo,ranklv="Genus",y="BMIgt25",Is_y_discrete=T)
}
