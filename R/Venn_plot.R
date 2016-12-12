#
#
#
#
#
#

Venn_plot=function(phylo,rankname,groupvar,...){
  a2= phylo %>>% tax_glom(rankname) %>>% psmelt() %>>% filter(Abundance!=0) %>>%
    select_(rankname,groupvar) %>>% mutate(Y=1) %>>% unique() %>>%
    dcast(formula(sprintf("%s~%s",rankname,groupvar)),value.var="Y",fill=0)

  likes <- function(group) {
    ppl <- a2
    for (i in 1:length(group)) {
      ppl <- subset(ppl, ppl[group[i]] == T)
    }
    nrow(ppl)
  }


  plotVenn <- function(a, ...) {
    grid.newpage()
    if (length(a) == 1) {
      out <- draw.single.venn(likes(a), ...)
    }
    if (length(a) == 2) {
      out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
    }
    if (length(a) == 3) {
      out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]),
                              likes(a[2:3]), likes(a[c(1, 3)]), likes(a), ...)
    }
    if (length(a) == 4) {
      out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]),
                            likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]),
                            likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2,
                                                                                       4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), likes(a), ...)
    }
    if (length(a) == 5) {
      out <- draw.quintuple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), likes(a[5]), likes(a[1:2]), likes(a[c(1,3)]), likes(a[c(1,4)]), likes(a[c(1,5)]), likes(a[2:3]), likes(a[c(2,4)]), likes(a[c(2,5)]), likes(a[3:4]), likes(a[c(3,5)]), likes(a[c(4,5)]), likes(a[1:3]), likes(a[c(1,2,4)]), likes(a[c(1,2,5)]), likes(a[c(1,3,4)]), likes(a[c(1,3,5)]), likes(a[c(1,4,5)]), likes(a[2:4]), likes(a[c(2,3,5)]), likes(a[c(2,4,5)]), likes(a[c(3,4,5)]), likes(a[c(1,2,3,4)]), likes(a[c(1,2,3,5)]), likes(a[c(1,2,4,5)]), likes(a[c(1,3,4,5)]), likes(a[c(2,3,4,5)]), likes(a[1:5]), ...)
    }
    if (!exists("out"))
      out <- "Oops"
    return(out)
  }
  get_variable(phylo,groupvar) %>>% unique() %>>% {
    if(is.factor(.)){
      levels(.)
    } else {
      unique(.)
    }
  } %>>% plotVenn(...)
  return(a2)
}
