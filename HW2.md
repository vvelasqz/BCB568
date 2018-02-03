---
title: "HW2_BCB568"
author: "Valeria Velasquez"
date: "January 29, 2018"
output: html_document
---


```{r}
setwd("~/iowa state/semester 4/BCB568/hw2")
ngs_site <- read.table('~/iowa state/semester 4/BCB568/hw2/ngs_site_data.Rtxt', header = TRUE)
```

We create a fucntion in R that depenst of P and ngs_site data 

```{r}
log.likelihood<- function(x, ngs_site){ 
  sum_total=0
  for (i in 1:length(ngs_site$i)){
    sum_site=0
    if (ngs_site$a[i]==1)
   {
      sum_site = (10**(-ngs_site$q[i]/10))/3 *(x**2) + (1/2*(1-10**(-ngs_site$q[i]/10))+ 1/2* ((10**(-ngs_site$q[i]/10))/3))* (2*x*(1-x))+ (1-10**(-ngs_site$q[i]/10))*(1-x)**2
   
   }
    else
   {
     sum_site =  (1-10**(-ngs_site$q[i]/10))*(x**2) + (1/2*(1-10**(-ngs_site$q[i]/10))+ 1/2* ((10**(-ngs_site$q[i]/10))/3))* (2*x*(1-x))+ (10**(-ngs_site$q[i]/10))/3 *(1-x)**2
   
   }
    sum_total= sum_total + log(sum_site)
  }
return(sum_total)
}
```

# http://www.montana.edu/rotella/documents/502/MarkdownEqnExamples.Rmd


```{r}
P=0
a <- log.likelihood(P, ngs_site)
```


now a loop with P[0,1]


```{r}
P <- seq(0, 1, length=101)
likelihood_vector = integer(101)
 for (i in 1:length(P)){
   likelihood_vector[i] = log.likelihood(P[i], ngs_site)
 }
plot(P, likelihood_vector, xlab="Psi", ylab="Log.Likelihood", main="Log likelhood function")

```



