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
for (ind in 1:max(ngs_site$i)){
    sum_products_ind=0
    prod_0=1
    prod_1=1
    prod_2=1
    for (read in 1:length(ngs_site$i[ngs_site$i==ind]))
      {
      if (ngs_site$a[(ind-1)*length(ngs_site$i[ngs_site$i==ind-1])+read]==1)
        {
        prod_0 = prod_0 * (10**(-ngs_site$q[i]/10))/3*(1-x)**2
        prod_1 = prod_1 * (1/2*(1-10**(-ngs_site$q[i]/10))+ 1/2* ((10**(-ngs_site$q[i]/10))/3))* (2*x*(1-x))
        prod_2= prod_2* (1-10**(-ngs_site$q[i]/10))*(x)**2
        }
      else
        {
          prod_0 = prod_0 * (1-10**(-ngs_site$q[i]/10))*(x-1)**2
          prod_1= prod_1 * (1/2*(1-10**(-ngs_site$q[i]/10))+ 1/2* ((10**(-ngs_site$q[i]/10))/3))* (2*x*(1-x))
          prod_2= prod_2 * (10**(-ngs_site$q[i]/10))/3 *(x)**2
          }
          
    
  }

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



