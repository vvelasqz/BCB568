---
title: "HW2_BCB568"
author: "Valeria Velasquez"
date: "January 29, 2018"
output:
  html_document: default
  pdf_document: default
---

##Part II: Combining multiple samples

In this part, we will extend the analysis to n individuals. Suppose the individuals have been randomly
sampled from a population where the reference nucleotide relative frequency is (unknown).

6. What is Hardy-Weinberg equilibrium (HWE), and what does it imply about the genotype Gi of the ith
individual?

In population genetics this principle states that without disturbances, the allele frequences will remain constant over the generations. This means that without any evolutionary force making an effect and with random mating, the genetic variation will remain constant in a population. We can express the frequencias as a multinomial distribution, and for a bialelic case it becomes into a binomial distribution. If we define m as the ploidy, G=g as the genotype, and $\psi$ as the allele frequency we have:

\[
 P(G=g|\psi) =   \binom{m}{g} \cdot \psi^g (1-\psi)^{m-g}
\]

This equation is telling us that that the probability of the genotype of the individual Gi=g can be modelled as a binomial density function of \psi as p, m and g. 

7. If we observe data $D_1=d_1,D_2=d_2,ldots,D_n=d_n$ and all the previous assumptions are true along with HWE, what is the likelihood function $L(\psi \mid d_1,d_2,ldots,d_n)$?

Using the likehood definition, with n as the total number of individuals and the law of total probability for G=g we can write:

\[
 L(\psi|di) = \prod_{i=1}^n fx(di|\psi)  =  \prod_{i=1}^n P(di|\psi) \cdot P(\psi)= \prod_{i=1}^n \sum_{g=0}^m P(di|G=g) \cdot P(G=g|\psi)
\]

Taking the log we end up wth:

\[
 log[L(\psi|di)] = l(\psi|di)= \sum_{i=1}^n log(\sum_{g=0}^m P(di|G=g) \cdot P(G=g|\psi))
\]

8. Consider the data, ngs_site_data.Rtxt, available on Canvas, and shown partially below. The columns
are individual i, read j, quality score q, and allele a, where allele 1 indicates the reference allele r and 0 the non-reference allele.

a. Write an R function, e.g.log.likelihood(psi, data), that computes the log likelihood.
Use it to plot the log likelihood as a function of $\psi$.

we start with importing the dataset into a data frame

```{r}
#setwd("~/iowa state/semester 4/BCB568/hw2")
data <- read.table('ngs_site_data.Rtxt', header = TRUE)
```


We create a fucntion log.likelihood in R that depends of x which represents $\psi$ and ngs_site data set

```{r}
log.likelihood<- function(x, ngs_site){ 
  read_index=0
  sum_products_ind=0
  for (ind in unique(ngs_site$i)){
    prod_0=1
    prod_1=1
    prod_2=1
    for (read in 1:length(ngs_site$i[ngs_site$i==ind]))
      {
      read_index= read_index +1 
      if (ngs_site$a[read_index]==1){
        prod_0 = prod_0 * (10**(-ngs_site$q[read_index]/10))/3*(1-x)**2
        prod_1 = prod_1 * (1/2*(1-10**(-ngs_site$q[read_index]/10))+ 1/2* ((10**(-ngs_site$q[read_index]/10))/3))* (2*x*(1-x))
        prod_2= prod_2* (1-10**(-ngs_site$q[read_index]/10))*(x)**2
      }
      else{
        prod_0 = prod_0 * (1-10**(-ngs_site$q[read_index]/10))*(x-1)**2
        prod_1= prod_1 * (1/2*(1-10**(-ngs_site$q[read_index]/10))+ 1/2* ((10**(-ngs_site$q[read_index]/10))/3))* (2*x*(1-x))
        prod_2= prod_2 * (10**(-ngs_site$q[read_index]/10))/3 *(x)**2
      }
      }
    sum_products_ind= sum_products_ind + log(prod_0+prod_1+prod_2)
  }
return(sum_products_ind)
}
```

After this we can do our graph by creating a vector for $\psi$, whcich we will call P between [0,1], and evaluate the log likelihood function on each of those values, to finally get a plot

```{r}
P <- seq(0, 1, length=101)
likelihood_vector = integer(101)
 for (i in 1:length(P)){
   likelihood_vector[i] = log.likelihood(P[i], data)
 }
plot(P, likelihood_vector, xlab="Psi", ylab="Log.Likelihood", main="Log likelhood function")

```

b. Use R optim() function to find the maximum likelihood estimate $\hat\psi$ for the population reference relative allele frequency.
You can use method="Brent" as suggested by (Li, 2011).

Now we optimize, using the optim function and Brent method

```{r}
psi_estimator=optim(0, lower = 0, upper = 1, log.likelihood, ngs_site=data, method = 'Brent', control=list(fnscale=-1))
psi_estimator
```

The results show that the maximum likehood of the function is given for a value of the parameter  $\hat\psi= 0.7595865$ with a value $l(\hat\psi)= -123.718$

c. Do you know of any way to estimate the variance in the estimate, $\text{Var}(\hat\psi)$?

We can calculate the variance of the MLE estimate $\text{Var}(\hat\psi)$. By using the information function I, the variance of $\hat\psi$ is asymptotically equivalent to:

\[
	{Var}(\hat\psi) = \frac{1}{nI(\hat\psi)}
\]

With $I(\psi)$ calculated as:

\[
	I(\psi) = -E(\frac{\partial^2 l(\psi|di)}{\partial \psi^2})
\]

##Part III: Handling uncertainty

The data from Part II are nested, reads (low) within individuals (high).
It is not obvious (to anyone) how to perform bootstrap resampling in this case (Ren, 2010) did some work to show that it is best to sample with replacement at the highest level (the usual bootstrap), but to sample without replacement at the lower levels. In other words, simply sample individuals with replacement, leaving their entire data intact.

a. Plan and implement this approach to obtain bootstrap estimates $\hat\psi^{(b)}$ from these data.

For doing this we can do a sample at the individual level and compute the log likelihood funtion from there

```{r}
bootstrap<- function(ngs_site){ 
  bootstrap_resample <- sample(unique(ngs_site$i), replace = TRUE)
  new_ngs_site <- data.frame(old_i=numeric(), j=numeric(), q=numeric(), a=numeric(),i=numeric(), stringsAsFactors=FALSE)

  for (n in 1:length(bootstrap_resample)){
    b_ngs_site <- ngs_site[ngs_site$i==bootstrap_resample[n], ]
    b_ngs_site <- cbind(b_ngs_site, i=n)
    new_ngs_site = rbind(new_ngs_site, b_ngs_site)}
  row.names(new_ngs_site) <- seq(1,length(new_ngs_site$i),1)
  colnames(new_ngs_site)<- c('old_i', 'j', 'q', 'a','i')
  return(new_ngs_site)
}


```

Now we use the log.likelihood function to the bootstrap resample and repeat for 1000 times

```{r}
likelihood_bootstrap <- data.frame(psi=numeric(), likelihood=numeric())
for (b in 1:1000){
  w= bootstrap(data)
  psi_bootstrap <- optim(0, lower = 0, upper = 1, log.likelihood, ngs_site=w, method = 'Brent', control=list(fnscale=-1))
  likelihood_bootstrap = rbind(likelihood_bootstrap, c(psi_bootstrap$par,psi_bootstrap$value))
}
colnames(likelihood_bootstrap)<- c('psi','likehood')
likelihood_bootstrap<- likelihood_bootstrap[order(likelihood_bootstrap$psi),]
```

b. Make a histograph of the bootstrap estimates. Obtain a bootstrap confidence interval for $\psi$.
Any concerns with the analysis?

After ordering the bootstrap calculations we can plot the histogram 
```{r}
hist(likelihood_bootstrap$psi, breaks=20)
```

For a 1-$\alpha$ confidence interval we can take the bootstrap-t method that states that if we have a plug-in estimator T(x) of a parameter t(F) and an accompanying variance estimator V(x), we can calculate the confidence interval of t(F) as:

\[
	P\left[t(F)\geq T(x) -\sqrt{V(x)} \xi^\star_{1-\alpha}\right]\geq 1-\alpha
\]

Then we can calculate $T(x) -\sqrt{V(x)} \xi^\star_{1-\alpha}$ as an approximate level $1-\alpha$ lower confidence bound for t(F)

Now we have as T(x) our maximum likehood estimator $\hat\psi$, as the variance the variance of the bootstrap estimates and $\xi^\star_{1-\alpha}$ as the quantiles that correspond to that confidence value. For example if we take $\alpha=0.05$ we will take the ordered vector of bootstrap estimates indexes 25 and 975

```{r}
Variance_psi= var(likelihood_bootstrap$psi)
b_t= numeric(1000)
for (j in 1:1000) {
b_t[j] = (psi_estimator$par - likelihood_bootstrap$psi[j])/sqrt(Variance_psi)
}
b_quantiles<-quantile(b_t, probs = c(0.025, 0.975))
lower_limit= psi_estimator$par - sqrt(Variance_psi)*b_quantiles[2]
upper_limit=psi_estimator$par - sqrt(Variance_psi)*b_quantiles[1]
cat("bootstrap-t lower_limit:", lower_limit, "\n")
cat("bootstrap-t upper_limit:", upper_limit, "\n")
```

This gives us the confidence interval for $\psi=[0.5668, 0.9091]$

c. An alternative to the above procedure for obtaining a confidence interval is to utilize the likelihood ratio test of null hypothesis
\[
	H_0: \psi = \psi_0
\]
The likelihood ratio, with notation for this model and these data, is
\[
	\Lambda(d_1, d_2, \ldots, d_n) = \frac{L(\psi_0\mid d_1, d_2, \ldots, d_n)}{\text{sup}\{L(\psi \mid d_1, d_2, \ldots, d_n): \psi\in(0,1)\}},
\]

The likelihood ratio test statistic $\lambda := -2\ln\Lambda \stackrel{\cdot}{\sim} \chi^2_1$ has an asymptotic chi-squared distribution with one degree of freedom.
Large values of $\lambda$ indicate against $H_0$, and the approximation by the $\chi_1^2$ distribution becomes better for large $n$.
One can construct a confidence interval from those $\psi_0$ values for which $H_0$ is not rejected, which implies the confidence interval contains all $\psi_0$ where
\[
	-2\left[l(\psi_0\mid d_1, d_2, \ldots, d_n) - l(\hat\psi\mid d_1, d_2, \ldots, d_n)\right] \leq \xi_{1-\alpha},
\]
$\xi_p$ is the $p$th quantile of $\chi_1^2$ and $1-\alpha$ is the confidence level of the interval.
The boundaries of this interval can be found by solving for the two $\psi_0$ for which the above inequality becomes an equality.
This numerical problem is easily handled by the function uniroot() in R.

Find this likelihood ratio test-based confidence interval.Any concerns with the analysis?

We can define a function findroot, to find the roots as:
\[
	f(\psi_0)=-2\left[l(\psi_0\mid d_1, d_2, \ldots, d_n) - l(\hat\psi\mid d_1, d_2, \ldots, d_n)\right] - \xi_{1-\alpha}
\]

This function has teo parameters, the maximum likelihood estimation for the data at $\hat\psi$ which we calculated as $l(\hat\psi)= -123.718$, and Po that is the corresponding $\psi_0$ over which we will calculate a log.likehood function and then, we can use it to calculate the limits of the confidence interval
```{r}
find.root<-function(max_lik, Po) {
  return(-2*(log.likelihood(Po, data) - max_lik) - qchisq(p=0.95, df = 1))
    }
```

Then we use uniroot to calculate the limits $\psi_0=Po$

```{r}
save.lower<-uniroot(f = find.root, interval = c(0, 0.75), max_lik = psi_estimator$value)
save.upper<-uniroot(f = find.root, interval = c(0.75, 1), max_lik = psi_estimator$value)
```

This procedure gives us a confidence interval for $\psi=[0.6984, 0.8143]$




