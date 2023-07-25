# Median bias reduction in cumulative link models 
This package offers methods for fitting cumulative link models using adjusted score function for median bias reduction (Kenne Pagui et al., 2017, 2020). The proposal is illustrated in Gioia et al. (2023). 

To install the package use 
``` r
devtools::install_github("VinGioia90/mbrclm@master")
```

Here, we illustrate the proposal by using a modified version of the `wine` dataset available in the `ordinal` package and related to Randall (1989). In this respect, the central categories are merged into a single one. 

``` r
library(mbrclm)
data(wine)

n <- nrow(wine)
y <- matrix(NA, nrow = n, ncol = length(levels(wine$rating_comb3)))

for(i in 1 : n){
  if(wine$rating_comb3[i] == "1") y[i,]<-c(1, 0, 0)
  if(wine$rating_comb3[i] == "2 - 4") y[i,] <- c(0, 1, 0)
  if(wine$rating_comb3[i] == "5") y[i,] <- c(0, 0, 1)
}


x <- wine[, c("temp", "contact")]
X <- cbind(as.numeric(x[, 1]) - 1, as.numeric(x[,2])-1)
``` 

``` r
## maximum likelihood estimates ##
fit_ML <- mbrclm(x = X, y = y, type = "AS_ml", link = "logit")
summary(fit_ML)
#Coefficients (with logit link):
# Estimate Std. Error z value Pr(>|z|)
#[1,] -1.322e+00  5.304e-01  -2.492   0.0127 *
#[2,]  2.801e+01  1.023e+05   0.000   0.9998
#[3,] -2.581e+01  1.023e+05   0.000   0.9998
#[4,] -1.307e+00  7.175e-01  -1.821   0.0686 .
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Type of estimator: AS_ml (maximum likelihood)
#Number of iterations in the quasi-Fisher scoring: 27
#
#Convergence status: TRUE
```


The median bias reduced estimates can be obtained via
``` r
## median bias-reduced estimates ##
fit_medianBR <- mbrclm(x = X, y = y, type = "AS_median", link = "logit")
summary(fit_medianBR)
#Coefficients (with logit link):
#Estimate Std. Error z value Pr(>|z|)
#[1,]  -1.2893     0.5194  -2.482  0.01305 *
#[2,]   6.4622     2.3239   2.781  0.00542 **
#[3,]  -4.4784     2.2862  -1.959  0.05012 .
#[4,]  -1.2369     0.6757  -1.831  0.06716 .
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Type of estimator: AS_median (median bias-reduced)
#Number of iterations in the quasi-Fisher scoring: 18
#
#Convergence status: TRUE
```


# References

- Christensen, R. H. B. 2022. Ordinal - Regression Models for Ordinal Data. R package version 2022.11-16. http://CRAN.R-project.org/package=ordinal

- Kenne Pagui, E. C., A. Salvan, and N. Sartori. 2017. Median bias reduction of maximum likelihood estimates. Biometrika 104 (4):923–38.

- Kenne Pagui, E. C., A. Salvan, and N. Sartori. 2020. Efficient implementation of median bias reduction with applications to general regression models. arXiv: 2004.08630, https://arxiv.org/abs/2004.08630

- Randall, J. H. 1989. The analysis of sensory data by generalized linear model. Biometrical Journal 31 (7):781–93. 

