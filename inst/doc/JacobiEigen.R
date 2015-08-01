## ---- comment = ""-------------------------------------------------------
imod <- aov(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, iris)
(R <- cor(resid(imod)))

library(JacobiEigen)
suppressPackageStartupMessages(library(dplyr))
rEig <- JacobiR(R) 
cEig <- JacobiCpp(R)
identical(rEig, cEig)  ## the R and Rcpp implementations are identical
cEig
(eEig <- eigen(R))
all.equal(eEig$values, cEig$values)  ## eigenvalues are (practically) identical
crossprod(eEig$vectors, cEig$vectors) %>% ## eigenvectors differ in signs
  round(10) 

## ---- comment = ""-------------------------------------------------------
library(microbenchmark)
microbenchmark(JacobiR(R), JacobiCpp(R), eigen(R))

## ---- comment = ""-------------------------------------------------------
suppressPackageStartupMessages(library(tidyr))
set.seed(1234)
N <- 100
iseq <- seq(5, 50, by = 5)
res <- lapply(iseq,  function(n) {
  S <- crossprod(matrix(rnorm(N*n), N, n))/N
  runTime <- microbenchmark(JacobiCpp(S), eigen(S), times = 20)
  c(n = n, with(runTime, tapply(time, expr, median))/1000)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>% 
  gather(key = expr, value = time, `JacobiCpp(S)`, `eigen(S)`)

suppressPackageStartupMessages(library(ggplot2))
ggplot(res) + aes(x = n, y = log10(time), colour = expr) + geom_line() + geom_point() +
  theme(legend.position = "top") + xlab("matrix size") +
  ylab("log10(median run time in milliseconds)")

