This R package, diffShrinkHDR, is for the following published paper:

Drikvandi R. (2025). "High dimensional regression with many nuisance parameters: both cases of specified and unspecified parameters of interest". Electronic Journal of Statistics 19 (1), 2923–2957. https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-19/issue-1/High-dimensional-regression-with-many-nuisance-parameters--Both-cases/10.1214/25-EJS2401.full

You can install this R package using the following two commands in R:

library(devtools)

install_github("rezadrikvandi/diffShrinkHDR")

and then load it using library(diffShrinkHDR)

Details about this package can be found in the DESCRIPTION file. The main R function of the package is called "diffShrinkHDR" which can be applied for estimation and inference with the proposed method for high dimensional data. There is a test data to try, see below for details.

Example:

Test data: the test data are simulated with n=100 and p=200, where the first column is the response Y and the others are covariates X. Below is the command for applying the R package to this test data along with the output:

diffShrinkHDR(X=as.matrix(testdata[,-1]),Y=c(testdata[,1]),Index_XI=NULL,alpha=0.05)

Here is the ourput:


