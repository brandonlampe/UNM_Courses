library(matrixcalc)
library(Matrix)

B <- 5.905e-6
A <- .45

m <- 2
n <- 10


num <- (m+1)*A*m
den <- (n+1)*B*n
R<- (num/den)^(1/(m-n))
R

F <- (m*A)/R^(m+1) -(n*B)/R^(n+1)
F

S <- R^4/(A*(n-1))
S
