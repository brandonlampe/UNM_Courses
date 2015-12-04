source("MatrixOp.R")
library(Matrix)

a <- .1
b <- .1
g <- 1

m <- g/(a*b)

phi <- .1

c <- cos(phi)
x <- sin(phi)

E.u <- 0.5 * matrix(c(2 *(a*c-1) + (a*c-1)^2 + (a*c-b*s)^2,
                a*s + a*c -b*s + a*s*(a*c - 1) + (a*c - b*s)*(a*s+b*c-1),
                0,
                a*s+a*c-b*s + a*s *(a*c -1 ) + (a*c - b*s)*(g*s + b*c -1),
                2*(b*s + b*c -1) +(a*s)^2 + (g*s + b*c -1)^2,
                0,
                0,0,0), ncol = 3)

E.F <- 0.5 * matrix(c(a^2 + b^2 -1, g*b, 0,
                g*b, b^2 -1, 0,
                0,0,0), ncol = 3)

E.F
E.u

F.eE <- matrix(c(a,0,0,g,b,0,0,0,1),ncol = 3)
F.eE

check <- t(F.eE) %*% e.F %*% F.eE
check

e.F <- 0.5 * matrix(c(1 - 1/a^2, g/(a^2 * b), 0,
                      g/(a^2*b), 1 - m^2 - 1/b^2,0,
                      0,0,0), ncol = 3)

e.F

e.u <- 0.5 * matrix(c(), ncol = 3)
e.u

#####################
# problem 5
#####################
par(mfrow=c(3,1))

lambda <- seq(0, 2, 0.025)
E.11 <- 0.5 * (lambda^2 - 1)
plot(lambda,E.11 , xlab = "Stretch_1", ylab = "E_11")
lines(lambda, E.11)

e.11 <- 0.5 * (1 - 1/lambda^2)
plot(lambda, e.11, xlab = "Stretch_1", ylab = "e_11")
lines(lambda, e.11)

U.ln <- log(lambda)
plot(lambda, U.ln, xlab = "Stretch_1", ylab = "ln(U)")
lines(lambda, U.ln)

mfpar()
