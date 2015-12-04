#### HW 02
library(ggplot2)
############# PROBLEM 3.1
epsilon <- matrix(c(0.02, 0.1))
epsilon[2,1]

sigma <- matrix(c(175, 185))
sigma

n <- (log(sigma[2,1]) - log(sigma[1,1])) /
    ((log(epsilon[2,1]) - log(epsilon[1,1])))

k <- sigma[1,1]/ epsilon[1,1]^n

hollomon <- function(e, n, k){
          sigma <- k * e ^ n
          }

check <- hollomon(epsilon, n, k)

e.test <- seq(0.0,0.5, 0.005)
s.check <- hollomon(e.test, n, k)
plot(epsilon, sigma)
lines(e.test, s.check)

######## PROBLEM 3.16

# compressive stress / tensile stress [MPa]
sc.st <- matrix(c(0.9, .8, .7, .6, .5, .4, .35, .325, .32))

# x 10^(-2)
ep <- matrix(c(0.1, .25, .4, .65, .95, 1.5, 2, 2.5, 3))

e <- ep * 10 ^ (-2)

plot(ep, sc.st)
lines(ep, sc.st)

fsigma.T <- function (e){
        300 + 450 * e ^ 0.5 }

sigma.C <- sc.st * fsigma.T(e)
sigma.C

sigma.T <- fsigma.T(e)

## final Plot
par(mfrow = c(2,1))
plot(e, sigma.C, xlab = "Strain",
     ylab = "Compressive Stress [MPa]")
lines(e, sigma.C)

plot(e, sigma.T, xlab = "Strain",
     ylab = "Tensile Stress [MPa]")
lines(e,sigma.T)
