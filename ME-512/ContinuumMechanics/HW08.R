source("MatrixOp.R")
library(Matrix)

alpha <- 2
beta <- 3
gamma <- 1

phi <- pi/4

x <- c(alpha, beta, 1)
X <- 1/x

F.eE <- x * Diagonal(3)
F.Ee <- X * Diagonal(3)
F.eE %*% F.Ee

F.eE

c <- cos(phi)
s <- sin(phi)

# a.eE <- matrix(c(c, -s, 0, s, c, 0, 0, 0, 1), nrow = 3)
# a.Ee <- t(a.eE)
U.eE <- (t(F.eE) %*% F.eE)^(1/2)
U.eE

U.eE.inv <- solve(U.eE)
U.eE.inv %*% U.eE # should == I
U.eE.inv

R.eE <- F.eE %*% U.eE.inv
R.eE

R.eE %*% U.eE # should == F.eE

R.eE %*% t(R.eE) # should == I
################################
# problem 3

P.ee <- matrix(c(1, 3, 0, 0, 4, 0, 0, 0, 1), nrow = 3)
P.ee

P.ee.inv <- solve(P.ee)

test <-P.ee.inv %*% Lambda %*% P.ee

P <- eigen(test)$vectors
P.inv <- solve(P)

test.2 <- P %*% Lambda %*% P.inv

eigen(test.2)

alpha <- 2
beta <- 3^(1/2)
gamma <- 1

F.eE <- matrix(c(alpha, 0, 0, gamma, beta, 0, 0 ,0,1), nrow =3)
F.eE

U.eE.2 <- t(F.eE) %*% F.eE
U.eE.2

eigen(U.eE.2)
U.pp <- t((eigen(U.eE.2)$values)^(1/2)) %*% Diagonal(3)
U.pp <- matrix(c(U.pp[1], 0, 0, 0, U.pp[2], 0, 0, 0, U.pp[3]),nrow = 3)
U.pp
U.pp.vect <- eigen(U.eE.2)$vectors

U.eE <- t(U.pp.vect) %*% U.pp %*% U.pp.vect
U.eE

U.pp.inv <- solve(U.pp)
t(U.pp.vect) %*% U.pp.inv %*% U.pp.vect


U.eE.inv


R.eE <- F.eE %*% solve(U.eE)
R.eE
R.eE %*% t(R.eE)

R.eE %*% U.eE

V.eE.2 <- (F.eE) %*% t(F.eE)
V.eE.2

eigen(V.eE.2)
V.pp <- t((eigen(V.eE.2)$values)^(1/2)) %*% Diagonal(3)
V.pp
V.pp.vect <- eigen(V.eE.2)$vectors

V.pp.vect %*% U.pp.vect
