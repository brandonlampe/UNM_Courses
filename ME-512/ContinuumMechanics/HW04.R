source("MatrixOp.R")
library(Matrix)

Tee <- matrix(c(3, 0, -1, 0, 1, 0, -1, 0, 3),ncol = 3)

Tee %*% t(Tee)

# ==== 1
Tee2 <- Tee %*% Tee
Tee3 <- Tee2 %*% Tee

# === 2
# the trace
IT <- sum(diag(Tee))
IIT <- sum(diag(Tee2))
IIIT <- sum(diag(Tee3))

# ====

ITS <- IT
IITS <-0.5*(IIT - IT^2)
IIITS <- (1/6)*(IT^3 - 3*IT * IIT + 2* IIIT)

#---- check, should equal IIITS ----
det(Tee)
p <- eigen(Tee)$vectors
lambda <- eigen(Tee)$values

zero <- matrix(c(0,0,0), ncol = 1)
lambda1 <- diag(4,nrow = 3, ncol = 3)
lambda1 <- diag(2,nrow = 3, ncol = 3)
lambda1 <- diag(1,nrow = 3, ncol = 3)

p1 <- matrix(c(sqrt(2)/2,0,-sqrt(2)/2), ncol = 1)
p2 <- matrix(c(sqrt(2)/2,0,sqrt(2)/2), ncol = 1)
CrossProduct3D(p1, p2)

(Tee - lambda1) %*% p1
lambda * p[,1]
Tlambda <- diag(lambda, 3,3)

p <- eigen(Tee)$vectors
sum(p[1,] * p[1,])
eigen(Tee)

# ---- problem 4
# ---- transformation matrix from e-e to p-p
Ape <- t(p %*% matrix(c(-1,0,0,0,1,0,0,0,-1),ncol = 3))
det(Ape) # should == 1

Tpp <- Ape %*% Tee %*% t(Ape)
Tpp
Tpp2 <- Tpp %*% Tpp
Tpp3 <- Tpp %*% Tpp %*% Tpp

IITp <- sum(diag(Tpp))

# ---- problem 5 ----
Tee <- t(Ape) %*% Tpp %*% Ape
Tee

# check cayley-homilton theorem
Tee3 <- Tee %*% Tee %*% Tee
Tee2 <- Tee %*% Tee

# ---- should == zeros
Tee3 - ITS * Tee2 - IITS * Tee - IIITS*diag(1,3,3)

# ---- problem 7
Troot <- t(Ape) %*% Tlambda %*% Ape

Tinvee <- solve(Tee)

Tinvpp <- Ape %*% Tinvee %*% t(Ape)

sum(diag(Tinvpp))
sum(diag(Tinvee))
Tinvpp
