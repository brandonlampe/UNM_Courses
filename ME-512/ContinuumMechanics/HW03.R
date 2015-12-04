source("CrossProd.R")
library(Matrix)

# problem 1

T <- matrix(c(-1, 2, 4, 2, -2, 3, 3, 2, 4), ncol = 3)
T

u <- matrix(c(1, -2, 2), ncol = 1)
v <- matrix(c(-2, 1, -3), ncol = 1)
w <- matrix(c(-1, 1, 2), ncol = 1)

A <- matrix(c(3,2,1,2,4,2,1,2,1), ncol = 3)
A

T1 <- sum((A %*% u) * CrossProduct3D(v,w))
T2 <- sum(u * CrossProduct3D(sum(A %*% v),w))
T3 <- sum(u * CrossProduct3D(v,sum(A %*% w)))

L <- T1 + T2 + T3
L

trA <- 3 + 4 + 1

T4 <- sum(u * CrossProduct3D(v,w))

R <- L * T4

# ---- a ----
sum(u*v)

# ---- b ----
t(v) %*% T %*% u

# ---- d ----
 t(u) %*% T %*% v

# ---- iv
Tsym <- .5*(T + t(T))
Tsym

Tskw <- .5*(T - t(T))
Tskw

Tsym %*% Tskw

# ---- 2
v <- matrix(c(-2, 1, -3), ncol = 1)
v

n <- v/(sum(v^2))^0.5
n

u <- sum(v*n) # dot product

sum(n*v) * n

sum(n*n) * v

sum(n*n)


n %*% u

CrossProduct3D(v,n)

# --- problem 5

#
Ee <- matrix(c(pi/2, pi/4, pi/4, pi/4, pi/3, pi*2/3, pi*3/4, pi/3, pi*2/3),
               ncol = 3)
Ee

ve <- matrix(c(1, -2, 3), ncol = 1)
ve
Tee <- matrix(c(2, 0, -3, 0, 6, 0, -3, 0, 4), ncol = 3)
Tee
# direction cosine matrix
AeE <- cos(eE)
AeE

AEe <- t(AeE)

det(AEe)

AEe %*% AEe # yield the identity matrix if ortho normal

E1 <- AEe[,1]
E2 <- AEe[,2]
E3 <- AEe[,3]

sum(E1 * E1)
sum(E2 * E2)
sum(E3 * E3)

CrossProduct3D(E1,E2)
CrossProduct3D(E2,E3)
CrossProduct3D(E3,E1)


e1 <- AEe[1,0]
e2 <- AEe[2,0]
e3 <- AEe[3,0]

# ---- c ----
vE <- AEe %*% ve
vE

ve <- AeE %*% vE
ve

# ---- d ----

TEE <- AEe %*% Tee %*% AeE
TEE

# ---- e ----
TEe <- AEe %*% Tee
TEe
TEE %*% AEe # check = TEe

TeE <- Tee %*% AeE
TeE
AeE %*% TEE # check = TeE

# ---- 5 ----

# -- a --
aEe <- matrix(c(cos(.25*pi), -sin(.25 * pi), 0,
                sin(.25 * pi), cos(.25*pi), 0,
                0, 0, 1), ncol = 3)
aEe

aEe %*% t(aEe)

# -- b --
aEg <- matrix(c(1,0,0,0,cos(.25*pi), sin(.25 * pi), 0,
                -sin(.25 * pi), cos(.25*pi)), ncol = 3)
aEg
aEg %*% t(aEg)

# -- c --

aeg <- t(aEe) %*% aEg
aeg

aeg %*% t(aeg)
