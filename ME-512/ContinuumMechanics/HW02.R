source(CrossProd)
library(Matrix)

# ---- 1.i
v <- c(3,2,1)

# ---- 1.ii
# u <- c(-1, -2, 7)
u <- c(2, 2, -10)
# ---- check that u is orthogonal to v
sum(u*v) # must = 0

# ---- 1.iii
Mv <- sum(v^2)^(1/2)
v/Mv

E1 <- v / sum(v^2)^(1/2)
E1

E2 <-  u / sum(u^2)^(1/2)
E2

E3 <- CrossProduct3D(E1, E2)
E3
# ---- 1.iv
e <- Diagonal(3)
E <- rbind(E1, E2, E3)
E

sum(E1 * e[1,])
sum(E1 * e[2,])
sum(E1 * e[3,])

sum(E2 * e[1,])
sum(E2 * e[2,])
sum(E2 * e[3,])

sum(E3 * e[1,])
sum(E3 * e[2,])
sum(E3 * e[3,])

AEe <- E

# ===============
# problem 2

det(AEe)

AEe %*% t(AEe)

sum(E1*E2)
sum(E1*E3)
sum(E3*E2)

# ================
# Prob 3
wE <- matrix(c(3, 4, 5), nrow = 3)
wE.mag <- (sum(wE^2))^(0.5)
wE.mag

wE
AeE <- matrix(t(AEe), ncol = 3)
AeE

we <- AeE %*% wE

# ===============
# Prob 4
we
we.mag <- (sum(we^2))
we.mag

# ================
# Prob 5
Ze <- matrix(c(7,8,9),nrow = 3)
Ze

sum(we * Ze)


