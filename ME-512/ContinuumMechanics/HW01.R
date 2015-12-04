# CONTINUUM MECHANICS
# ASSIGNMENT 1
# DUE AUGUST 28, 2014 (THURSDAY)

# ---- PROBLEM 1 ----
T <- matrix(c(-4, 2, 3, 2, 3, 1, 4, 3, 1), ncol = 3, byrow = TRUE)

u <- matrix(c(0, -1, 2))

v <- matrix(c(1, -2, -3))

# ---- 1.b.i ----
w <- T %*% v # inner product
print(w)

# ---- 1.b.ii ----
u %*% t(v)

# ---- 1.b.iii ----
t(v) %*% T %*% u

# ---- 1.b.iv ----
t(u) %*% T

# ---- 1.b.v ----
t(u) %*% T %*% v

# ---- 1.b.vi ----
library(matrixcalc)
matrix.trace(T)

# ---- 1.b.vii ----
T %*% T

# ---- 1.b.viii ----
t(T) %*% T

# ---- 1.b.ix ----
matrix.trace(t(T) %*% T)
