# ==== ADA1 - FA2014 - BRANDON LAMPE
# ====Start HW 01 ======
# load ggplot
library(ggplot2)

#==== Problem 1 ====
#---- (a) ----
# download monthly precipitation and sulphur concentration
# data from Univ. of Stockholm
d1 <- read.csv("http://statacumen.com/teach/ADA1/ADA1_HW_01_F14-1.csv")

# assign columns in data frame to individual variables
Month <- d1$Month
Precip <- as.numeric(d1$Precip)   # change type to numeric
Sulphur <- as.numeric(d1$Sulphur) # change type to numeric

# stem and leaf plot of Precip
stem(Precip, scale = 2)

# histogram of Precip
Precip.hist <- ggplot(d1, aes(x = Precip))
Precip.hist <- Precip.hist + geom_histogram(binwidth = 5)
Precip.hist <- Precip.hist + labs(title = "Monthly Precipitation")
print(Precip.hist)

# boxplot of Precip
Precip.box <- ggplot(d1, aes(x = "in", y = Precip)) # boxplot of Precip
Precip.box <- Precip.box + geom_boxplot()
Precip.box <- Precip.box + coord_flip()
Precip.box <- Precip.box + labs(title = "Monthly Precipitation")
print(Precip.box)

# ---- (b) ----
mean(Precip)    # mean of precip
median(Precip)  # median of precip
sd(Precip)      # standard deviation of precip
IQR(Precip)     # inter quartile range (range for middle half of data)

Precip.den <- ggplot(d1, aes(x = Precip))
Precip.den <- Precip.den + geom_histogram(aes(y = ..density..),
                                          binwidth = 5,
                                          color = "black",
                                          fill = "white")
Precip.den <- Precip.den + labs(title = "Monthly Precipitation")
Precip.den <- Precip.den + geom_density(alpha = 0.1)
print(Precip.den)

library(ggplot2)
mass.hist <- ggplot(mass.df, aes(x = mass))
mass.hist <- mass.hist + geom_histogram(binwidth = 10000)
mass.hist <- mass.hist + labs(title = "mass [gram]")
mass.hist
