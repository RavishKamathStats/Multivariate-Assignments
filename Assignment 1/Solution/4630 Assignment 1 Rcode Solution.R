##############
# Question 1 #
##############

#input data into A
A = matrix(c(2, 1 , 0, 1, 4, 1, 0, 1, 4), ncol=3, byrow=T)

#(a1) calculate eigenvalues and eigenvectors
result = eigen(A)
result

#(a2) calculate Cholesky decomposition
L = t(chol(A))
L

#(a3) calculate spectral decomposition

#D is a diagonal matrix with diagonal being the eigenvalues
D = diag(result$values,  ncol=3)

#P is the matrix with the i-th column being the i-th eigenvalues
P = result$vectors

#spectral decomposition
sqrtA = P %*% sqrt(D) %*% t(P)
sqrtA



##############
# Question 4 #
##############

#set up the necessary packages
library(xlsx2dfs)
library(DescTools) 

#set up the data matrix X
x1 = c(3, 3, 4, 5, 6, 8)
x2 = c(17.95, 15.54, 14.00, 12.95, 8.94, 7.49)
X = matrix(c(x1, x2), ncol=2, byrow=F)
n = nrow(X)
p = ncol(X)

#calculate the mean and variance of X
xbar = apply(X, 2, mean)
S = var(X)

#calculate the distance vector
centerx = sweep(X, 2, xbar, "-")
d = c(diag(centerx %*% solve(S) %*% t(centerx)))
d

#calculate the table value for a 50% confidence contour of the mean
qchisq(0.5, p)


#set up the hypothesised value mu_0 and perform Hotelling T^2 test using R
mu0 = c(3, 10)
HotellingsT2Test(X, mu=mu0)

#performing the calculations by hand
barx = matrix(xbar, ncol=1)
testmu = matrix(mu0, ncol=1)

Tsq = t(barx-testmu) %*% (n*solve(S)) %*% (barx - testmu)
Tsq

obsteststat = (n-p)/((n-1)*p)*Tsq
obsteststat

pvalue = 1 - pf(obsteststat,p, n-p)
pvalue

##############
# Question 5 #
##############

#start the library such that you can read in the data
library(xlsx2dfs)
library(DescTools) 
data <- read.xlsx(file.choose(), 1) 

#I prefer to rename the data
North = data$North
South = data$South
East = data$East
West = data$West

#put 4 normal Q-Q plots on same page 
par(mfrow=c(2, 2))
qqnorm(North, main="Q-Q plot of North")
qqnorm(South, main="Q-Q plot of South")
qqnorm(East, main="Q-Q plot of East")
qqnorm(West, main="Q-Q plot of West")

#use Shapiro-Wilks test to test if the normality assumption is violated
shapiro.test(North)
shapiro.test(South)
shapiro.test(East)
shapiro.test(West)

#obtain the t-test for each variable separately
t.test(North)
t.test(South)
t.test(East)
t.test(West)

#form the multivariate data matrix X
X = cbind(North, South, East, West)
n = nrow(X)
p = ncol(X)

#calculate the mean vector, variance matrix and correlation matrix
databar = apply(data, 2, mean)
databar
S = var(X)
S
cor(X)

#d_i = (x_i - xbar)' S^{-1} (x_i - xbar)
centerdata = sweep(X, 2, databar, "-")
d = c(diag(centerdata %*% solve(S) %*% t(centerdata)))

#arrange d_i in ascending order and obtain the theoretical quantiles
orderd = sort(d)
quantile = qchisq((seq(1, n, 1)-0.5)/n, p)

#obtain the multivariate normal Q-Q-plot
par(mfrow=c(1,1))
plot(quantile, orderd, main="Multivariate normal Q-Q plot", 
     xlab="Theoretial quantile", ylab="Sample quantile")

#set up the hypothesised value mu_0 and perform Hotelling T^2 test 
mu0 <- c(1450, 1900, 1700, 1700)
HotellingsT2Test(X, mu=mu0)

#check if mu_0 falls within the 95% confidence region of mu
critF = qf(0.95, p, n-p)
critF
comparevalue = (n-p)/((n-1)*p)*n*t(databar-mu0)%*%solve(S)%*%(databar-mu0)
comparevalue


