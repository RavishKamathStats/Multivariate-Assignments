#### Question 1 ####
A = matrix(c(2,1,0,1,4,1,0,1,2), nrow = 3, 
           ncol = 3, byrow = TRUE)
#EIGEN#
eigen(A)

#CHOLESKY#
t(chol(A))

#SPECTRAL#
ev = eigen(A)
L = ev$values
L
V = ev$vectors
V
D = diag(L)
D
sqrtD = sqrt(D)
sqrtD
sqrtA = V%*%sqrtD%*%t(V)
sqrtA
all.equal(A, zapsmall(sqrtA%*%t(sqrtA)) )



##### Question 4 ####
X = matrix(c(3, 17.95, 3, 15.54, 4, 14, 5, 12.95, 6, 8.94,
             8, 7.49), nrow = 6, ncol = 2, byrow = TRUE)
n = 6
p = 2

#Part A
vec1 = matrix(1, 6, 1)
xbar = 1/6*t(X)%*%vec1
xbar

#Part B
M = t(X)%*%X
L = xbar%*%t(xbar)
N = 6*L
S = 1/5*(M-N)
S

#Part C
S_inv = solve(S)
vec1 = matrix(1, 6, 1)
r = X[,1] - xbar[1,]
t = X[,2] - xbar[2,]
centered_mat = cbind(r,t)
distance = centered_mat%*%S_inv%*%t(centered_mat)
diag(distance)

#Part D.1
plot(xbar[1], xbar[2], type="p", xlim=c(0, 10), 
     ylim=c(5, 20), xlab="mu1", ylab="mu2")

mu1 = matrix(seq(-5, 20, 0.05), ncol=1, byrow=T)
nmu1 = nrow(mu1)
mu2 = matrix(seq(5, 20, 0.05), ncol=1, byrow=T)
nmu2 = nrow(mu2)


for (i in 1:nmu1) {
  for (j in 1:nmu2) {
    mu = matrix(c(mu1[i, 1], mu2[j, 1]), ncol=1, byrow=T)
    Fcomp = c((n-p)/((n-1)*2)*(n*t(xbar-mu)%*%solve(S)%*%(xbar-mu)))
    Fcrit = qf(0.50, p, n-p)
    if (Fcomp < Fcrit) points(mu1[i, 1], mu2[j, 1], pch="*")
  }
}
points(xbar[1], xbar[2], col='red')





#Part D.2
mu0 = matrix(c(3, 10), ncol=1, byrow=T)
Tobs = n*t(xbar-mu0)%*%S_inv%*%(xbar-mu0)
Tobs

Fcriticalvalue = (n-1)*p/(n-p)*qf(p = 0.05, df1 = p, df2 = n-p, lower.tail = FALSE)
Fcriticalvalue

pvalue = 1-pf((n-p)/((n-1)*2)*Tobs, p, n-p)
pvalue

#As we can see that since our observed Hotelling squared statistic is larger than the critical value,
#we can say that we will reject H_0 and say that there is evidence that the population mean is different 
#from mu_0




#### Question 5 ####
df = read.csv('/Users/ravishkamath/Desktop/University/2. York Math/1 MATH/1. Statistics /MATH 4630/3. Assessments/Assignments/Assignment 1/Dataset/MATH4630_a1data.csv')
View(df)
df = data.frame(df)
X = data.matrix(df)
n = dim(df)[1]
p = dim(df)[2]


#Part A
par(mfrow = c(2,2))

qqnorm(df$East, main = 'EAST')
qqline(df$East)

qqnorm(df$South, main = 'South')
qqline(df$South)

qqnorm(df$West, main = 'West')
qqline(df$West)

qqnorm(df$North, main = 'North')
qqline(df$North)


# Part B
xbar =1/n*t(X)%*%onemat
xbar
alpha = 0.05
degrees.freedom = n - 1
t.score= qt(p = alpha/2, df = degrees.freedom, lower.tail = F)
t.score

#North C.I.
sample.sd = sd(df$North)
sample.se = sample.sd/sqrt(n)
sample.se

lower_bound = xbar[1] - t.score*sample.se
upper_bound = xbar[1] + t.score*sample.se
c(lower_bound, upper_bound)
#Therefore the C.I. for North would be (1321.866, 1606.034)

#South C.I.
sample.sd = sd(df$South)
sample.se = sample.sd/sqrt(n)
sample.se

lower_bound = xbar[2] - t.score*sample.se
upper_bound = xbar[2] + t.score*sample.se
c(lower_bound, upper_bound)
#Therefore the C.I. for South would be (1726.799, 2050.401)

#South C.I.
sample.sd = sd(df$East)
sample.se = sample.sd/sqrt(n)
sample.se

lower_bound = xbar[3] - t.score*sample.se
upper_bound = xbar[3] + t.score*sample.se
c(lower_bound, upper_bound)
#Therefore the C.I. for East would be (1574.315, 1894.484)

#South C.I.
sample.sd = sd(df$West)
sample.se = sample.sd/sqrt(n)
sample.se

lower_bound = xbar[4] - t.score*sample.se
upper_bound = xbar[4] + t.score*sample.se
c(lower_bound, upper_bound)
#Therefore the C.I. for West would be (1542.455, 1861.445)

#Part C
#Sample Mean Vector
onemat = matrix(1, n, 1)
xbar =1/n*t(X)%*%onemat
xbar

#Sample Variance-Covariance Matrix
M = t(X)%*%X
L = xbar%*%t(xbar)
N = n*L
S = 1/(n - 1)*(M-N)
S

#Sample Correlation Matrix
variances = diag(S)
D = matrix(diag(variances),ncol=4)
D_sqrt = sqrt(D)
D_sqrt_inv = solve(D_sqrt)
samp_cor = D_sqrt_inv%*%S%*%D_sqrt_inv
samp_cor

#Can also use correlation function
cor(X)

#Part D
library(RVAideMemoire)
mqqnorm(X, main = 'Multi-normal Q-Q plot')

#Part E
S_inv = solve(S)
S_inv

#Part F
mu0 = matrix(c(1450,1900,1700,1700), ncol=1, byrow=T)
Tobs = n*t(xbar-mu0)%*%S_inv%*%(xbar-mu0)
Tobs

Fcriticalvalue = (n-1)*p/(n-p)*qf(p = 0.05, df1 = p, df2 = n-p, lower.tail = FALSE)
Fcriticalvalue

pvalue = pf((n-p)/((n-1)*p)*Tobs, p, n-p,lower.tail = FALSE)
pvalue


#Part G
#Based of the R code for Part E, since the p-value is greater than 0.05, we would say that the vector
#mu = (1450, 1900, 1700, 1700) would fall within the 95% confidence region. Furthermore, we can say that since the 
#Hotelling T^2 observed statistic is not greater than the F critical value, we cannot reject H_0 and say that there is
#no evidence to show that population mean vector is different from mu_0 = (1450, 1900, 1700, 1700). 


