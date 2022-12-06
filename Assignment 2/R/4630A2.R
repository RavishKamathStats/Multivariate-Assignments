##### 4630 Assignment 2 ####
install.packages('mvtnorm')
install.packages('pacman')
install.packages('Mass')
pacman::p_load("car", "psych", "RVAideMemoire")
library("car")
library("RVAideMemoire")
library("psych")
library(mvtnorm)
library(MASS)

#####DATA#####
df = read.csv(file.choose(), header =  TRUE)
print(df)
View(df)


#### Question 1 ####
#Part B
SSE = matrix(c(3.1532, -14.9535, -14.9535, 535.9725), nrow = 2, ncol = 2, byrow = T)
SSLub = matrix(c(1.6927, -9.6989, -9.6989, 56.3233), nrow = 2, ncol = 2, byrow = T)
SSVel = matrix(c(0.6240, 14.8834, 14.8834, 354.9704), nrow = 2, ncol = 2, byrow = T)
SSInt = matrix(c(0.0290, -0.102, -0.102, 4.7233), nrow = 2, ncol = 2, byrow = T)
SStr = SSLub + SSVel + SSInt
SST = SSE + SStr

wilkslam = det(SSE)/det(SSE + SStr)
wilkslam

a = 3
b = 2
p = 2
n = 4

#chi-squared observed
chi_obs = -(a*b*(n-1)- ((p+1) - (a*b - 1))/2)*log(wilkslam)
chi_obs

#P-Value
pchisq(chi_obs, df = 10, lower.tail = F)

#### Question 2 ####
# Part A
n1_indices = which(df$Lubricant == 'L1')
n1 =  dim(df[n1_indices,])[1]
n2_indices = which(df$Lubricant == 'L2')
n2 =  dim(df[n2_indices,])[1]
n3_indices = which(df$Lubricant == 'L3')
n3 =  dim(df[n3_indices,])[1]
grp_obs = c(n1, n2 ,n3)
n = sum(grp_obs)
g = nlevels(factor(df$Lubricant))
p = 2

#Finding the means for each group and overall mean
L1mean = as.vector(colMeans(df[n1_indices[1]:tail(n1_indices, n = 1),3:4]))
L2mean = as.vector(colMeans(df[n2_indices[1]:tail(n2_indices, n = 1),3:4]))
L3mean = as.vector(colMeans(df[n3_indices[1]:tail(n3_indices, n = 1),3:4]))
lub_group_mean = cbind(L1mean, L2mean, L3mean)
mean = (n1*L1mean + n2*L2mean + n3*L3mean)/(n1 + n2 + n3)
mean

#Finding the sample covariance for each lubricant
X1 = as.matrix(df[n1_indices[1]:tail(n1_indices, n = 1),3:4])
X2 = as.matrix(df[n2_indices[1]:tail(n2_indices, n = 1),3:4])
X3 = as.matrix(df[n3_indices[1]:tail(n3_indices, n = 1),3:4])
S1 = 1/(n1 - 1)*(t(X1)%*%X1 - n1*colMeans(X1)%*%t(colMeans(X1)))
S2 = 1/(n2 - 1)*(t(X2)%*%X2 - n2*colMeans(X2)%*%t(colMeans(X2)))
S3 = 1/(n3 - 1)*(t(X3)%*%X3 - n3*colMeans(X3)%*%t(colMeans(X3)))
S = list(S1, S2, S3)

#Calculate the Within Matrix
W = matrix(0,p,p)
Wnew = matrix(0,p,p)
for (i in 1:g){
  Wnew = as.matrix(lapply(S[i], '*', (grp_obs[i] -1)))
  Wnew = matrix(unlist(Wnew), ncol = 2, byrow = T)
  W = W + Wnew
}
print(W)

#Calculating the Between Matrix
B = matrix(0,p,p)
Bnew = matrix(0,p,p)
for (i in 1:g){
  Bnew = grp_obs[i]*(lub_group_mean[,i] - mean)%*%t(lub_group_mean[,i] - mean)
  B = B + Bnew
}
print(B)

#Wilks Lambda and finding p-value
wilkslam = det(W)/det(W + B)
wilkslam
fobs = (n - g - 1)/(g - 1)*((1 - sqrt(wilkslam))/sqrt(wilkslam))
fobs
pvalue = pf(fobs, df1 = 2*(g-1),df2 = 2*(n - g - 1), lower.tail = F)
pvalue

#P-value is 0.07. If we have our alpha to be 0.05, then there is no evidence to reject H_0
#Hence this implies that lubricant effect does not exist. 

# Part B
#two sample difference
tau_vec = lub_group_mean
A = c(1, 0, -1)
tau_hat = tau_vec%*%A
tau_hat

#Finding the pooled variance for L1 and L3 
Spooled_L1andL3 = ( (n1 - 1)*S1 +(n3 - 1)*S3 )/(n1 + n3 - 2 )
Spooled_L1andL3

#Variance for tau_hat
var_tau_hat = (1/n1 + 1/n3)*Spooled_L1andL3
var_tau_hat

#95% Ellipsoid
plot(tau_hat[1], tau_hat[2] , type="p", xlim=c(-0.2, 0.2), 
     ylim=c(-2, 2), xlab="L1mean", ylab="L3mean", 
     main = '95% C.I. for mean difference between L1 and L3')

tau1 = matrix(seq(-2, 2, 0.005), ncol=1, byrow=T)
ntau1 = nrow(tau1)
tau3 = matrix(seq(-2, 2, 0.005), ncol=1, byrow=T)
ntau3 = nrow(tau3)

for (i in 1:ntau1) {
  for (j in 1:ntau3) {
    tau = matrix(c(tau1[i, 1], tau3[j, 1]), ncol=1, byrow=T)
    Tsq_obs = t((tau_hat - tau))%*%solve(var_tau_hat)%*%(tau_hat-tau)
    
    Fcomp = c( ( (n1 + n3 - 2) - p + 1)/((n1 + n3 - 2)*p )* Tsq_obs)
    Fcrit = qf(0.05, p, n1 + n3 - p - 1)
    if (Fcomp < Fcrit) points(tau1[i, 1], tau3[j, 1], pch="*")
  }
}
points(tau_hat[1], tau_hat[2], col='red')

#IGNORE CODE
#t critical value
# alpha = 0.05
# tcrit = qt(p = alpha/(p*g*(g - 1)), df = (n - g), lower.tail=F)
# tcrit
# 
# #C.I.
# lower = rep(0,2)
# upper = rep(0,2)
# for (i in 1 : p ){
# lower[i] = tau_hat[i] - tcrit*sqrt(W[i,i]/(n-g)*(1/n1 + 1/n3))  
# upper[i] = tau_hat[i] + tcrit*sqrt(W[i,i]/(n-g)*(1/n1 + 1/n3))
# }
# 
# CR = cbind(lower, upper)
# CR


# Part C
#Taking the determinants of the sample co variances
detS1 = det(S1)
detS2 = det(S2)
detS3 = det(S3)
detS = c(detS1, detS2, detS3)

#Finding the Spooled and its determinant
Spooled = W/(n1 + n2 + n3 - g)
detSpooled = det(Spooled)
detSpooled

#Doing the Box Test for equality of co variances
grp_obs = c(n1, n2 ,n3)
lambda = rep(1,1)
lambda_new = rep(0,1)
for (i in 1:g){
  lambda_new = (detS[i]/detSpooled)^((grp_obs[i]-1)/2)
  lambda = lambda*lambda_new
}
print(lambda)

M = -2*log(lambda)
M
u = (sum(1/(grp_obs - 1)) - 1/(sum(grp_obs - 1)))*((2*p^2+3*p-1)/(6*(p+1)*(g-1)))
u
C = (1 - u)*M
pchisq(C, df = ((1+p)*p*(g-1))/2, lower.tail = F)

#### Question 3 ####
df = read.csv(file.choose(), header = T)
print(df)

n1_indices = which(df$Coating == 1)
n1 =  dim(df[indices,])[1]
n2_indices = which(df$Coating == 2)
n2 =  dim(df[indices,])[1]
grp_obs = c(n1, n2)
n = grp_obs[1]
p = 2
g = nlevels(factor(df$Coating))

#Difference of Mean
C1 = df[n1_indices, 2:3]
C2 = df[n2_indices, 2:3]
d = C1 - C2
d_bar = colMeans(d)
d_bar

# sample Variance covariance matrix
S = cov(d)
S

HotellingT = n*t(d_bar)%*%solve(S)%*%d_bar
HotellingT
F_Obs = (n-p)/((n-1)*p)*HotellingT
F_Obs
pvalue = pf(F_Obs, df1 = p, df2 = 15 - p, lower.tail = F)
pvalue

##### Wong's Example #####



df2 = matrix(c(1,9,3,1,6,2,1,9,7,2,0,4,2,2,0,3,3,8,3,1,9,3,2,7),ncol =3, byrow = T)
View(df2)
X1 = df2[1:3,2:3]
X2 = df2[4:5,2:3]
X3 = df2[6:8,2:3]

n1 = length(df2[1:3,1])
n2 = length(df2[4:5,1])
n3 = length(df2[6:8,1])
n = n1 + n2 + n3 
p = 2

S1 = 1/(n1 - 1)*(t(X1)%*%X1 - n1*colMeans(X1)%*%t(colMeans(X1)))
S2 = 1/(n2 - 1)*(t(X2)%*%X2 - n2*colMeans(X2)%*%t(colMeans(X2)))
S3 = 1/(n3 - 1)*(t(X3)%*%X3 - n3*colMeans(X3)%*%t(colMeans(X3)))

W = (n1 - 1)*S1 + (n2-1)*S2 + (n3 - 1)*S3





df1 = read.csv('/Users/ravishkamath/Desktop/University/2. York Math/1 MATH/1. Statistics /MATH 4630/1. Lectures/4. Comparisons of Several Multivariate Means/SAS Code/Ex6.13_PlasticFlim.csv', header = TRUE)
View(df1)




print(df1)
dep_var = cbind(df1$x1, df1$x2, df1$x3)
df1$A = as.factor(df1$A)
df1$B = as.factor(df1$B)

model = lm(dep_var ~ A + B + A*B - 1, data = df1)
manova = Manova(model)
summary(manova)

sse = c(1.76, -0.02, -3.07,
        -0.02, 2.63, -0.55, 
        -3.07, -0.55, 64.92)
ssint = c(0.005, 0.02, 0.04, 
          0.02, 0.54, 1.47, 
          0.04, 1.47, 3.96)

ss1 = c(1.74, -1.5, 0.86,
        -1.5, 1.3, -0.74,
        1.93, -0.74, 0.42)
ss2 = c(0.76, 0.68, 1.93,
        0.68, 0.61, 1.73,
        1.94, 1.73, 4.9)
SSE = matrix(data = sse, nrow = 3, ncol = 3, byrow = TRUE)
SSint = matrix(data = ssint, nrow = 3, ncol = 3, byrow = TRUE)
SS1 = matrix(data = ss1, nrow = 3, ncol = 3, byrow = TRUE)
SS2 = matrix(data = ss2, nrow = 3, ncol = 3, byrow = TRUE)


SStr = SSint + SS1 + SS2

wilksobs = det(SSE)/det(SSE + SStr)





