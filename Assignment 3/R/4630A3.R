##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             4630: Assignment 3                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data sets ----------------------------------------------------------------
library(xlsx2dfs)

q1df = read.xlsx(file.choose(), sheet = 1)
q2df = read.xlsx(file.choose(), sheet = 2)
q3df = q2df
q5df = read.xlsx(file.choose(), sheet = 3)
q5df = q5df[,-1]



# Functions Built ---------------------------------------------------------
xbar = function(matrix){
  n = dim(matrix)[1]
  onevec = rep(1,n)
  xbar = 1/n*t(matrix)%*%onevec
  return(xbar)
}
Spool = function(dat){
  #Getting the n lengths and parameters
  nbar = rep(0,length(dat)) 
  for (i in 1: length(dat)){
    nbar[i] = dim(dat[[i]])[1]
  }
  n = sum(nbar)
  params = dim(dat[[1]])[2]
  g = length(dat)
  
  #Getting the variance matrices
  S = list()
  for (i in 1:g){
    S = append(S,list(var(dat[[i]])))
  }
  
  #Calculating the Within Matrix
  W = matrix(0,params,params)
  Wnew = matrix(0,params,params)
  for (i in 1:g){
    Wnew = as.matrix(lapply(S[i], '*', (nbar[i] -1)))
    Wnew = matrix(unlist(Wnew), ncol = params, byrow = T)
    W = W + Wnew
  }
  #s_pooled
  result = W/(n - g)
  return(result)
}
spectral = function(matrix){
  eig = eigen(matrix)
  D   = diag(eig$values)
  P   = eig$vectors
  matsqrt = P%*%sqrt(D)%*%t(P)
  output = list(D, P, matsqrt)
  names(output) = c('D', 'P', 'SqrtMat')
  return(output)
}

# Question 1 --------------------------------------------------------------

#Part A
#Check Notes

#Part B
Y = as.matrix(q1df[,1:2])
X = as.matrix(q1df[,3:4])
onevec = rep(1, dim(X)[1])
Xnew = cbind(onevec,X)
beta_coef = solve(crossprod(Xnew))%*%crossprod(Xnew, Y)

#OR 
fit = lm(Y ~ X)
summary(fit)


#Part C
Yhat = Xnew%*%beta_coef
ehat = Y - Yhat
n = dim(X)[1]
Sigmahat = 1/n*t(ehat)%*%ehat
E = n*Sigmahat

ybar = colMeans(Y)
n = nrow(Y)
m = ncol(Y)
r = ncol(X)
q = 1

Ybar = matrix(ybar, n ,m, byrow = T)
SST = crossprod(Y - Ybar)



#Part D
X1 = Xnew[,1:2]
Betahat1 = solve(crossprod(X1))%*%crossprod(X1,Y)
Sigmahat1 = 1/n*crossprod(Y - X1%*%Betahat1)

E = n*Sigmahat
H = n*(Sigmahat1 - Sigmahat)

Lambda = det(E)/det(E + H)


# Question 2 --------------------------------------------------------------
#Part A
X = data.matrix(q2df)
n = dim(X)[1]
params = dim(X)[2]
onevec = rep(1, n)

#Mean Vector
xbar = 1/n*t(X)%*%onevec
xbar

#Variance Matrix
S = 1/(n - 1)*(t(X)%*%X - n*xbar%*%t(xbar))
S

#OR
var(X)

#Getting the eigenvalues and vectors
eig = eigen(S)
eigvalues = eig$values

cumeig = rep(0,length(eigvalues))
j = 0
jnew = 0

for (i in 1: length(eig$values)){
  jnew = eigvalues[i]/sum(eigvalues)
  j = j + jnew
  cumeig[i] = j
}

cumeig


cumsum(eigvalues)/sum(eigvalues)

#PCA equations

#Y1 = -0.4728X1 - 0.3918X2 - 0.4875X3 - 0.4677X4 - 0.4080X5
#Y2 =  0.5763X1 + 0.1082X2 - 0.0960X3 + 0.1206X4 - 0.7952X5
#Y3 =  0.4168X1 - 0.4524X2 + 0.4794X3 - 0.6195X4 + 0.0887X5
#Y4 =  0.2286X1 + 0.6558X2 - 0.3690X3 - 0.5803X4 + 0.2115X5
#Y5 =  0.4672X1 - 0.4472X2 - 0.6222X3 + 0.2146X4 + 0.3854X5

#OR 

result = prcomp(X, scale = F)
print(result)

#Part B
var_exp = result$sdev^2/sum(result$sdev^2)

plot(var_exp, type = 'b', col = 'green', yaxt = 'n', ylab = '', xlab = '')
par(new = T)
plot(cumsum(var_exp), xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot Using Covariance Matrix",
     ylim = c(0, 1), type = "b", col = 'blue')
legend(4.2, 0.4, legend=c("Proportion", "Cumulative"),
       col=c("green", "blue"), lty=1:1, cex=0.8)

# I would use 3 principal components. Since Timm needs to cover at least 90%, we
#cannot have less than 3 PC. By adding the 4th or 5th PC, it barely adds to the 
#explanation of variance. 



#Part C
#Calculate Correlation Matrix
sd = sqrt(diag(diag(S),5))
invsd = solve(sd)
cor = invsd%*%S%*%t(invsd)
cor

cor(X)

eig = eigen(cor)
eigvalues = eig$values

#OR 

result = prcomp(scale(X), scale = F )
print(result)

cumeig = rep(0,length(eigvalues))
j = 0
jnew = 0

for (i in 1: length(eig$values)){
  jnew = eigvalues[i]/sum(eigvalues)
  j = j + jnew
  cumeig[i] = j
}
cumeig

var_exp = result$sdev^2/sum(result$sdev^2)
var_exp


plot(var_exp, type = 'b', col = 'green', yaxt = 'n', ylab = '', xlab = '')
par(new = T)
plot(cumsum(var_exp), xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot Using Correlation Matrix",
     ylim = c(0, 1), type = "b", col = 'blue')
legend(4.2, 0.4, legend=c("Proportion", "Cumulative"),
       col=c("green", "blue"), lty=1:1, cex=0.8)



# Question 3 --------------------------------------------------------------
X = data.matrix(q3df)
col_order = c('x4', 'x5', 'x1', 'x2', 'x3')
X = X[,col_order]

#Part A
cor(X)

# x4        x5        x1        x2        x3
# x4 1.0000000 0.5238876 0.5750730 0.7497770 0.6052716
# x5 0.5238876 1.0000000 0.4130573 0.5476595 0.6918927
# x1 0.5750730 0.4130573 1.0000000 0.6143902 0.7571850
# x2 0.7497770 0.5476595 0.6143902 1.0000000 0.5473897
# x3 0.6052716 0.6918927 0.7571850 0.5473897 1.0000000


#Covariance Matrices
cor11 = matrix(c(1.0000000, 0.5238876,
                  0.5238876, 1.0000000),
                ncol = 2, byrow = T)

cor11_sqrt = spectral(cor11)
cor11_invsqrt = solve(cor11_sqrt$SqrtMat)

cor12 = matrix(c(0.5750730, 0.7497770, 0.6052716,
                  0.4130573, 0.5476595, 0.6918927),
                ncol = 3, byrow = T)

cor21 = t(cor12)

cor22 = matrix(c(1.0000000, 0.6143902, 0.7571850,
                  0.6143902, 1.0000000, 0.5473897,
                  0.7571850, 0.5473897, 1.0000000),
                ncol = 3, byrow = T)

cor22_inv = solve(cor22)


A = cor11_invsqrt%*%cor12%*%cor22_inv%*%cor21%*%cor11_invsqrt
A

canon_cor = sqrt((eigen(A)$values))
canon_cor

#Part B
n = dim(X)[1]
p = min(3,2)
q = max(3,2)

Bartlett = -(n - 1 - 1/2*(p + q + 1))*log(prod((1-canon_cor^2)))
Bartlett

df = p*q

1-pchisq(Bartlett, df)

#P-value is quite large, hence we cannot reject H0.


#Part C
#Testing p_1 not equal to  0 
Bartlett = -(n - 1 - 1/2*(p + q + 1))*log((1-canon_cor[2]^2))
Bartlett

df = (p - 1)*(q-1)

1-pchisq(Bartlett, df)

#Testing p_1 and p_2 not equal to 0
Bartlett = -(n - 1 - 1/2*(p + q + 1))*log((1-canon_cor[2]^2))
Bartlett

df = (p - 2)*(q - 2)

1-pchisq(Bartlett, df)

#??????

# Question 4 --------------------------------------------------------------

#Part B.1
curve(expr = (1 - abs(x)), from = -1, to = 1, col = 'green', xlim = c(-1, 2),
      main = "Probability Density Function Graph",
      ylab = 'Density',
      xlab = 'Values of x')
curve(expr = 1 - abs(x - 1/2), from = -0.5, to = 1.5, col = 'blue', add = TRUE)



abline(h = 0.2)
abline(v = -1/3)

abline(h = 0.8, lty = 1, col='green')
abline(v = -1/3)

#Part B.2
x = seq(-1, 1, 0.01)

density1 = density(1 - abs(x))

x = seq(-0.5, 1.5, 0.01)

density2 = density(1 - abs(x - 1/2))


# Question 5 --------------------------------------------------------------
A = as.matrix(q5df[q5df$Group == "A",1:4])
B =as.matrix(q5df[q5df$Group =="B",1:4])
dat = list(A, B)

#Part A
# Mean Vectors
xbar1 = xbar(A)
xbar2 = xbar(B)
spooled = Spool(dat)

ahat = t(xbar1 - xbar2)%*%solve(spooled)

ybar1 = ahat%*%xbar1
ybar2 = ahat%*%xbar2

midpoint = 1/2*(ybar1 + ybar2)
midpoint

# Hence the classification is -7.062536

#Part B

cutoff = log(1/2 * (0.75/0.25))
cutoff

#Part C



#Part D

#When equal cost and prior discriminants
newobs = c(50, 48, 47, 49)
ahat%*%newobs

#Since the observation is bigger than the cutoff point, we will classify
#this obs. to population A. 

#When p1 = 0.25, and c(1|2) is half of c(2|1)

ahat%*%newobs - midpoint

#Since the cutoff value for this classification is 0.4055, we will also
#classify this obs. to population A. 


#When we use Bayesian Rule






