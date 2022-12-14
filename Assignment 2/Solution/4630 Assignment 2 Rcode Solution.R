#Question 2
L1 = matrix(c(7.80, 90.4, 7.12, 85.1, 7.10, 88.9, 7.06, 89.0,
              7.89, 85.9, 7.45, 75.9, 7.82, 88.8, 7.45, 77.9), ncol=2, byrow=T)
L2 = matrix(c(9.00, 82.5, 8.19, 66.0, 8.43, 92.4, 8.25, 74.5,
              7.65, 82.4, 7.45, 83.1, 7.70, 87.4, 7.45, 86.4), ncol=2, byrow=T)
L3 = matrix(c(7.60, 94.1, 7.06, 81.2, 7.00, 86.6, 7.04, 79.9,
              7.82, 85.9, 7.52, 86.4, 7.80, 88.8, 7.70, 76.4), ncol=2, byrow=T)
E = matrix(c(3.8062, -0.1721, -0.1721, 895.6663), ncol=2, byrow=T)
Spooled = E/21

S1 = var(L1)
S2 = var(L2)
S3 = var(L3)

n1 = n2 = n3 = 8
p = 2
a = 3

lambda = (det(S1)/det(Spooled))^((n1-1)/2)*
  (det(S2)/det(Spooled))^((n2-1)/2)*
  (det(S3)/det(Spooled))^((n3-1)/2)

M = -2 * log(lambda)

U = ((1/(n1-1)+1/(n2-1)+1/(n3-1)) - 1/((n1-1)+(n2-1)+(n3-1)))*
  (2*p*p + 3*p - 1)/(6*(p+1)*(a-1))

obsteststat = (1 - U)*M

df = (1+p)*p*(a-1)/2

1-pchisq(obsteststat, df)





#Question 3
library(DescTools)
c1x1 = c(73, 43, 47, 53, 58, 47, 52, 38, 61, 56, 56, 34, 55, 65, 75)
c1x2 = c(31, 19, 22, 26, 36, 30, 29, 36, 34, 33, 19, 19, 26, 15, 18)
c2x1 = c(51, 41, 43, 41, 47, 32, 24, 43, 53, 52, 57, 44, 57, 40, 68)
c2x2 = c(35, 14, 19, 29, 34, 26, 19, 37, 24, 27, 14, 19, 30,  7, 13)
d = matrix(c(c1x1-c2x1, c1x2-c2x2), ncol=2, byrow=F)
mu0 = c(0, 0)
HotellingsT2Test(d, mu=mu0)


