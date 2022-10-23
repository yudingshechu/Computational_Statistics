################# Exercise 1 #######################
################### a, training sample
set.seed(55)
N <- 1000
X1 <- rep(1,N)
X2 <- rnorm(N,0,sqrt(1.5))
e <- rnorm(N,0,sqrt(10))
beta <- c(5,-0.5)
X <- cbind(X1,X2)
Y <- X%*%beta + e
training <- data.frame(X,Y)

####################### b, test sample
X10 <- rep(1,N)
X20 <- rnorm(N,0,sqrt(1.5))
e0 <- rnorm(N,0,sqrt(10))
X0 <- cbind(X10,X20)
Y0 <- X0%*%beta + e0
testdata <- data.frame(X0,Y0)

############ c, beta_hat is calculated by training sample
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
beta_hat

##################### d, 
## training MSE
fX_hat <- X%*%beta_hat ## predicted by training data 
deviation <- Y-fX_hat
MSE <- sum(deviation^2)/N
MSE
## Average prediction error
fX_hat0 <- X0%*%beta_hat # predicted by test data
deviation <- Y0-fX_hat0
PE <- mean(deviation^2)
PE

########################  e, 
## at first calculate polynomial for X2 and X20
X2_2 <- X2^2
X2_3 <- X2^3
X2_4 <- X2^4
X20_2 <- X20^2
X20_3 <- X20^3
X20_4 <- X20^4
Polyno <- data.frame(X1,X2,X2_2,X2_3,X2_4) #training data
Polyno0 <- data.frame(X1,X20,X20_2,X20_3,X20_4)#test data

X_loop <- c()

#because if we include X2^0 and X1 into function, there
# is multicolinearity, thus for the case of X2^0, I
# choose to drop constant term X1
X_loop0 <- c()
targets <- matrix(0,5,2)
beta_hats <- c()

# this loop first computes new beta_hat vector, then 
# use this to compute MSE and avg.prediction error
for (i in 1:5) {
  X_loop <- cbind(X_loop,Polyno[,i])
  beta_hats <- solve(t(X_loop)%*%X_loop)%*%t(X_loop)%*%Y
  
  fX_hat <- X_loop%*%beta_hats
  deviation <- Y-fX_hat
  targets[i,1] <- sum(deviation^2)/N
  
  X_loop0 <- cbind(X_loop0,Polyno0[,i])
  fX_hat0 <- X_loop0%*%beta_hats 
  deviation0 <- Y0-fX_hat0
  targets[i,2] <-  mean(deviation0^2)
}
colnames(targets) <- c("MSE","avg.P.E")
rownames(targets) <- c("X^0","X^1","X^2","X^3","X^4")
targets
################# exercise 2  #####################
##################### a & b, 
set.seed(100)
N <- 1000
T <- 1000
result <- matrix(NA,T,2)
X1 <- rep(1,N)
X10 <- rep(1,N)
beta <- c(5,-0.5)
# loop for 1000 times simulation creating 1000 MSE and PE
for (i in 1:T) {
X2 <- rnorm(N,0,sqrt(1.5))
e <- rnorm(N,0,sqrt(10))
X <- cbind(X1,X2)
Y <- X%*%beta + e
X20 <- rnorm(N,0,sqrt(1.5))
e0 <- rnorm(N,0,sqrt(10))
X0 <- cbind(X10,X20)
Y0 <- X0%*%beta + e0
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
fX_hat <- X%*%beta_hat 
deviation <- Y-fX_hat
result[i,1] <- sum(deviation^2)/N
fX_hat0 <- X0%*%beta_hat 
deviation <- Y0-fX_hat0
result[i,2] <- mean(deviation^2)
}
colnames(result) <- c("MSE","avg.P.E")
## the answer of b
avg <- c(mean(result[,1]),mean(result[,2]))
avg

######################## c, 
###################################################
set.seed(100)
N <- 1000
T <- 1000
targets_MSE <- matrix(NA,T,5)
targets_PE <- matrix(NA,T,5)
X1 <- rep(1,N)
X10 <- rep(1,N)
beta <- c(5,-0.5)

for (i in 1:T) {
  X2 <- rnorm(N,0,sqrt(1.5))
  e <- rnorm(N,0,sqrt(10))
  X <- cbind(X1,X2)
  Y <- X%*%beta + e
  X20 <- rnorm(N,0,sqrt(1.5))
  e0 <- rnorm(N,0,sqrt(10))
  X0 <- cbind(X10,X20)
  Y0 <- X0%*%beta + e0
  X2_2 <- X2^2
  X2_3 <- X2^3
  X2_4 <- X2^4
  X20_2 <- X20^2
  X20_3 <- X20^3
  X20_4 <- X20^4
  Polyno <- data.frame(X1,X2,X2_2,X2_3,X2_4) 
  Polyno0 <- data.frame(X1,X20,X20_2,X20_3,X20_4)
  X_loop <- c()
  X_loop0 <- c()
  beta_hats <- c()
  for (z in 1:5) {
    X_loop <- cbind(X_loop,Polyno[,z])
    beta_hats <- solve(t(X_loop)%*%X_loop)%*%t(X_loop)%*%Y
    
    fX_hat <- X_loop%*%beta_hats
    deviation <- Y-fX_hat
    targets_MSE[i,z] <- sum(deviation^2)/N
    
    X_loop0 <- cbind(X_loop0,Polyno0[,z])
    fX_hat0 <- X_loop0%*%beta_hats 
    deviation0 <- Y0-fX_hat0
    targets_PE[i,z] <-  mean(deviation0^2)
  }}
colnames(result) <- c("MSE","avg.P.E")
colnames(targets_MSE) <- c("X^0","X^1","X^2","X^3","X^4")
colnames(targets_PE) <- c("X^0","X^1","X^2","X^3","X^4")
result_final <- matrix(0,5,2)
rownames(result_final) <- c("X^0",
                            "X^1","X^2","X^3","X^4")
colnames(result_final) <- c("MSE","avg.P.E")
for (i in 1:5) {
  result_final[i,1] <- mean(targets_MSE[,i])
  result_final[i,2] <- mean(targets_PE[,i])
}
result_final
plot(c(0:4),result_final[,1],pch=18,col='red',
     main='average MSE',xlab = 'power over X',
     ylab = 'average MSE')
plot(c(0:4),result_final[,2],col='blue',pch=19,
     main='average Prediction Error',
     xlab = 'power over X',
     ylab = 'average Prediction Error',
     ylim = c(9.9,10.5))
points(c(0:4),result_final[,1],pch=18,col='red')

# the initial setting of std in error term and sample size
# influences the deviation between MSE and PE 
# the fewer N, the higher overfitting
# the fewer e std, the lower overfitting
# lower x std, more precise estimating
# lower beta means lower explained varity of Y

# to show the correct objective of function, we need to
# compare different e term under different std to show it
# (by some standardized scaling e)
# in plot showing, we should keep everything in the same 
# scale 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


