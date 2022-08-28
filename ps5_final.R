
rm(list = ls())
library(MASS)

# Parameters
set.seed(777)
n <- 100
beta.true <- c(0.5, 0.5, -0.5)
sigma <- matrix(c(2,0.1,0.1,
                  0.1,3,0.1,
                  0.1,0.1,4) ,
                nrow= 3, ncol= 3, byrow=TRUE)
mu <- rep(0,3)

# Data generator

data.generator <- function(n, sigma, mu, beta){
  x <- mvrnorm(n, mu, sigma)
  e <- rnorm(n, 0, sqrt(10))
  y <- x %*% beta + e
  x.1 <- x[,1] / sd(x[,1])
  x.2 <- x[,2]/ sd(x[,2])
  x.3 <- x[,3]/ sd(x[,3])
  x.sd <- cbind(x.1, x.2, x.3)
  data <- data.frame ("y" = y, "x.sd" = x.sd, "x" = x)
  return(data)
}

# Training data for the excerise

data <- data.generator(n, sigma, mu, beta.true)
# Storing the centered and uncentered covariates in a matrix for later convenience 

x.sd <- cbind(data$x.sd.x.1, data$x.sd.x.2, data$x.sd.x.3)
x <- cbind(data$x.1, data$x.2, data$x.3)

# a)

# Ridge-regression function

ridge <- function(x,y,lambda, p = 3){
  beta.ridge <- solve(t(x) %*% x + 
  diag(lambda, nrow = p, ncol = p)) %*% t(x) %*% y
  return(beta.ridge)
}

# b) 

# Grid for the lambdas 

grid <- 10^ seq (3, -3, length = 100)

# Storage matrix for the regression outputs 

beta.ridge <- matrix(NA, nrow = 3, ncol= length(grid))

# Ridge regression over the grid

for (i in 1:length(grid)){
  lam <- grid[i]
  beta.ridge[,i] <- ridge(x.sd, data$y, lam)
}

# Plot


ylimits = c(-1.5, 1.5)
plot(x=grid, y=beta.ridge[1,], col = "red", ylim = ylimits, 
     xlab = expression(lambda), ylab = "", 
     main = expression(hat(beta) ~ "for different" ~ lambda), type = "n")
grid()
points(x=grid, y=beta.ridge[1,], col = "orange", lwd = 1)
points(x=grid, y=beta.ridge[2,], col = "red", lwd = 1)
points(x=grid, y=beta.ridge[3,], col = "darkblue", lwd = 1)
abline(h = 0, col = "black")
legend("topright", c(expression(beta[1]), expression(beta[2]), expression(beta[3])), 
       col = c("orange", "red", "darkblue"), pch = 1)

# c)

set.seed(123)

# Creating a new grid for the values of lambda

grid <- 10^ seq (2, -10, length = 100)

# Generating test data

test.data <- data.generator(n, sigma, mu, beta.true)

# Uncentered test data
x.t <- cbind(test.data$x.1, test.data$x.2, test.data$x.3)

# Centered test data

x.sd.t <- cbind(test.data$x.sd.x.1, test.data$x.sd.x.2, test.data$x.sd.x.3)

# Containers for the train- and test errors of the ridge regression 

train.error <- c()
test.error <- c()

for (i in 1:length(grid)){
  y.hat <- x.sd %*% beta.ridge[,i]
  train.error[i] <- mean((data$y - y.hat)^2)
  y.hat.t <- x.sd.t %*% beta.ridge[,i]
  test.error[i] <- mean((test.data$y - y.hat.t)^2)
}

# Train - and test errors for the OLS function

lm.obj <- lm(data$y ~ x.sd.x.1 + x.sd.x.2 + x.sd.x.3 -1, data = data)
lm.obj.test <- x.sd.t %*% lm.obj$coefficients
#lm.obj.test <- predict(lm.obj, newdata= data.frame(x.sd.t))

ols.train.error <- mean((lm.obj$fitted.values - data$y)^2)
ols.test.error <- mean((lm.obj.test - test.data$y)^2)



# Plot

plot(x = log(grid), y = train.error, col = "red", 
     ylim = c(min(train.error, test.error), max(train.error, test.error)), 
     type = "n", 
     xlab = expression(log(lambda)), ylab = "Error Rate")
grid()
points(x = log(grid),  y = test.error, col = "blue")
points(x = log(grid),  y = train.error, col = "red")
abline(h = ols.train.error, col = "darkgreen", lwd = 3)
abline(h = ols.test.error, col = "orange", lwd = 3)
legend("topleft", legend = c("Ridge Training", "Ridge Test", "OLS Training", 
                             "OLS Test"), pch = c(1,1, NA, NA), 
       lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3),
       col = c("red", "blue", "darkgreen", "orange"), 
       cex = 1)


#d)

x.sd2 <- x.sd
x.sd2[,2] <- rep(-10, n)

grid <- 10^ seq (3, -3, length = 100)

# Storage matrix for the regression outputs 

beta.ridge.d <- matrix(NA, nrow = 3, ncol= length(grid))

# Ridge regression over the grid

for (i in 1:length(grid)){
  lam <- grid[i]
  beta.ridge.d[,i] <- ridge(x.sd2, data$y, lam)
}

# OLS regression


lm.obj <- lm(data$y ~ x.sd2 -1)

par(las = 1)
ylimits = c(-1.5, 1.5)
plot(x=grid, y=beta.ridge.d[1,], col = "red", ylim = ylimits, 
     xlab = expression(lambda), ylab = "", 
     main = expression(hat(beta) ~ "for different" ~ lambda), type = "n")
grid()
points(x=grid, y=beta.ridge.d[1,], col = "orange", lwd = 1)
points(x=grid, y=beta.ridge.d[2,], col = "red", lwd = 1)
points(x=grid, y=beta.ridge.d[3,], col = "darkblue", lwd = 1)
points(x=grid, y=beta.ridge[1,], col = "coral", lwd = 0.5, pch = 20)
points(x=grid, y=beta.ridge[2,], col = "firebrick2", lwd = 0.5, pch = 20)
points(x=grid, y=beta.ridge[3,], col = "blue", lwd = 0.5, pch= 20)
abline(h = 0, col = "black")
legend("bottomright", c(expression(beta[1]), expression(beta[2]), expression(beta[3])), 
       col = c("orange", "red", "darkblue"), pch = 1, title =" With constant")
legend("topright", c(expression(beta[1]), expression(beta[2]), expression(beta[3])), 
       col = c("coral", "firebrick2", "blue"), pch = 20, title = "Without constant")


# e)

library(glmnet)

x.sd <- cbind(data$x.sd.x.1, data$x.sd.x.2, data$x.sd.x.3)

ridge.mod <- glmnet(x.sd, data$y, alpha = 0, lambda = grid, tresh = 1e-12, intercept = FALSE)

set.seed(1)

cv.ridge <- cv.glmnet(x.sd, data$y, alpha = 0, intercept = FALSE)
plot(cv.ridge)

lam.star <- cv.ridge$lambda.min
lam.star


# Excerise 2

# a)

set.seed(527)

rep <- 100


MSE <- matrix(NA,rep,2)
colnames(MSE) <- c("Ridge", "OLS")

for (i in 1:rep){
  
  # Data generating process 
  train.data <- data.generator(n, sigma, mu, beta.true)
  test.data <- data.generator(n,sigma, mu, beta.true)
  
  # Storing the x variables into matrixes for convinience 
  
  # Train covariates
  x <- cbind(train.data$x.1, train.data$x.2, train.data$x.3)
  x.sd <- cbind(train.data$x.sd.x.1, train.data$x.sd.x.2, train.data$x.sd.x.3)
  
  #Test covariates
  
  x.t <- cbind(test.data$x.1, test.data$x.2, test.data$x.3)
  x.sd.t <- cbind(test.data$x.sd.x.1, test.data$x.sd.x.2, test.data$x.sd.x.3)
  
  # Ridge regession and prediction of the optimal lambda 
  
  ridge.mod <- glmnet(x.sd, train.data$y, alpha = 0,
                      lambda = grid, tresh = 1e-12, intercept = FALSE)
  
  cv.ridge <- cv.glmnet(x.sd, train.data$y, alpha = 0, intercept = FALSE)
  
  lam.star <- cv.ridge$lambda.min
  
  ridge.test <- predict(ridge.mod, s= lam.star, newx = x.sd.t)
  
  # Ridge MSE
  
  MSE[i,1] <- mean(( ridge.test - test.data$y)^2)
  
  # OLS fit 
  
  lm.obj <- lm(data$y ~ x.sd.x.1 + x.sd.x.2 + x.sd.x.3 -1, data = data)
  
  lm.obj.test <- x.sd.t %*% lm.obj$coefficients
  
  #OLS MSE
  
  MSE[i,2] <- mean( (lm.obj.test - test.data$y)^2)
}

# Plot OLS MSE vs. Ridge MSE

# Percentage of values above the 45 degree line

above <- mean(MSE[,1]<=MSE[,2])*100


plot(MSE[,1],MSE[,2], xlim = c(5,20), ylim = c(5,20),
     main = "Test MSE According to Ridge vs OLS",
     xlab = "Ridge", ylab = "OLS", type = "n")
grid()
points(MSE[,1],MSE[,2], col = "darkblue")
abline(a = 0, b = 1, col = "red")
text(10, 15, above)
text(10, 6, 100-above)



# Boxplot

boxplot(MSE, ylab = "MSE", col = c("cyan1", "lightsalmon"),  main = "Test MSE According to Ridge vs OLS")


#b )

##### change to a highly correlated covmatrix

set.seed(527)
n <- 100

# n <- 1000 #given a bigger n, with highly correlated vcov
# ridge performs well better 

beta.true <- c(0.5, 0.5, -0.5)

# a,b,c is the diagonal 

covmatrix <- function(a,b,c,d){
  matrix(c(a,sqrt(a*b)-d,sqrt(a*c)-3*d,
           sqrt(a*b)-d,b,sqrt(b*c)-d,
           sqrt(a*c)-3*d,sqrt(b*c)-d,c),
         nrow= 3, ncol= 3, byrow=TRUE)
}

sigma2 <- covmatrix(2,3,4,0.01)
sigma2 <- covmatrix(11,15,11,0.01)#works

mu <- rep(0,3)
rep <- 100
grid <- 10^ seq (5, -2, length = 100)

MSE2 <- matrix(NA,rep,2)
colnames(MSE2) <- c("Ridge*", "OLS*")

for (i in 1:rep){
  train.data <- data.generator(n, sigma2, mu, beta.true)
  test.data <- data.generator(n,sigma2, mu, beta.true)
  
  x <- cbind(train.data$x.1, train.data$x.2, train.data$x.3)
  x.sd <- cbind(train.data$x.sd.x.1, train.data$x.sd.x.2, train.data$x.sd.x.3)
  
  x.t <- cbind(test.data$x.1, test.data$x.2, test.data$x.3)
  x.sd.t <- cbind(test.data$x.sd.x.1, test.data$x.sd.x.2, test.data$x.sd.x.3)
  
  ridge.mod <- glmnet(x.sd, train.data$y, alpha = 0,
                      lambda = grid, tresh = 1e-12,intercept = F
  )
  cv.ridge <- cv.glmnet(x.sd, train.data$y, alpha = 0,intercept = F)
  
  lam.star <- cv.ridge$lambda.min
  
  ridge.test <- predict(ridge.mod, s= lam.star, newx = x.sd.t)
  
  lm.obj <- lm(train.data$y ~ train.data$x.sd.x.1 + train.data$x.sd.x.2 +
                 train.data$x.sd.x.3 -1)
  
  lm.obj.test <- x.sd.t %*% lm.obj$coefficients
  
  MSE2[i,1] <- mean(( ridge.test - test.data$y)^2)
  MSE2[i,2] <- mean( (lm.obj.test - test.data$y)^2)
}

# Plot OLS MSE vs. Ridge MSE


above2 <- mean(MSE2[,1]<=MSE2[,2])*100

plot(MSE2[,1],MSE2[,2], col="blue",xlim = c(8.5,11.5),ylim = c(8.5,11.5),
     xlab = "Ridge",ylab = "OLS", main = "Test MSE According to Ridge vs OLS (high covariance)")
lines(c(0:20),c(0:20), col = "red")
text(9, 10, above2)
text(10, 9.5, 100-above2)


#plot(1:100,MSE[,1]-MSE[,2])
#points(1:100,MSE2[,1]-MSE2[,2],col="red")
#lines(1:100,rep(0,100),col="blue")

# Boxplot

boxplot(cbind(MSE,MSE2),  ylab = "MSE", col = c("cyan1", "lightsalmon"), main = "Test MSE According to Ridge vs OLS, various scenarios" )
lines(0:7,rep(10,8),col="red", lwd= 2.5)


#######################################################
# now assume we have some noise variables
set.seed(527)
n <- 100
beta.true <- c(0.5, 0.5, -0.5)
mu <- rep(0,3)
rep <- 100
grid <- 10^ seq (5, -2, length = 100)

MSE3 <- matrix(NA,rep,2)

colnames(MSE3) <- c("Ridge?", "OLS?")

sigma <- matrix(c(2,0.1,0.1,0.1,3,0.1,0.1,0.1,4) ,
                nrow= 3, ncol= 3, byrow=TRUE)

for (i in 1:rep){
  
  train.data <- data.generator(n, sigma, mu, beta.true)
  test.data <- data.generator(n,sigma, mu, beta.true)
  
  # Noise variables are already centered 
  x4 <- rnorm(n,2,1)
  x5 <- runif(n,1-sqrt(3),1+sqrt(3))
  x6 <- runif(n,2-sqrt(3),2+sqrt(3))
  x7 <- rnorm(n,4,1)
  x8 <- runif(n,-1-sqrt(3),-1+sqrt(3))
  x9 <- rnorm(n,-2,1)
  x10 <- runif(n,-sqrt(3),sqrt(3))
  x4t <- rnorm(n,2,1)
  x5t <- runif(n,1-sqrt(3),1+sqrt(3))
  x6t <- runif(n,2-sqrt(3),2+sqrt(3))
  x7t <- rnorm(n,4,1)
  x8t <- runif(n,-1-sqrt(3),-1+sqrt(3))
  x9t <- rnorm(n,-2,1)
  x10t <- runif(n,-sqrt(3),sqrt(3))
  
  
  x.sd <- cbind(train.data$x.sd.x.1, train.data$x.sd.x.2,
                train.data$x.sd.x.3,x4,x5,x6,x7,x8,x9,x10)
  x.sd.t <- cbind(test.data$x.sd.x.1, test.data$x.sd.x.2, 
                  test.data$x.sd.x.3,x4t,x5t,x6t,x7t,x8t,x9t,x10t)
  
  
  ridge.mod <- glmnet(x.sd, train.data$y, alpha = 0,
                      lambda = grid, tresh = 1e-12,intercept = F)
  
  cv.ridge <- cv.glmnet(x.sd, train.data$y, alpha = 0,intercept = F)
  
  lam.star <- cv.ridge$lambda.min
  
  ridge.test <- predict(ridge.mod, s= lam.star, newx = x.sd.t)
  
  lm.obj <- lm(train.data$y ~ train.data$x.sd.x.1 + train.data$x.sd.x.2 +
                 train.data$x.sd.x.3+x4+x5+x6+x7+x8+x9+x10 -1)
  lm.obj.test <- x.sd.t %*% lm.obj$coefficients
  
  MSE3[i,1] <- mean(( ridge.test - test.data$y)^2)
  MSE3[i,2] <- mean( (lm.obj.test - test.data$y)^2)
}

# Plot OLS MSE vs. Ridge MSE

above3 <- mean(MSE3[,1]<=MSE3[,2])*100

plot(MSE3[,1],MSE3[,2],col="blue",
     xlab = "Ridge",ylab = "OLS",xlim = c(7,15),
     ylim = c(7,15), main = "Test MSE According to Ridge vs OLS (noise variables)")
lines(c(0:20),c(0:20), col = "red")
text(10,12, above2)
text(10,8, 100-above2)



#plot(1:100,MSE[,1]-MSE[,2])
#points(1:100,MSE3[,1]-MSE3[,2],col="red")
#lines(1:100,rep(0,100),col="blue")

# Boxplot
boxplot(cbind(MSE,MSE2,MSE3), ylab = "MSE", col = c("cyan1", "lightsalmon"), main = "Test MSE According to Ridge vs OLS, various scenarios" )
lines(0:7,rep(10,8),col="red", lwd= 2.5)


plot(MSE[,1],MSE[,2],col="darkblue",
     xlab = "Ridge",ylab = "OLS",xlim = c(7,15),
     ylim = c(7,15), main = "Test MSE According to Ridge vs OLS, various scenarios", pch= 16)
grid()
points(MSE2[,1],MSE2[,2],col="orange", pch = 15)
points(MSE3[,1],MSE3[,2],col="red", pch = 17)
lines(c(0:20),c(0:20), col = "black")
legend("bottomright", c("Basline", "High covariance", "Noise variables"),
       col = c("darkblue", "orange", "red"), pch = c(16,15,17))




# Reduce n to below p

n_low <- 2
set.seed(55) # same as above 

data_n_low <- data.generator(n = n_low, sigma, mu, beta = beta.true)
x.sd_n_low <- cbind(data_n_low$x.sd.x.1, data_n_low$x.sd.x.2, data_n_low$x.sd.x.3)
x_n_low <- cbind(data_n_low$x.1, data_n_low$x.2, data_n_low$x.3)

# OLS fails: 

try(beta.OLS_n_low <- solve(t(x_n_low) %% x_n_low) %% t(x_n_low) %*% y)

# Ridge works: 

beta.ridge_n_low <- matrix(NA, ncol = length(grid), nrow = 3)

for (i in 1:length(grid)){
  lam <- grid[i]
  beta.ridge_n_low[,i] <- ridge(x.sd_n_low, data_n_low$y, lam)
}

head(t(beta.ridge_n_low))



