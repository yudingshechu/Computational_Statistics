# test the shrinkage path of ridge and Lasso
rm(list = ls())
library(glmnet)
library(MASS)
set.seed(777)
n <- 100
beta.true <- c(0.5, 0.5, -0.5)
sigma <- matrix(c(2,0.1,0.1,
                  0.1,2,0.1,
                  0.1,0.1,2) ,
                nrow= 3, ncol= 3, byrow=TRUE)

sigma2 <- matrix(c(2,0.1,0.1,
                  0.1,4,0.1,
                  0.1,0.1,5) ,
                nrow= 3, ncol= 3, byrow=TRUE)
sigma3 <- matrix(c(2,1.9,0.1,
                  1.9,2,0.1,
                  0.1,0.1,2) ,
                nrow= 3, ncol= 3, byrow=TRUE)
sigma4 <- matrix(c(2,0.1,1.9,
                  0.1,2,0.1,
                  1.9,0.1,2) ,
                nrow= 3, ncol= 3, byrow=TRUE)
sigma5 <- matrix(c(2,1.9,1.9,
                  1.9,2,1.9,
                  1.9,1.9,2) ,
                nrow= 3, ncol= 3, byrow=TRUE)
mu <- rep(0,3)
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

set.seed(777)
data <- data.generator(n, sigma, mu, beta.true)
set.seed(777)
data <- data.generator(n, sigma2, mu, beta.true)
set.seed(777)
data <- data.generator(n, sigma3, mu, beta.true)
set.seed(777)
data <- data.generator(n, sigma4, mu, beta.true)
set.seed(777)
data <- data.generator(n, sigma5, mu, beta.true)

x.sd <- cbind(data$x.sd.x.1, data$x.sd.x.2, data$x.sd.x.3)
x <- cbind(data$x.1, data$x.2, data$x.3)
ridge <- function(x,y,lambda, p = 3){
  beta.ridge <- solve(t(x) %*% x + 
                        diag(lambda, nrow = p, ncol = p)) %*% t(x) %*% y
  return(beta.ridge)
}
grid <- 10^ seq (3, -3, length = 100)
beta.ridge <- matrix(NA, nrow = 3, ncol= length(grid))
beta.ridge.nonsd <- matrix(NA, nrow = 3, ncol= length(grid))
for (i in 1:length(grid)){
  lam <- grid[i]
  beta.ridge[,i] <- ridge(x.sd, data$y, lam)
}
for (i in 1:length(grid)){
  lam <- grid[i]
  beta.ridge.nonsd[,i] <- ridge(x, data$y, lam)
}

### Lasso
lasso.mod<-glmnet(x.sd,data$y, alpha=1, lambda = grid,intercept=F)
lassobeta <- as.matrix(lasso.mod$beta)
plot(lasso.mod)

ylimits = c(-1.5, 1.5)
plot(x=log(grid), y=lassobeta[1,], col = "red", ylim = ylimits, 
     xlab = expression(Log~lambda), ylab = "", 
     main = expression(hat(beta) ~ "for different" ~ lambda), type = "n")
grid()
points(x=log(grid), y=lassobeta[1,], col = "orange", lwd = 1)
points(x=log(grid), y=lassobeta[2,], col = "red", lwd = 1)
points(x=log(grid), y=lassobeta[3,], col = "darkblue", lwd = 1)
abline(h = 0, col = "black")
legend("topright", c(expression(beta[1]), expression(beta[2]), expression(beta[3])), 
       col = c("orange", "red", "darkblue"), pch = 1)
mtext("LASSO")


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
mtext("Ridge")

ylimits = c(-1.5, 1.5)
plot(x=grid, y=beta.ridge.nonsd[1,], col = "red", ylim = ylimits, 
     xlab = expression(lambda), ylab = "", 
     main = expression(hat(beta) ~ "for different" ~ lambda), type = "n")
grid()
points(x=grid, y=beta.ridge.nonsd[1,], col = "orange", lwd = 1)
points(x=grid, y=beta.ridge.nonsd[2,], col = "red", lwd = 1)
points(x=grid, y=beta.ridge.nonsd[3,], col = "darkblue", lwd = 1)
abline(h = 0, col = "black")
legend("topright", c(expression(beta[1]), expression(beta[2]), expression(beta[3])), 
       col = c("orange", "red", "darkblue"), pch = 1)
mtext("non-standardized data")


qqnorm(data$x.sd.x.1)
qqline(data$x.sd.x.1)
qqnorm(data$x.sd.x.2)
qqline(data$x.sd.x.2)
qqnorm(data$x.sd.x.3)
qqline(data$x.sd.x.3)
hist(data$x.sd.x.1)
hist(data$x.sd.x.2)
hist(data$x.sd.x.3)
olscoe <- matrix(NA,1000,3)
for (i in 1:1000) {
  data <- data.generator(n, sigma, mu, beta.true)
  x <- cbind(data$x.1, data$x.2, data$x.3)
   x.sd <- cbind(data$x.sd.x.1, data$x.sd.x.2, data$x.sd.x.3)
  olscoe[i,] <- coef(lm(data$y~x.sd-1))
}
plot(1:1000,olscoe[,1]-olscoe[,2])
boxplot(olscoe[,1]-olscoe[,2])
matplot(olscoe)
boxplot(olscoe)
