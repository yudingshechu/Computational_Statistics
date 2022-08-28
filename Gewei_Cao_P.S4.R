#install.packages("boot")
library(boot)
rm(list=ls())
##############################################
### 1.a
n <- 1000
beta <- matrix(c(1,1.5,-1.5,1.5,0.5),c(5,1))
sigma <- 1.2
data.gen <- function(beta,sigma){
  e <- rnorm(n,0,1)
  X <- rnorm(n,0,sigma)
  constant <- rep(1,n)
  X <- cbind(constant,X,X^2,X^3,X^4)
  Y <- X%*%beta+e
  data <- data.frame(Y,X)
}
sim <- data.gen(beta,1.2)

f_hat <- lm(Y~.-1,data=sim)
library(sandwich)
C_hat <- vcovHC(f_hat)
x0 <- seq(-1,1,0.1)
CI <- matrix(NA,length(x0),2)
Y_hat <- c()
coe <- as.matrix(coef(f_hat),c(5,1))
for (i in 1:length(x0)) {
  l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
  Y_hat[i] <- t(coe)%*%l
  var_f <- t(l)%*%C_hat%*%l
  se <- sqrt(var_f)
  CI[i,] <- c(Y_hat[i]-2*se,Y_hat[i]+2*se)
}
CI

######### b
boot.fn <- function(data,index){
  coef(lm(Y~.,data=sim,subset = index))
}
boot.result <- boot(sim,boot.fn,R=100)
boot.result <- boot.result$t
boot.result <- boot.result[,-2]
boot.result <- as.matrix(boot.result)
CI_data.boot <- matrix(NA,length(x0),100)
# rows for each x0; columns for 100 times resample
for (i in 1:length(x0)) {
  l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
  CI_data.boot[i,] <- t(boot.result%*%l)
}
## sort CI data 
for (i in 1:length(x0)) {
  CI_data.boot[i,] <- sort(CI_data.boot[i,])
}
## find corresponding critical value in bootstrap data 
# according to theoretical CI
###################################################
# this part is fancy shit, useless
precent.critical <- matrix(NA,length(x0),2)
for (i in 1:length(x0)) {
  select_l <- max(CI_data.boot[i,CI[i,1]>=CI_data.boot[i,]])
  precent.critical[i,1] <- select_l
  select_h <- min(CI_data.boot[i,CI[i,2]<=CI_data.boot[i,]])
  precent.critical[i,2] <- select_h
  if(is.finite(precent.critical[i,1])==F){
    precent.critical[i,1] <- min(CI_data.boot[i,])}
  if(is.finite(precent.critical[i,2])==F){
    precent.critical[i,2] <- max(CI_data.boot[i,])}
}
##################################################

## now find naive CI in bootstrap data
CI_bootstrap <- matrix(NA,length(x0),2)
for (i in 1:length(x0)) {
  q <- unname(quantile(CI_data.boot[i,],
                       c(0.025,0.975)))
  CI_bootstrap[i,1] <- q[1]
  CI_bootstrap[i,2] <- q[2]
}
CI_bootstrap

plot(x0,CI_bootstrap[,1],col="red","l")
lines(x0,CI_bootstrap[,2],col="red")
lines(x0,CI[,1],col="blue")
lines(x0,CI[,2],col="blue")

plot(rep(x0[1],100),CI_data.boot[1,],col="red",
     ylim = c(-4,4),xlim = c(-1,1),
     xlab="x0",ylab = "f(x0)",
     main = "Bootstrap result & 2SE CI",pch=19)
for (i in 2:length(x0)) {
  points(rep(x0[i],100),CI_data.boot[i,],col="red")
}
lines(x0,CI[,1],col="blue",lwd=2)
lines(x0,CI[,2],col="blue",lwd=2)
grid()
# c
# coverage probability
####################################################
################################################
################# here we adjust 
c.p <- c()
x0 <- seq(-1,1,0.1)
boot.fn.cp <- function(data,index){
  f_hat <- lm(Y~.,data=sim,subset = index)
  C_hat <- vcovHC(f_hat)
  x0 <- seq(-1,1,0.1)
  CI <- matrix(NA,length(x0),2)
  Y_hat <- c()
  coe <- as.matrix(coef(f_hat),c(5,1))
  for (i in 1:length(x0)) {
    l0 <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
                c(5,1))
    Y_hat[i] <- t(coe[-2])%*%l0
    var_f <- t(l0)%*%C_hat%*%l0
    se <- sqrt(var_f)
    CI[i,] <- c(Y_hat[i]-2*se,Y_hat[i]+2*se)
  }
  return(CI)
}
boot.result <- boot(sim,boot.fn.cp,R=100)
boot.result <- boot.result$t
## now boot.result is 100by42, 100 rows is 100 resampling
# and for 42 columes, each i and i+21 is a pair of CI
CI.boot.stor <- matrix(NA,length(x0),2)
judge <- matrix(NA,length(x0),100)
for (i in 1:100) {
  for (j in 1:length(x0)) {
    CI.boot.stor[j,1] <- boot.result[i,j]
    CI.boot.stor[j,2] <- boot.result[i,j+length(x0)]
    judge[j,i] <- (Y_hat[j]>= CI.boot.stor[j,1]&
      Y_hat[j]<= CI.boot.stor[j,2]   )              
    }
}
for (i in 1:length(x0)) {
  c.p[i] <- mean(judge[i,])
}
c.p

# d
interval_length <- precent.critical[,2]-precent.critical[,1]
interval_length

##################################################
################  ex 2 ###########################
##################################################
set.seed(527)
n <- 1000
beta <- matrix(c(1,1.5,-1.5,1.5,0.5),c(5,1))
sigma <- 1.2
x0 <- seq(-1,1,0.1)
S <- 100 # simulation times
## for simplicity, we create some functions

# compute theoretical CI with 2SE
theoretical.CI_fun <- function(sim){
  f_hat <- lm(Y~.-1,data=sim)
  C_hat <- vcovHC(f_hat)
  CI <- matrix(NA,length(x0),2)
  Y_hat <- c()
  coe <- as.matrix(coef(f_hat),c(5,1))
  for (i in 1:length(x0)) {
   l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
   Y_hat[i] <- t(coe)%*%l
   var_f <- t(l)%*%C_hat%*%l
   se <- sqrt(var_f)
   CI[i,] <- c(Y_hat[i]-2*se,Y_hat[i]+2*se)
 }
 return(CI)
}
# compute bootstrap data
bootstrap.data_fun <- function(sim){
  boot.result <- boot(sim,boot.fn,R=100)
  boot.result <- boot.result$t
  boot.result <- boot.result[,-2]
  boot.result <- as.matrix(boot.result)
  CI_data.boot <- matrix(NA,length(x0),100)
  for (i in 1:length(x0)) {
   l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
   CI_data.boot[i,] <- t(boot.result%*%l)
 }
  return(CI_data.boot)
}

## compute bootstrap CI at 95%
CI_bootstrap_fun <- function(CI_data.boot){
  CI_bootstrap <- matrix(NA,length(x0),2)
 for (i in 1:length(x0)) {
   q <- unname(quantile(CI_data.boot[i,],
                        c(0.025,0.975)))
   CI_bootstrap[i,1] <- q[1]
   CI_bootstrap[i,2] <- q[2]
 }
  return(CI_bootstrap)
}

### all we need finally are stored in here
final_result <- matrix(NA,length(x0),8)
colnames(final_result) <- c("theo.CI.low","theo.CI.up",
      "boot.CI.low","boot.CI.up","coverage.Prob.SE",
      "coverage.Prob.boot","length.theo.CI","length.boot.CI")
frc1 <- matrix(NA,length(x0),S)
frc2 <- matrix(NA,length(x0),S)
frc3 <- matrix(NA,length(x0),S)
frc4 <- matrix(NA,length(x0),S)
frc7 <- matrix(NA,length(x0),S)
frc8 <- matrix(NA,length(x0),S)

for (i in 1:S) {
 sim <- data.gen(beta,1.2)
 CI <- theoretical.CI_fun(sim)
 frc1[,i] <- CI[,1]
 frc2[,i] <- CI[,2]
 CI_data.boot <- bootstrap.data_fun(sim)
 CI_bootstrap <- CI_bootstrap_fun(CI_data.boot)
 frc3[,i] <- CI_bootstrap[,1]
 frc4[,i] <- CI_bootstrap[,2]
 frc7[,i] <- CI[,2]-CI[,1]
 frc8[,i] <- CI_bootstrap[,2]-CI_bootstrap[,1]
}
# stupid value given procedure
for (i in 1:length(x0)) {
  final_result[i,1] <- mean(frc1[i,])
}
for (i in 1:length(x0)) {
  final_result[i,2] <- mean(frc2[i,])
}
for (i in 1:length(x0)) {
  final_result[i,3] <- mean(frc3[i,])
}
for (i in 1:length(x0)) {
  final_result[i,4] <- mean(frc4[i,])
}
for (i in 1:length(x0)) {
  final_result[i,7] <- mean(frc7[i,])
}
for (i in 1:length(x0)) {
  final_result[i,8] <- mean(frc8[i,])
}
y_true <- c()
for (i in 1:length(x0)) {
  l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
  y_true[i] <- t(beta)%*%l
}
c.p <- c()
c.p.boot <- c()
judge <- matrix(NA,length(x0),100)
judge2 <- matrix(NA,length(x0),100)
for (j in 1:length(x0)) {
  judge[j,] <- (y_true[j]>= frc1[j,]&
                 y_true[j]<= frc2[j,])
  judge2[j,] <- (y_true[j]>= frc3[j,]&
                   y_true[j]<= frc4[j,])
}
for (i in 1:length(x0)) {
  c.p[i] <- mean(judge[i,])
  c.p.boot[i] <- mean(judge2[i,]) 
}
final_result[,5] <- c.p
final_result[,6] <- c.p.boot
final_result
# final result shows that 2SE over estimate 95%CI
# change 2 to 1.96 will improve, but doesn't matter
# now let's change dgp to improve bootstrap

##
### what is a good CI?? 
# if we say good CI is which closer to real 1-alpha
# then less n 
# non-normal error term
# and if we think the shorter CI length, the better
# both above should hold
# change to hetroscho....error 
n <- 1000
#n <- 100
beta <- matrix(c(1,1.5,-1.5,1.5,0.5),c(5,1))
sigma <- 1.2
x0 <- seq(-1,1,0.1)
S <- 100 # simulation times
set.seed(527)
data.gen.v2 <- function(beta,sigma){
  e <- runif(n,-sqrt(5),sqrt(5))
  #e <- rnorm(n,0,1)
  X <- rnorm(n,0,sigma)
  constant <- rep(1,n)
  X <- cbind(constant,X,X^2,X^3,X^4)
  Y <- X%*%beta+e
  data <- data.frame(Y,X)
}

final_result1 <- matrix(NA,length(x0),8)
colnames(final_result1) <- c("theo.CI.low","theo.CI.up",
                            "boot.CI.low","boot.CI.up","coverage.Prob.SE",
                            "coverage.Prob.boot","length.theo.CI","length.boot.CI")
frc1 <- matrix(NA,length(x0),S)
frc2 <- matrix(NA,length(x0),S)
frc3 <- matrix(NA,length(x0),S)
frc4 <- matrix(NA,length(x0),S)
frc7 <- matrix(NA,length(x0),S)
frc8 <- matrix(NA,length(x0),S)
## this loop may have warning message, ignore it
for (i in 1:S) {
  sim <- data.gen.v2(beta,1.2)
  CI <- theoretical.CI_fun(sim)
  frc1[,i] <- CI[,1]
  frc2[,i] <- CI[,2]
  CI_data.boot <- bootstrap.data_fun(sim)
  CI_bootstrap <- CI_bootstrap_fun(CI_data.boot)
  frc3[,i] <- CI_bootstrap[,1]
  frc4[,i] <- CI_bootstrap[,2]
  frc7[,i] <- CI[,2]-CI[,1]
  frc8[,i] <- CI_bootstrap[,2]-CI_bootstrap[,1]
}
for (i in 1:length(x0)) {
  final_result1[i,1] <- mean(frc1[i,])
}
for (i in 1:length(x0)) {
  final_result1[i,2] <- mean(frc2[i,])
}
for (i in 1:length(x0)) {
  final_result1[i,3] <- mean(frc3[i,])
}
for (i in 1:length(x0)) {
  final_result1[i,4] <- mean(frc4[i,])
}
for (i in 1:length(x0)) {
  final_result1[i,7] <- mean(frc7[i,])
}
for (i in 1:length(x0)) {
  final_result1[i,8] <- mean(frc8[i,])
}
y_true <- c()
for (i in 1:length(x0)) {
  l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
  y_true[i] <- t(beta)%*%l
}
c.p <- c()
c.p.boot <- c()
judge <- matrix(NA,length(x0),100)
judge2 <- matrix(NA,length(x0),100)
for (j in 1:length(x0)) {
  judge[j,] <- (y_true[j]>= frc1[j,]&
                  y_true[j]<= frc2[j,])
  judge2[j,] <- (y_true[j]>= frc3[j,]&
                   y_true[j]<= frc4[j,])
}
for (i in 1:length(x0)) {
  c.p[i] <- mean(judge[i,])
  c.p.boot[i] <- mean(judge2[i,]) 
}
final_result1[,5] <- c.p
final_result1[,6] <- c.p.boot

################################################
## small sample 
set.seed(527)
data.gen.v3 <- function(beta,sigma){
  #e <- runif(n,-sqrt(5),sqrt(5))
  e <- rnorm(n,0,1)
  X <- rnorm(n,0,sigma)
  constant <- rep(1,n)
  X <- cbind(constant,X,X^2,X^3,X^4)
  Y <- X%*%beta+e
  data <- data.frame(Y,X)
}
n <- 100
final_result2 <- matrix(NA,length(x0),8)
colnames(final_result2) <- c("theo.CI.low","theo.CI.up",
                            "boot.CI.low","boot.CI.up","coverage.Prob.SE",
                            "coverage.Prob.boot","length.theo.CI","length.boot.CI")
frc1 <- matrix(NA,length(x0),S)
frc2 <- matrix(NA,length(x0),S)
frc3 <- matrix(NA,length(x0),S)
frc4 <- matrix(NA,length(x0),S)
frc7 <- matrix(NA,length(x0),S)
frc8 <- matrix(NA,length(x0),S)

for (i in 1:S) {
  sim <- data.gen.v3(beta,1.2)
  CI <- theoretical.CI_fun(sim)
  frc1[,i] <- CI[,1]
  frc2[,i] <- CI[,2]
  CI_data.boot <- bootstrap.data_fun(sim)
  CI_bootstrap <- CI_bootstrap_fun(CI_data.boot)
  frc3[,i] <- CI_bootstrap[,1]
  frc4[,i] <- CI_bootstrap[,2]
  frc7[,i] <- CI[,2]-CI[,1]
  frc8[,i] <- CI_bootstrap[,2]-CI_bootstrap[,1]
}
for (i in 1:length(x0)) {
  final_result2[i,1] <- mean(frc1[i,])
}
for (i in 1:length(x0)) {
  final_result2[i,2] <- mean(frc2[i,])
}
for (i in 1:length(x0)) {
  final_result2[i,3] <- mean(frc3[i,])
}
for (i in 1:length(x0)) {
  final_result2[i,4] <- mean(frc4[i,])
}
for (i in 1:length(x0)) {
  final_result2[i,7] <- mean(frc7[i,])
}
for (i in 1:length(x0)) {
  final_result2[i,8] <- mean(frc8[i,])
}
y_true <- c()
for (i in 1:length(x0)) {
  l <- matrix(c(1,x0[i],x0[i]^2,x0[i]^3,x0[i]^4),
              c(5,1))
  y_true[i] <- t(beta)%*%l
}
c.p <- c()
c.p.boot <- c()
judge <- matrix(NA,length(x0),100)
judge2 <- matrix(NA,length(x0),100)
for (j in 1:length(x0)) {
  judge[j,] <- (y_true[j]>= frc1[j,]&
                  y_true[j]<= frc2[j,])
  judge2[j,] <- (y_true[j]>= frc3[j,]&
                   y_true[j]<= frc4[j,])
}
for (i in 1:length(x0)) {
  c.p[i] <- mean(judge[i,])
  c.p.boot[i] <- mean(judge2[i,]) 
}
final_result2[,5] <- c.p
final_result2[,6] <- c.p.boot

#View(final_result1)

# the error term changed case
plot(x0,final_result[,5],"l",lwd=2,
     ylim=c(0.9,1),main = "Coverage Probability",
     ylab = "Coverage Probability")
lines(x0,final_result1[,5],lty=2,lwd=2)
lines(x0,final_result[,6],lty=1,col='red',lwd=2)
lines(x0,final_result1[,6],lty=2,col='red',lwd=2)
legend("bottomright",lty = c(1,2,1,2),lwd = c(2,2,2,2),
       col=c("black","black",'red','red'),
       legend = c("Normal.2SE","Uniform.2SE",
                  "Normal.boot","Uniform.boot"))
lines(x0,rep(0.95,length(x0)),col='blue',lwd=2)
# the small sample case 
plot(x0,final_result[,5],"l",lwd=2,
     ylim=c(0.9,1),main = "Coverage Probability",
     ylab = "Coverage Probability")
lines(x0,final_result2[,5],lty=2,lwd=2)
lines(x0,final_result[,6],lty=1,col='red',lwd=2)
lines(x0,final_result2[,6],lty=2,col='red',lwd=2)
legend("bottomright",lty = c(1,2,1,2),lwd = c(2,2,2,2),
       col=c("black","black",'red','red'),
       legend = c("1000.2SE","100.2SE",
                  "1000.boot","100.boot"))
lines(x0,rep(0.95,length(x0)),col='blue',lwd=2)

# plot the difference of coverage prob.
plot(x0,rep(0.95,length(x0))-final_result[,5],"l",lwd=2,
 ylim=c(-0.05,0.05),ylab = "95% - Coverage.Prob.",
 main = "Deviation of Coverage Prob. from 95%")
lines(x0,rep(0.95,length(x0))-final_result2[,5],lty=2,lwd=2)
lines(x0,rep(0.95,length(x0))-final_result[,6],lty=1,col='red',lwd=2)
lines(x0,rep(0.95,length(x0))-final_result2[,6],lty=2,col='red',lwd=2)
legend("topright",lty = c(1,2,1,2),lwd = c(2,2,2,2),
       col=c("black","black",'red','red'),
       legend = c("1000.2SE","100.2SE",
                  "1000.boot","100.boot"))
lines(x0,rep(0,length(x0)),col='blue',lwd=2)
# plot the deviation of interval length
# within groups
plot(x0,final_result[,7]-final_result1[,7],"l",
     ylim=c(-0.056,-0.045),
     main = "Advantage of Bootstrap CI over 2SE CI",
     ylab = "(Normal case) - (Uniform case)",lwd=2)
lines(x0,final_result[,8]-final_result1[,8],lty=2,lwd=2)
legend("topright",lty = 1:2,
       lwd = c(2,2),legend = c("Theoretical",
                                       "Bootstrap"))

plot(x0,final_result[,7]-final_result2[,7],"l",
     ylim=c(-0.52,-0.35),
     main = "Advantage of Bootstrap CI over 2SE CI",
     ylab = "(n=1000) - (n=100)",lwd=2)
lines(x0,final_result[,8]-final_result2[,8],lty=2,lwd=2)
legend("topright",lwd = c(2,2),lty = 1:2,legend = c("Theoretical",
                                       "Bootstrap"))


# between groups
plot(x0,final_result1[,7]-final_result1[,8],
   ylim = c(0.008,0.016),lwd=2,
     ylab = "2SE CI - Bootstrap CI",
     main = "Advantage of Bootstrap CI over 2SE CI","l")
# ylim for error setting
lines(x0,final_result[,7]-final_result[,8],lty=2,lwd=2)
legend("bottomleft",lty = 1:2,lwd = c(2,2),
       legend = c("Uniform Error","Normal Error"))


plot(x0,final_result2[,7]-final_result2[,8],lwd=2,
     ylim = c(0.01,0.07),ylab = "2SE CI - Bootstrap CI",
     main = "Advantage of Bootstrap CI over 2SE CI","l")
# ylim for small n setting
lines(x0,final_result[,7]-final_result[,8],lty=2,lwd=2)
legend("topleft",lwd = c(2,2),
       lty = 1:2,legend = c("n = 100","n = 1000"))
