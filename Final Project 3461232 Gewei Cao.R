library(MASS)
n <- 100
beta <- c(1,1,1)
sigmaz <- matrix(NA,7,7)
mu <- rep(0,7)
for (i in 1:7) {
  for (j in 1:7) {
    sigmaz[i,j] <- 0.6^abs(i-j)
  }
}
dgp_low <- function(n,mu,sigmaz){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  z <- mvrnorm(n,mu,sigmaz)
  x1 <- z[,1]+z[,2]+z[,3]+z[,4]+z[,5]+z[,6]+z[,7]+error_c
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+error_y
  data <- data.frame(y,x1,z)
  return(data)
}
training <- dgp_low(n,mu,sigmaz)
colnames(training) <- c("y","x","z1","z2","z3",
                        "z4","z5","z6","z7")
# column 1 is y, 2 is x, 3 to 9 are z

library(glmnet)
library(pls)
library(AER)
T <- 200
Coef <- matrix(NA,T,5)
colnames(Coef) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")

set.seed(527)
for (i in 1:T) {
training <- dgp_low(n,mu,sigmaz)
colnames(training) <- c("y","x","z1","z2","z3","z4","z5","z6","z7")
Instruments_tr <- as.matrix(training[,3:9])
## benchmark OLS and IV
NaiveOLS <- lm(y~x,data=training)
Coef[i,1] <- coef(NaiveOLS)[2]

IVbenchmark <- ivreg(training$y~training$x|Instruments_tr)
Coef[i,2] <- coef(IVbenchmark)[2]

## lasso
grid <- 10^seq(-3,3,length=200)
Lasso <- glmnet(Instruments_tr,training$x,
                lambda = grid,alpha=1)
cv.out<-cv.glmnet(Instruments_tr, training$x,alpha=1)
best_lambda<-cv.out$lambda.min
lasso.xhat<-predict(Lasso, s=best_lambda,
                    newx = Instruments_tr)
Lasso.2sls <- lm(training$y~lasso.xhat)
Coef[i,3] <- coef(Lasso.2sls)[2]

## post-lasso
cv.out.post<-cv.glmnet(Instruments_tr, training$x,
                       alpha=1)
best_lambda<-cv.out.post$lambda.min
Lasso.post <- glmnet(Instruments_tr,training$x,
                lambda = best_lambda,alpha=1)
coelasso <- coef(Lasso.post)
select <- ifelse(coelasso==0,FALSE,TRUE)[-1]
Instruments_tr2 <- as.matrix(training[,3:9][,select])
IVpostLasso <- ivreg(training$y~training$x|Instruments_tr2)
Coef[i,4] <- coef(IVpostLasso)[2]

# PCA
PCR<-pcr(training$x~Instruments_tr,
             scale=TRUE,validation="CV")
cverr <- RMSEP(PCR)$val[1,,-1]
imin<-which.min(cverr)
PCA_IV <- PCR$scores[,1:imin]
IVPCR <- ivreg(training$y~training$x|PCA_IV)
Coef[i,5] <- coef(IVPCR)[2]

}

# write summary table
Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef[,i])
  Coef_summ[2,i] <- var(Coef[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 1
plot(density(Coef[,1]),lwd=2,xlim=c(0.8,1.25),ylim=c(0,9.5),
     col="paleturquoise4",main = "Distribution of beta (Figure 1)")
lines(density(Coef[,2]),col="red",lwd=3)
lines(density(Coef[,3]),col="blue",lwd=2)
lines(density(Coef[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topright",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,11),seq(0,10,by=1),lwd=2)
# Figure 2
boxplot(Coef,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        main="Beta  (Figure 2)",lwd=2)
abline(1,0,col="darkred",lwd=2)



################################################
################################################
################# Case 2########################
################################################
################################################
#rm(list=ls())
n <- 100
beta <- c(1,1,1)
sigmaz <- matrix(NA,7,7)
mu <- rep(0,7)
for (i in 1:7) {
  for (j in 1:7) {
    sigmaz[i,j] <- 0.6^abs(i-j)
  }
}
beta_tild <- c(1,0.8,0.6,0.5,0.3,0.2,0.02)
dgp_low2 <- function(n,mu,sigmaz){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  z <- mvrnorm(n,mu,sigmaz)
  x1 <- z%*%beta_tild+error_c
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+error_y
  data <- data.frame(y,x1,z)
  return(data)
}
training <- dgp_low2(n,mu,sigmaz)
colnames(training) <- c("y","x","z1","z2","z3",
                        "z4","z5","z6","z7")
# Case 2--weak instruments 
T <- 200
Coef2 <- matrix(NA,T,5)
colnames(Coef2) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
for (i in 1:T) {
  training <- dgp_low2(n,mu,sigmaz)
  colnames(training) <- c("y","x","z1","z2","z3","z4","z5","z6","z7")
  Instruments_tr <- as.matrix(training[,3:9])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x,data=training)
  Coef2[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x|Instruments_tr)
  Coef2[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  Lasso <- glmnet(Instruments_tr,training$x,
                  lambda = grid,alpha=1)
  cv.out<-cv.glmnet(Instruments_tr, training$x,alpha=1)
  best_lambda<-cv.out$lambda.min
  lasso.xhat<-predict(Lasso, s=best_lambda,
                      newx = Instruments_tr)
  Lasso.2sls <- lm(training$y~lasso.xhat)
  Coef2[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out.post<-cv.glmnet(Instruments_tr, training$x,
                         alpha=1)
  best_lambda<-cv.out.post$lambda.min
  Lasso.post <- glmnet(Instruments_tr,training$x,
                       lambda = best_lambda,alpha=1)
  coelasso <- coef(Lasso.post)
  select <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  Instruments_tr2 <- as.matrix(training[,3:9][,select])
  IVpostLasso <- ivreg(training$y~training$x|Instruments_tr2)
  Coef2[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x~Instruments_tr,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  IVPCR <- ivreg(training$y~training$x|PCA_IV)
  Coef2[i,5] <- coef(IVPCR)[2]
}

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef2[,i])
  Coef_summ[2,i] <- var(Coef2[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 3
plot(density(Coef2[,1]),lwd=2,xlim=c(0.6,1.6),
     col="paleturquoise4",main = "Distribution of beta (Figure 3)")
lines(density(Coef2[,2]),col="red",lwd=3)
lines(density(Coef2[,3]),col="blue",lwd=2)
lines(density(Coef2[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef2[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topright",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,11),seq(0,10,by=1),lwd=2)
# Figure 4
boxplot(Coef2,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 4)")
abline(1,0,col="darkred",lwd=2)

################################################
################################################
################# Case 3########################
################################################
################################################
#rm(list=ls())
beta <- c(1,1,1)
sigmaz1 <- matrix(NA,7,7)
mu <- rep(0,7)
for (i in 1:7) {
  for (j in 1:7) {
    sigmaz1[i,j] <- 0.9^abs(i-j)
  }
}

training <- dgp_low(n,mu,sigmaz1)
colnames(training) <- c("y","x","z1","z2","z3",
                        "z4","z5","z6","z7")
# Case 3--weak instruments 
T <- 200
Coef3 <- matrix(NA,T,5)
colnames(Coef3) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
for (i in 1:T) {
  training <- dgp_low(n,mu,sigmaz1)
  colnames(training) <- c("y","x","z1","z2","z3","z4","z5","z6","z7")
  Instruments_tr <- as.matrix(training[,3:9])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x,data=training)
  Coef3[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x|Instruments_tr)
  Coef3[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  Lasso <- glmnet(Instruments_tr,training$x,
                  lambda = grid,alpha=1)
  cv.out<-cv.glmnet(Instruments_tr, training$x,alpha=1)
  best_lambda<-cv.out$lambda.min
  lasso.xhat<-predict(Lasso, s=best_lambda,
                      newx = Instruments_tr)
  Lasso.2sls <- lm(training$y~lasso.xhat)
  Coef3[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out.post<-cv.glmnet(Instruments_tr, training$x,
                         alpha=1)
  best_lambda<-cv.out.post$lambda.min
  Lasso.post <- glmnet(Instruments_tr,training$x,
                       lambda = best_lambda,alpha=1)
  coelasso <- coef(Lasso.post)
  select <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  Instruments_tr2 <- as.matrix(training[,3:9][,select])
  IVpostLasso <- ivreg(training$y~training$x|Instruments_tr2)
  Coef3[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x~Instruments_tr,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  IVPCR <- ivreg(training$y~training$x|PCA_IV)
  Coef3[i,5] <- coef(IVPCR)[2]
}

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef3[,i])
  Coef_summ[2,i] <- var(Coef3[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 3
plot(density(Coef3[,1]),lwd=2,xlim=c(0.8,1.2),
     col="paleturquoise4",main = "Distribution of beta (Figure 5)")
lines(density(Coef3[,2]),col="red",lwd=3)
lines(density(Coef3[,3]),col="blue",lwd=2)
lines(density(Coef3[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef3[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topright",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,14),seq(0,13,by=1),lwd=2)
# Figure 4
boxplot(Coef3,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 6)")
abline(1,0,col="darkred",lwd=2)


################################################
################################################
################# Case 4########################
################################################
################################################
#rm(list=ls())
set.seed(527)
beta <- c(1,1,1)
sigmaz3 <- matrix(NA,60,60)
mu <- rep(0,60)
for (i in 1:60) {
  for (j in 1:60) {
    sigmaz3[i,j] <- 0.6^abs(i-j)
  }
}
alpha <- c()
for (i in 1:60) {
  alpha[i] <- 0.7^i
} #two options this one is sparse one
#alpha <- sample(seq(0,1,by=0.001),60) #non-sparse one
dgp_low4 <- function(n,mu,sigmaz){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  x1 <- c()
  z <- mvrnorm(n,mu,sigmaz)
  for (i in 1:n) {
      x1[i] <- z[i,1:60]%*%alpha+error_c[i]
  }
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+error_y
  data <- data.frame(y,x1,z)
  return(data)
}
training <- dgp_low4(n,mu,sigmaz3)

barplot(alpha,ylim = c(0,1),names.arg = seq(1,60,by=1),
        xlab = "z",ylab = "alpha",
        main = "alpha for each z  (Figure 7)")
barplot(alpha,ylim = c(0,0.75),names.arg = seq(1,60,by=1),
        xlab = "z",ylab = "alpha",
        main = "alpha for each z  (Figure 19)")
# Case 4--High Dimensional 
T <- 200
Coef4 <- matrix(NA,T,5)
colnames(Coef4) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
for (i in 1:T) {
  training <- dgp_low4(n,mu,sigmaz3)
  Instruments_tr <- as.matrix(training[,3:62])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x1,data=training)
  Coef4[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x1|Instruments_tr)
  Coef4[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  Lasso <- glmnet(Instruments_tr,training$x1,
                  lambda = grid,alpha=1)
  cv.out<-cv.glmnet(Instruments_tr, training$x1,alpha=1)
  best_lambda<-cv.out$lambda.min
  lasso.xhat<-predict(Lasso, s=best_lambda,
                      newx = Instruments_tr)
  Lasso.2sls <- lm(training$y~lasso.xhat)
  Coef4[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out.post<-cv.glmnet(Instruments_tr, training$x1,
                         alpha=1)
  best_lambda<-cv.out.post$lambda.min
  Lasso.post <- glmnet(Instruments_tr,training$x1,
                       lambda = best_lambda,alpha=1)
  coelasso <- coef(Lasso.post)
  select <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  Instruments_tr2 <- as.matrix(training[,3:62][,select])
  IVpostLasso <- ivreg(training$y~training$x1|Instruments_tr2)
  Coef4[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x1~Instruments_tr,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  IVPCR <- ivreg(training$y~training$x1|PCA_IV)
  Coef4[i,5] <- coef(IVPCR)[2]
}

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef4[,i])
  Coef_summ[2,i] <- var(Coef4[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 8
plot(density(Coef4[,1]),lwd=2,xlim=c(0.8,1.2),
     col="paleturquoise4",main = "Distribution of beta (Figure 8)")
lines(density(Coef4[,2]),col="red",lwd=3)
lines(density(Coef4[,3]),col="blue",lwd=2)
lines(density(Coef4[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef4[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topright",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)
# Figure 9
boxplot(Coef4,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 9)")
abline(1,0,col="darkred",lwd=2)

# figure 20 with sparse alpha
plot(density(Coef4[,1]),lwd=2,xlim=c(0.4,2),
     col="paleturquoise4",main = "Distribution of beta (Figure 20)")
lines(density(Coef4[,2]),col="red",lwd=3)
lines(density(Coef4[,3]),col="blue",lwd=2)
lines(density(Coef4[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef4[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)


################################################
################################################
################# Case 5########################
################################################
################################################
#rm(list=ls())
set.seed(527)
beta2 <- c(0.7,0.5,0.35,0.4,0.2,0.1)
mu <- rep(0,7)
alpha <- sample(seq(0,1,by=0.001),7)
muc <- rep(0,3)
sigmac <- matrix(c(1,0.7,0.49,0.7,1,0.7,0.49,0.7,1),3,3,byrow = TRUE)
dgp_low5 <- function(n,mu,sigmaz){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  x3 <- rnorm(n,0,2)+error_c
  x4 <- runif(n,-1,1)+error_c
  x5 <- rnorm(n,0,1.5)+error_c
  X678 <- mvrnorm(n,muc,sigmac)
  x1 <- c()
  X <- cbind(x3,x4,x5,X678)
  z <- mvrnorm(n,mu,sigmaz)
  for (i in 1:n) {
    x1[i] <- z[i,1:7]%*%alpha+X[i,]%*%beta2+error_c[i]
  }
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+X%*%beta2+error_y
  data <- data.frame(y,x1,z,X)
  return(data)
}
training <- dgp_low5(n,mu,sigmaz)


# Case 5--Control Variables
set.seed(527)
T <- 200
Coef5 <- matrix(NA,T,5)
colnames(Coef5) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
for (i in 1:T) {
  training <- dgp_low5(n,mu,sigmaz)
  Instruments_tr <- as.matrix(training[,3:9])
  Controls <- as.matrix(training[,10:15])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x1+x3+x4+x5+V4+V5+V6,data=training)
  Coef5[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x1+Controls|Instruments_tr+Controls)
  Coef5[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  Lasso <- glmnet(cbind(Instruments_tr,Controls),training$x1,
                  lambda = grid,alpha=1)
  cv.out<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,alpha=1)
  best_lambda<-cv.out$lambda.min
  lasso.xhat<-predict(Lasso, s=best_lambda,
                      newx = cbind(Instruments_tr,Controls))
  Lasso.2sls <- lm(training$y~lasso.xhat)
  Coef5[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out.post<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,
                         alpha=1)
  best_lambda<-cv.out.post$lambda.min
  Lasso.post <- glmnet(cbind(Instruments_tr,Controls),training$x1,
                       lambda = best_lambda,alpha=1)
  coelasso <- coef(Lasso.post)
  select <- ifelse(coelasso==0,FALSE,TRUE)[2:8]
  Instruments_tr2 <- as.matrix(training[,3:9][,select])
  IVpostLasso <- ivreg(training$y~training$x1+Controls|Instruments_tr2+Controls)
  Coef5[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x1~Instruments_tr+Controls,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  ttt <- lm(training$x1~PCA_IV)
  xp <- predict(ttt)
  IVPCR <- lm(training$y~xp+Controls)
  Coef5[i,5] <- coef(IVPCR)[2]
}

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef5[,i])
  Coef_summ[2,i] <- var(Coef5[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 10
plot(density(Coef5[,1]),lwd=2,xlim=c(0.7,2.1),
     col="paleturquoise4",main = "Distribution of beta (Figure 10)")
lines(density(Coef5[,2]),col="red",lwd=3)
lines(density(Coef5[,3]),col="blue",lwd=2)
lines(density(Coef5[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef5[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("top",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)

# figure 12
plot(density(Coef5[,1]),lwd=2,xlim=c(0.7,1.3),
     col="paleturquoise4",main = "Distribution of beta (Figure 12)")
lines(density(Coef5[,2]),col="red",lwd=3)
lines(density(Coef5[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef5[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","darkgreen","darkgoldenrod2"),
       lty = c(1,1,2,4),lwd=c(2,3,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)
# Figure 11
boxplot(Coef5,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 11)")
abline(1,0,col="darkred",lwd=2)




################################################
################################################
################# Case 6########################
################################################
################################################
#rm(list=ls())
###############################################
set.seed(527)
beta2 <- c(0.7,0.5,0.35,0.4,0.2,0.1)
mu <- rep(0,7)
alpha <- sample(seq(0,1,by=0.001),7)
muc <- rep(0,3)
sigmac <- matrix(c(1,0.7,0.49,0.7,1,0.7,0.49,0.7,1),3,3,byrow = TRUE)
capitalpie <- matrix(NA,6,7)
for (i in 1:6) {
  capitalpie[i,] <- sample(seq(0,1,by=0.001),7)
}
# exogenous controls
dgp_lowCh <- function(n,mu,sigmaz,capitalpie){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  x3 <- rnorm(n,0,2)
  x4 <- runif(n,-1,1)
  x5 <- rnorm(n,0,1.5)
  X678 <- mvrnorm(n,muc,sigmac)
  x1 <- c()
  X <- cbind(x3,x4,x5,X678)
  xi <- mvrnorm(n,mu,sigmaz)
  z <- X%*%capitalpie+xi
  for (i in 1:n) {
    x1[i] <- z[i,1:7]%*%alpha+X[i,]%*%beta2+error_c[i]
  }
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+X%*%beta2+error_y
  data <- data.frame(y,x1,z,X)
  return(data)
}

# endogenous controls
dgp_lowChendo <- function(n,mu,sigmaz,capitalpie){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  x3 <- rnorm(n,0,2)+error_c
  x4 <- runif(n,-1,1)+error_c
  x5 <- rnorm(n,0,1.5)+error_c
  X678 <- mvrnorm(n,muc,sigmac)
  x1 <- c()
  X <- cbind(x3,x4,x5,X678)
  xi <- mvrnorm(n,mu,sigmaz)
  z <- X%*%capitalpie+xi
  for (i in 1:n) {
    x1[i] <- z[i,1:7]%*%alpha+X[i,]%*%beta2+error_c[i]
  }
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+X%*%beta2+error_y
  data <- data.frame(y,x1,z,X)
  return(data)
}
training <- dgp_lowCh(n,mu,sigmaz,capitalpie)


# Case 6
set.seed(527)
T <- 200
Coef6 <- matrix(NA,T,5)
colnames(Coef6) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
for (i in 1:T) {
  training <- dgp_low5(n,mu,sigmaz)
  #training <- dgp_lowCh(n,mu,sigmaz,capitalpie)
  # you can choose to use setting in 4.6 or 4.5
  Instruments_tr <- as.matrix(training[,3:9])
  Controls <- as.matrix(training[,10:15])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x1+x3+x4+x5+V4+V5+V6,data=training)
  Coef6[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x1+Controls|Instruments_tr+Controls)
  Coef6[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  cv.out<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso <- glmnet(cbind(Instruments_tr,Controls),training$x1,
                  lambda = best_lambda,alpha=1,intercept = F)
  gammadelta <- coef(Lasso) #step 1 (above)
  cv.out<-cv.glmnet(Controls, training$y,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  theta <- coef(glmnet(Controls,training$y,
                        lambda = best_lambda,alpha=1,intercept = F))
  # step 2 (above)
  dhat <- cbind(Instruments_tr,Controls)%*%gammadelta[-1]
  cv.out<-cv.glmnet( Controls,dhat,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  vhat <- coef(glmnet( Controls,dhat,
            lambda = best_lambda,alpha=1,intercept = F))
  # step 3 (above)
  rouy <- training$y-Controls%*%theta[-1]
  roud <- training$x1-Controls%*%vhat[-1]
  vvv <- dhat-Controls%*%vhat[-1]
  Lasso.2sls <- ivreg(rouy~roud|vvv)
  # step 4 (above)
  Coef6[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet(cbind(Instruments_tr,Controls),training$x1,
      lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_1 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  i_t <- as.matrix(training[,3:15][,select_1])
  gammadelta <- coef(lm(training$x1~i_t-1))
 # step 1 (above)
  cv.out<-cv.glmnet(Controls, training$y,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet(Controls,training$y,
                       lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_2 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  CC <- as.matrix(training[,10:15][,select_2])
  theta <- coef(lm(training$y~CC-1))
 #step 2 (above)
  dhat <-cbind(Instruments_tr,Controls)[,select_1]%*%gammadelta
  cv.out<-cv.glmnet( Controls,dhat,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet( Controls,dhat,
   lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_3 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  CCd <- as.matrix(training[,10:15][,select_3])
  vhat <- coef(lm(dhat~CCd-1)) #step 3 (above)
  if(length(theta)==1){
    rouy <- training$y-Controls[,select_2]*theta
  }else{
     rouy <- training$y-Controls[,select_2]%*%theta
  }
  if(length(vhat)==1){
    roud <- training$x1-Controls[,select_3]*vhat
    vvv <- dhat-Controls[,select_3]*vhat
  }else{
    roud <- training$x1-Controls[,select_3]%*%vhat
    vvv <- dhat-Controls[,select_3]%*%vhat
    }
  IVpostLasso <- ivreg(rouy~roud|vvv) #step 4 (above)
  Coef6[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x1~Instruments_tr+Controls,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  ttt <- lm(training$x1~PCA_IV)
  xp <- predict(ttt)
  IVPCR <- lm(training$y~xp+Controls)
  Coef6[i,5] <- coef(IVPCR)[2]
 }

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef6[,i])
  Coef_summ[2,i] <- var(Coef6[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 14
plot(density(Coef6[,1]),lwd=2,xlim=c(0.5,1.8),
     col="paleturquoise4",
     main = "Distribution of beta (Figure 14)")
lines(density(Coef6[,2]),col="red",lwd=3)
lines(density(Coef6[,3]),col="blue",lwd=2)
lines(density(Coef6[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef6[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)
# Figure 16
plot(density(Coef6[,1]),lwd=2,xlim=c(0.7,1.4),
     col="paleturquoise4",
     main = "Distribution of beta (Figure 16)")
lines(density(Coef6[,2]),col="red",lwd=3)
lines(density(Coef6[,3]),col="blue",lwd=2)
lines(density(Coef6[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef6[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)

# Figure 15
boxplot(Coef6,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 15)")
abline(1,0,col="darkred",lwd=2)



################################################
################################################
################# Case 4.6.4####################
################################################
################################################
#rm(list=ls())
###############################################
set.seed(527)
beta2 <- c(0.7,0.5,0.35,0.4,0.2,0.1)
mu <- rep(0,60)
alpha <- c()
for (i in 1:60) {
  alpha[i] <- 0.7^i
} # for simulation in 4.6.4, sparse case
alpha <- sample(seq(0,1,by=0.001),60) # for figure 21, non-sparse case
muc <- rep(0,3)
sigmac <- matrix(c(1,0.7,0.49,0.7,1,0.7,0.49,0.7,1),3,3,byrow = TRUE)
capitalpie60 <- matrix(NA,6,60)
for (i in 1:6) {
  capitalpie60[i,] <- sample(seq(0,1,by=0.001),60)
}
sigmaxi <- matrix(NA,60,60)
for (i in 1:60) {
  for (j in 1:60) {
    sigmaxi[i,j] <- 0.6^abs(i-j)
  }
}
dgp_lowCh60 <- function(n,mu,sigmaxi,capitalpie){
  error_c <- rnorm(n,0,1)
  delta <- runif(n,-1,1)
  x3 <- rnorm(n,0,2)
  x4 <- runif(n,-1,1)
  x5 <- rnorm(n,0,1.5)
  X678 <- mvrnorm(n,muc,sigmac)
  x1 <- c()
  X <- cbind(x3,x4,x5,X678)
  xi <- mvrnorm(n,mu,sigmaxi)
  z <- X%*%capitalpie+xi
  for (i in 1:n) {
    x1[i] <- z[i,]%*%alpha+X[i,]%*%beta2+error_c[i]
  }
  x2 <- 1+error_c
  error_y <- x2+delta
  y <- 1+x1+x2+X%*%beta2+error_y
  data <- data.frame(y,x1,z,X)
  return(data)
}
training <- dgp_lowCh60(n,mu,sigmaxi,capitalpie60)


# Case 6.n
set.seed(527)
T <- 200
Coef60 <- matrix(NA,T,5)
Coef601 <- matrix(NA,T,5) # only for complete figure 21
colnames(Coef60) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
colnames(Coef601) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
set.seed(527)
# if want to simulate 4.6.4, remember to change Coef601 into Coef60
for (i in 1:T) {
  training <- dgp_lowCh60(n,mu,sigmaxi,capitalpie60)
  Instruments_tr <- as.matrix(training[,3:62])
  Controls <- as.matrix(training[,63:68])
  ## benchmark OLS and IV
  NaiveOLS <- lm(y~x1+x3+x4+x5+V4+V5+V6,data=training)
  Coef601[i,1] <- coef(NaiveOLS)[2]
  IVbenchmark <- ivreg(training$y~training$x1+Controls|Instruments_tr+Controls)
  Coef601[i,2] <- coef(IVbenchmark)[2]
  ## lasso
  grid <- 10^seq(-3,3,length=200)
  cv.out<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso <- glmnet(cbind(Instruments_tr,Controls),training$x1,
                  lambda = best_lambda,alpha=1,intercept = F)
  gammadelta <- coef(Lasso) #step 1 (above)
  cv.out<-cv.glmnet(Controls, training$y,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  theta <- coef(glmnet(Controls,training$y,
                       lambda = best_lambda,alpha=1,intercept = F))
  # step 2 (above)
  dhat <- cbind(Instruments_tr,Controls)%*%gammadelta[-1]
  cv.out<-cv.glmnet( Controls,dhat,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  vhat <- coef(glmnet( Controls,dhat,
                       lambda = best_lambda,alpha=1,intercept = F))
  # step 3 (above)
  rouy <- training$y-Controls%*%theta[-1]
  roud <- training$x1-Controls%*%vhat[-1]
  vvv <- dhat-Controls%*%vhat[-1]
  Lasso.2sls <- ivreg(rouy~roud|vvv)
  # step 4 (above)
  Coef601[i,3] <- coef(Lasso.2sls)[2]
  ## post-lasso
  cv.out<-cv.glmnet(cbind(Instruments_tr,Controls), training$x1,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet(cbind(Instruments_tr,Controls),training$x1,
                       lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_1 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  i_t <- as.matrix(training[,3:68][,select_1])
  gammadelta <- coef(lm(training$x1~i_t-1))
  # step 1 (above)
  cv.out<-cv.glmnet(Controls, training$y,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet(Controls,training$y,
                       lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_2 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  CC <- as.matrix(training[,63:68][,select_2])
  theta <- coef(lm(training$y~CC-1))
  #step 2 (above)
  dhat <-cbind(Instruments_tr,Controls)[,select_1]%*%gammadelta
  cv.out<-cv.glmnet( Controls,dhat,alpha=1,intercept = F)
  best_lambda<-cv.out$lambda.min
  Lasso.post <- glmnet( Controls,dhat,
                        lambda = best_lambda,alpha=1,intercept = F)
  coelasso <- coef(Lasso.post)
  select_3 <- ifelse(coelasso==0,FALSE,TRUE)[-1]
  CCd <- as.matrix(training[,63:68][,select_3])
  vhat <- coef(lm(dhat~CCd-1)) #step 3 (above)
  if(length(theta)==1){
    rouy <- training$y-Controls[,select_2]*theta
  }else{
    rouy <- training$y-Controls[,select_2]%*%theta
  }
  if(length(vhat)==1){
    roud <- training$x1-Controls[,select_3]*vhat
    vvv <- dhat-Controls[,select_3]*vhat
  }else{
    roud <- training$x1-Controls[,select_3]%*%vhat
    vvv <- dhat-Controls[,select_3]%*%vhat
  }
  IVpostLasso <- ivreg(rouy~roud|vvv) #step 4 (above)
  Coef601[i,4] <- coef(IVpostLasso)[2]
  # PCA
  PCR<-pcr(training$x1~Instruments_tr+Controls,
           scale=TRUE,validation="CV")
  cverr <- RMSEP(PCR)$val[1,,-1]
  imin<-which.min(cverr)
  PCA_IV <- PCR$scores[,1:imin]
  ttt <- lm(training$x1~PCA_IV)
  xp <- predict(ttt)
  IVPCR <- lm(training$y~xp+Controls)
  Coef601[i,5] <- coef(IVPCR)[2]
}

Coef_summ <- matrix(NA,2,5)
for (i in 1:5) {
  Coef_summ[1,i] <- mean(Coef60[,i])
  Coef_summ[2,i] <- var(Coef60[,i])
}
colnames(Coef_summ) <- c("Naive OLS","2SLS","Lasso","Post-Lasso","PCA")
write.csv(Coef_summ, "coefficient.csv")

# Figure 17
plot(density(Coef60[,1]),lwd=2,xlim=c(0.2,2),
     col="paleturquoise4",
     main = "Distribution of beta (Figure 17)")
lines(density(Coef60[,2]),col="red",lwd=3)
lines(density(Coef60[,3]),col="blue",lwd=2)
lines(density(Coef60[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef60[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)

# Figure 18
boxplot(Coef60,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 18)")
abline(1,0,col="darkred",lwd=2)

# for checking Figure 21
plot(density(Coef601[,1]),lwd=2,xlim=c(0.8,1.2),
     col="paleturquoise4",
     main = "Distribution of beta (Figure 21)")
lines(density(Coef601[,2]),col="red",lwd=3)
lines(density(Coef601[,3]),col="blue",lwd=2)
lines(density(Coef601[,4]),col="darkgreen",lwd=2,lty=2)
lines(density(Coef601[,5]),col="darkgoldenrod2",lwd=3,lty=4)
legend("topleft",
       c("Naive OLS", "2SLS", "Lasso","Post-Lasso","PCA"),
       col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
       lty = c(1,1,1,2,4),lwd=c(2,3,2,2,3))
lines(rep(1,17),seq(0,16,by=1),lwd=2)

# for checking Figure 21
boxplot(Coef601,col = c("paleturquoise4","red","blue","darkgreen","darkgoldenrod2"),
        lwd=2,main="Beta  (Figure 21.2)")
abline(1,0,col="darkred",lwd=2)

