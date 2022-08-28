rm(list=ls())
library(MASS)
set.seed(1)
n<-100
p<-2
mu<-c(0,10)
sigma<-matrix(c(2,1,1,1.5),nrow=p,ncol=p,byrow = TRUE)
beta_range<-seq(0,0.1,by=0.001)
beta<-sample(beta_range,p)
e.mu=0
e.sigma=1
dgp = function(n, beta, mu, sigma){
  x<-mvrnorm(n,mu,sigma)
  error<-rnorm(n,e.mu,e.sigma)
  y <- x%*% beta + error
  data<-data.frame(y,x)
  return(data)
}
data.train<-dgp(n,beta,mu,sigma)
x1.train<-data.train[,2]
x2.train<-data.train[,3]
x.train<-cbind(x1.train,x2.train)
y.train<-data.train[,1]
data.test<-dgp(n,beta,mu,sigma)
x1.test<-data.test[,2]
x2.test<-data.test[,3]
x.test<-cbind(x1.test,x2.test)
y.test<-data.test[,1]
####for the question a
library(pls)
pcr.fit<-pcr(y.train~x.train,data=data.train,
             scale=TRUE,validation="CV")
first.com.cof<-c(pcr.fit$loadings[1,1],pcr.fit$loadings[2,1])
second.com.cof<-c(pcr.fit$loadings[1,2],pcr.fit$loadings[2,2])
x1.train.sd<-(x1.train)/sd(x1.train)
x2.train.sd<-(x2.train)/sd(x2.train)
x.train.sd <- cbind(x1.train.sd,x2.train.sd)
first.com<-first.com.cof[1]*(x1.train.sd-mean(x1.train.sd))+
  first.com.cof[2]*(x2.train.sd-mean(x2.train.sd))
second.com<-second.com.cof[1]*(x1.train.sd-mean(x1.train.sd))+
  second.com.cof[2]*(x2.train.sd-mean(x2.train.sd))
###another  way to verify the answer 
cbind(first.com,second.com,pcr.fit$scores)###and they are same 
###visualize the components
eigen_vectors<-eigen(cov(x.train))$vectors
x1.train.mean <- x1.train-mean(x1.train)
x2.train.mean <- x2.train-mean(x2.train)
plot(x1.train.mean,x2.train.mean,asp=1)
abline(a=0,b=(eigen_vectors[1,1] /eigen_vectors[2,1]),col="darkgreen")
abline(a=0,b=(eigen_vectors[1,2] /eigen_vectors[2,2]),col="red")
#####another way
# std data 
ei <- eigen(cov(x.train.sd))
ei <- ei$vectors
pl <- seq(-10,10,1)
ply <- pl*ei[2,1]/ei[1,1]+mean(x2.train.sd)
pl2 <- pl
ply2 <- pl2*ei[2,2]/ei[1,2]+mean(x2.train.sd)
###plot
plot(x1.train.sd,x2.train.sd,asp=1,main = "Standardized Data")
lines(pl,ply,col="darkgreen",lwd=2,xlim=c(-10,10))
lines(pl2,ply2,col='red',lty=2,xlim=c(-10,10))

### original data # wrong 
ei <- eigen(cov(x.train))
ei <- ei$vectors
pl <- seq(-10,10,1)
ply <- pl*ei[2,1]/ei[1,1]+mean(x2.train)
pl2 <- pl
ply2 <- pl2*ei[2,2]/ei[1,2]+mean(x2.train)
###plot
plot(x1.train,x2.train,asp=1,main = "original data")
lines(pl,ply,col="darkgreen",lwd=2,xlim=c(-10,10))
lines(pl2,ply2,col='red',lty=2,xlim=c(-10,10))

################################
plot(first.com,x1.train.sd,pch=15,col="red",
     main = "First Component v.s. X1")
plot(first.com,x2.train.sd,pch=16,col="red",
     main = "First Component v.s. X2")
plot(second.com,x1.train.sd,pch=17,col="red",
     main = "Second Component v.s. X1")
plot(second.com,x2.train.sd,pch=18,col="red",
     main = "Second Component v.s. X2")
#######################b
summary(pcr.fit)
validationplot(pcr.fit,val.type="MSEP")
###perform pcr on training data and evaluate test set performance
pcr.pred<-predict(pcr.fit,x.test,ncomp = 1)
##compute the test mse
mean((pcr.pred-y.test)^2)
###another one 
pcr.pred<-predict(pcr.fit,x.test,ncomp = 2)
##compute the test mse
mean((pcr.pred-y.test)^2)
####################################c
###############################for the ridge regression
library(glmnet)
I<-diag(1,p)
grid<-10^seq(-1,3,length=100)
cv.out<-cv.glmnet(x.train,y.train,lambda = grid,alpha=0,intercept=FALSE)
best_lambda<-cv.out$lambda.min
beta.hat.ridge<-solve((t(x.train)%*%x.train)+best_lambda*I)%*%t(x.train)%*%y.train
ri.pred<-x.test%*%beta.hat.ridge
mean((ri.pred-y.test)^2)
##############################for the lasso 
lasso.mod<-glmnet(x.train,y.train, alpha=1, lambda = grid,intercept=F)
##### we now perform cross-validation and compute the associated prediction error
cv.out<-cv.glmnet(x.train, y.train,alpha=1,intercept=F)
bestlambda<-cv.out$lambda.min
lasso.pred<-predict(lasso.mod, s=bestlambda, newx = x.test)
mean((lasso.pred-y.test)^2)
lasss <- glmnet(x.train,y.train, alpha=1, lambda = bestlambda,intercept=F)
summary(lasss)

abc <- ifelse(lasss$beta==0,0,1)
abc
####################### for PCR
###perform pcr on training data and evaluate test set performance
pcr.pred1<-predict(pcr.fit,x.test,ncomp = 1)
##compute the test mse
mean((pcr.pred1-y.test)^2)
###perform pcr on training data and evaluate test set performance
pcr.pred2<-predict(pcr.fit,x.test,ncomp = 2)
##compute the test mse
mean((pcr.pred2-y.test)^2)



#########################################################
#########################################################
##exercise 2
T<-100
p2<-10
I<-diag(1,p2,p2)
beta.hat.ols<-matrix(NA,nrow=T,ncol=p2)
beta.hat.ridge<-matrix(NA,nrow=T,ncol=p2)
mse.ols<-c()
mse.ri<-c()
mse.lasso<-c()
mse.pcr<-c()
y.pre<-matrix(NA,nrow=T,ncol=n)
y.pre.ri<-matrix(NA,nrow=T,ncol=n)
y.pre.lasso<-matrix(NA,nrow=T,ncol=n)
y.pre.pcr<-matrix(NA,nrow=T,ncol=n)
##for the sake of convenience,we design another dgp 
mu<-c(0,10,9,8,7,6,5,4,3,2)
sigma<-matrix(c(2,rep(0.01,10),2,rep(0.01,10),3,rep(0.01,10),3,
                  rep(0.01,10),1,rep(0.01,10),1,rep(0.01,10),2,rep(0.01,10),
                  1,rep(0.01,10),2,rep(0.01,10),1),nrow=p2,ncol=p2,byrow = TRUE)
sigma
beta.true<-sample(beta_range,10)
beta.true
e.mu<-0
e.sigma<-1
dgp2 = function(n, beta.true, mu, sigma){
  X<-mvrnorm(n,mu,sigma)
  error<-rnorm(n,e.mu,e.sigma)
  y <- X%*% beta.true + error
  data<-data.frame(y,X)
  return(data)
}
grid <- 10^seq(-1,3,length=100)
###the begin to simulation
set.seed(527)
for (i in 1:T){
  train.data<-dgp2(n,beta.true,mu,sigma)
  test.data<-dgp2(n,beta.true,mu,sigma)
  x.train<-cbind(train.data[,2],train.data[,3],train.data[,4],train.data[,5],train.data[,6],train.data[,7],train.data[,8],train.data[,9],train.data[,10],train.data[,11])
  x.test<-cbind(test.data[,2],test.data[,3],test.data[,4],test.data[,5],test.data[,6],test.data[,7],test.data[,8],test.data[,9],test.data[,10],test.data[,11])
  y.train<-train.data$y
  y.test<-test.data$y
  beta.hat.ols[i,]<-solve(t(x.train)%*%x.train)%*%t(x.train)%*%y.train
  y.pre[i,]<-x.test%*%(beta.hat.ols[i,])
  mse.ols[i]<-mean((y.pre[i,]-y.test)^2)
  ###for ridge regression
  cv.out1<-cv.glmnet(x.train,y.train,lambda = grid,alpha=0,intercept=FALSE)
  best_lambda1<-cv.out1$lambda.min
  beta.hat.ridge[i,]<-solve((t(x.train)%*%x.train)+best_lambda1*I)%*%t(x.train)%*%y.train
  y.pre.ri[i,]<-x.test%*%beta.hat.ridge[i,]
  mse.ri[i]<-mean((y.pre.ri[i,]-y.test)^2)
  ###for lasso
  lasso.mod<-glmnet(x.train,y.train, alpha=1, lambda = grid,intercept=F)
  cv.out<-cv.glmnet(x.train, y.train,alpha=1,intercept=F)
  bestlambda<-cv.out$lambda.min
  y.pre.lasso[i,]<-predict(lasso.mod, s=bestlambda, newx = x.test)
  mse.lasso[i]<-mean((y.pre.lasso[i,]-y.test)^2)
  ####for PCR 
  pcr.fit<-pcr(y.train~x.train,data=train.data,scale=TRUE,validation="CV")
  cverr <- RMSEP(pcr.fit)$val[1,,-1]
  imin<-which.min(cverr)
  y.pre.pcr[i,]<-predict(pcr.fit,x.test,ncomp = imin)
  mse.pcr[i]<-mean((y.pre.pcr[i,]-y.test)^2)
}
data.box<-cbind(mse.ols,mse.ri,mse.lasso,mse.pcr)
boxplot(data.box, names=c("OLS","Ridge","Lasso","PCR"),
        main="MSE")
abline(h=1,col="darkred")##### variance error term
######################################b
######first one :change beta
beta_range1<-seq(2,3,by=0.01)
beta.true1<-sample(beta_range1,10)
for (i in 1:T){
  train.data<-dgp2(n,beta.true1,mu,sigma)
  test.data<-dgp2(n,beta.true1,mu,sigma)
  x.train<-cbind(train.data[,2],train.data[,3],train.data[,4],train.data[,5],train.data[,6],train.data[,7],train.data[,8],train.data[,9],train.data[,10],train.data[,11])
  x.test<-cbind(test.data[,2],test.data[,3],test.data[,4],test.data[,5],test.data[,6],test.data[,7],test.data[,8],test.data[,9],test.data[,10],test.data[,11])
  y.train<-train.data$y
  y.test<-test.data$y
  beta.hat.ols[i,]<-solve(t(x.train)%*%x.train)%*%t(x.train)%*%y.train
  y.pre[i,]<-x.test%*%(beta.hat.ols[i,])
  mse.ols[i]<-mean((y.pre[i,]-y.test)^2)
  ###for ridge regression
  cv.out1<-cv.glmnet(x.train,y.train,lambda = grid,alpha=0,intercept=FALSE)
  best_lambda1<-cv.out1$lambda.min
  beta.hat.ridge[i,]<-solve((t(x.train)%*%x.train)+best_lambda1*I)%*%t(x.train)%*%y.train
  y.pre.ri[i,]<-x.test%*%beta.hat.ridge[i,]
  mse.ri[i]<-mean((y.pre.ri[i,]-y.test)^2)
  ###for lasso
  lasso.mod<-glmnet(x.train,y.train, alpha=1, lambda = grid,intercept=F)
  cv.out<-cv.glmnet(x.train, y.train,alpha=1,intercept=F)
  bestlambda<-cv.out$lambda.min
  y.pre.lasso[i,]<-predict(lasso.mod, s=bestlambda, newx = x.test)
  mse.lasso[i]<-mean((y.pre.lasso[i,]-y.test)^2)
  ####for PCR 
  pcr.fit<-pcr(y.train~x.train,data=train.data,scale=TRUE,validation="CV")
  cverr <- RMSEP(pcr.fit)$val[1,,-1]
  imin <- which.min(cverr)
  y.pre.pcr[i,]<-predict(pcr.fit,x.test,ncomp = imin)
  mse.pcr[i]<-mean((y.pre.pcr[i,]-y.test)^2)
}
data.box2<-cbind(mse.ols,mse.ri,mse.lasso,mse.pcr)
boxplot(data.box2, names=c("OLS","Ridge","Lasso","PCR"),
        main="MSE")
boxplot(data.box2[,-3], names=c("Ridge","Lasso","PCR"),
        main="MSE")
abline(h=1,col="darkred")##### variance error term
#########second way:change sigma
####construct a function to define the sigma matrix
si<-function(s){
  for (i in 1:10){
    sm[i,i]<-s[i]
  }
  for (i in 1:10){
    for (j in 1:10){
      if(j!=i){
        sm[i,j]<-sqrt(sm[i,i]*sm[j,j])*0.999999
        sm[j,i]<-sm[i,j]
      }
    }
  }
  return(sm)
}
sm<-matrix(NA,10,10)
s<-c(1,1,1,1,1,1,1,1,1,1)
si(s)
sigma<-si(s)
###simulation
for (i in 1:T){
  train.data<-dgp2(n,beta.true,mu,sigma)
  test.data<-dgp2(n,beta.true,mu,sigma)
  x.train<-cbind(train.data[,2],train.data[,3],train.data[,4],train.data[,5],train.data[,6],train.data[,7],train.data[,8],train.data[,9],train.data[,10],train.data[,11])
  x.test<-cbind(test.data[,2],test.data[,3],test.data[,4],test.data[,5],test.data[,6],test.data[,7],test.data[,8],test.data[,9],test.data[,10],test.data[,11])
  y.train<-train.data$y
  y.test<-test.data$y
  beta.hat.ols[i,]<-solve(t(x.train)%*%x.train)%*%t(x.train)%*%y.train
  y.pre[i,]<-x.test%*%(beta.hat.ols[i,])
  mse.ols[i]<-mean((y.pre[i,]-y.test)^2)
  ###for ridge regression
  cv.out1<-cv.glmnet(x.train,y.train,lambda = grid,alpha=0,intercept=FALSE)
  best_lambda1<-cv.out1$lambda.min
  beta.hat.ridge[i,]<-solve((t(x.train)%*%x.train)+best_lambda1*I)%*%t(x.train)%*%y.train
  y.pre.ri[i,]<-x.test%*%beta.hat.ridge[i,]
  mse.ri[i]<-mean((y.pre.ri[i,]-y.test)^2)
  ###for lasso
  lasso.mod<-glmnet(x.train,y.train, alpha=1, lambda = grid,intercept=F)
  cv.out<-cv.glmnet(x.train, y.train,alpha=1,intercept=F)
  bestlambda<-cv.out$lambda.min
  y.pre.lasso[i,]<-predict(lasso.mod, s=bestlambda, newx = x.test)
  mse.lasso[i]<-mean((y.pre.lasso[i,]-y.test)^2)
  ####for PCR 
  pcr.fit<-pcr(y.train~x.train,data=train.data,scale=TRUE,validation="CV")
  cverr <- RMSEP(pcr.fit)$val[1,,-1]
  imin <- which.min(cverr)
  y.pre.pcr[i,]<-predict(pcr.fit,x.test,ncomp = imin)
  mse.pcr[i]<-mean((y.pre.pcr[i,]-y.test)^2)
}
data.box4<-cbind(mse.ols,mse.ri,mse.lasso,mse.pcr)
boxplot(data.box4, names=c("OLS","Ridge","Lasso","PCR"),
        main="MSE")
abline(h=1,col="darkred")##### variance error term
# lasso for model selection in high correlated variables
# case may not play well, but for prediction it does 
# structural way: find how lasso and ridge shrink 
# parameters, it's about prediction and model selection
# 
# also could try lower cov 



