#install.packages("randomForest")
#install.packages("gbm")
#install.packages("caret")
#rm(list=ls())
library(randomForest)
library(caret)
library(gbm)
library(MASS)
library(tree)


### 






p <- 5
n <- 500
mu <- rep(0,p)
sigma <- matrix(c(2,0.1,0.1,0.1,0.1,
                  0.1,2.5,0.1,0.1,0.1,
                  0.1,0.1,3,0.1,0.1,
                  0.1,0.1,0.1,1.9,0.1,
                  0.1,0.1,0.1,0.1,2.3),nrow = p,byrow = T)

dgp <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-53*(x[,1]*x[,2]>=2) + 2.24^x[,3] + 1.05^(x[,4]*x[,5])+
  -exp(-abs(x[,4]))-25*(x[,3]/x[,4]>=1) + exp(x[,5]) +
  exp(x[,1])*2 +1/x[,2] +error
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}
set.seed(527)
data.train <- dgp(mu,sigma)
summary(data.train)
data.test <- dgp(mu,sigma)
summary(data.test)
# check conditional linearity
par(mfrow=c(2,3))
plot(data.train$X1,data.train$y,col=data.train$c)
plot(data.train$X2,data.train$y,col=data.train$c)
plot(data.train$X3,data.train$y,col=data.train$c)
plot(data.train$X4,data.train$y,col=data.train$c)
plot(data.train$X5,data.train$y,col=data.train$c)
hist(data.train$y)
par(mfrow=c(3,4))
plot(data.train$X1,data.train$X2,col=data.train$c)
plot(data.train$X1,data.train$X3,col=data.train$c)
plot(data.train$X1,data.train$X4,col=data.train$c)
plot(data.train$X1,data.train$X5,col=data.train$c)
plot(data.train$X2,data.train$X3,col=data.train$c)
plot(data.train$X2,data.train$X4,col=data.train$c)
plot(data.train$X2,data.train$X5,col=data.train$c)
plot(data.train$X3,data.train$X4,col=data.train$c)
plot(data.train$X3,data.train$X5,col=data.train$c)
plot(data.train$X4,data.train$X5,col=data.train$c)

## pruned tree
pvst <- matrix(NA,100,2)
colnames(pvst) <- c("Pruned","Tree")
set.seed(123)
for (i in 1:100) {
  data.train <- dgp(mu,sigma)
  data.test <- dgp(mu,sigma)
treebasic <- tree(c~.-y,data = data.train)
#summary(treebasic)
#plot(treebasic)
#text(treebasic,pretty = 0)
#treebasic
cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
#cv.tessbasic
#plot(cv.tessbasic$size , cv.tessbasic$dev, type = "b")
pruned <- prune.misclass(treebasic,best =
      cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
#plot(pruned)
#text(pruned,pretty=0,digits = 1)
c_tree_orig <- predict(treebasic,data.test,
                       type = "class")
c_tree_pruned <- predict(pruned,newdata = data.test,
                         type = "class")

a <- confusionMatrix(c_tree_pruned,data.test$c)
b <- confusionMatrix(c_tree_orig,data.test$c)
pvst[i,1] <- a$overall[1]
pvst[i,2] <- b$overall[1]
}
boxplot(pvst)
par(mfrow=c(2,1))
hist(pvst[,1],ylim = c(0,30))
hist(pvst[,2])
plot(pvst[,1],pvst[,2])
abline(0,1,col='red')
mean(pvst[,1]>pvst[,2])
mean(pvst[,1]==pvst[,2])
mean(pvst[,1]<pvst[,2])
## boosted
# 10 fold CV
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
# grid for CV
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.05,
                        n.minobsinnode = 10)
# the minimum number of training set samples in
# a node to commence splitting (n.minobsinnode)
set.seed(527)
gbmFit2 <- train(c ~.-y, data = data.train, 
                 method = "gbm", 
                 trControl = trControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 distribution = "bernoulli")
gbmFit2
summary(gbmFit2)
gbmFit2$bestTune$n.trees
plot(gbmFit2)
#boosted <- gbmFit2$finalModel
#summary(boosted)
c_boosted <- predict(gbmFit2,newdata=data.test)
confusionMatrix(c_boosted,data.test$c)

# alternative not works well
#set.seed(527)
#contrasts(data.train$c)
#C <- ifelse(data.train$c=="high",0,1)
#boosted2 <- gbm(c~.-y -c,data=data.train,
#               distribution = "bernoulli",
#               verbose = FALSE,
#               shrinkage = 0.05,
#               n.trees = gbmFit2$bestTune$n.trees,
#               n.minobsinnode = 10,
# interaction.depth =gbmFit2$bestTune$interaction.depth)
#summary(boosted2)
#plot(boosted2)
#c_boosted2 <- predict(boosted2,data.test,n.trees = 
#                        gbmFit2$bestTune$n.trees,"response")
#c_boosted2 <- 
#confusionMatrix(c_boosted2,data.test$c)

################################################
### exc 2
# a
p <- 5
mu <- rep(0,p)
sigma <- matrix(c(2,0.1,0.1,0.1,0.1,
                  0.1,2.5,0.1,0.1,0.1,
                  0.1,0.1,3,0.1,0.1,
                  0.1,0.1,0.1,1.9,0.1,
                  0.1,0.1,0.1,0.1,2.3),nrow = p,byrow = T)

dgp2a <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-x[,1]+2*x[,2]+1.5*x[,1]*x[,2]+x[,3]^2+
    1/(abs(x[,4])+0.2)+log(abs(x[,5])+1)+0.5^x[,4]*x[,5]
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}

dgp2control <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-1.23*ifelse(x[,1]>=0.49,1.5,-1.5)+
    2*ifelse(x[,2]<=0.1,2,-2)-0.8*x[,3]+
    1.9*x[,4]+2*x[,5]+0.5*x[,4]*x[,5]
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}
set.seed(527)
data.train <- dgp2control(mu,sigma)
summary(data.train)
data.test <- dgp2a(mu,sigma)
summary(data.test)
par(mfrow=c(2,3))
plot(data.train$X1,data.train$y,col=data.train$c)
plot(data.train$X2,data.train$y,col=data.train$c)
plot(data.train$X3,data.train$y,col=data.train$c)
plot(data.train$X4,data.train$y,col=data.train$c)
plot(data.train$X5,data.train$y,col=data.train$c)
hist(data.train$y)
par(mfrow=c(3,4))
plot(data.train$X1,data.train$X2,col=data.train$c)
plot(data.train$X1,data.train$X3,col=data.train$c)
plot(data.train$X1,data.train$X4,col=data.train$c)
plot(data.train$X1,data.train$X5,col=data.train$c)
plot(data.train$X2,data.train$X3,col=data.train$c)
plot(data.train$X2,data.train$X4,col=data.train$c)
plot(data.train$X2,data.train$X5,col=data.train$c)
plot(data.train$X3,data.train$X4,col=data.train$c)
plot(data.train$X3,data.train$X5,col=data.train$c)
plot(data.train$X4,data.train$X5,col=data.train$c)

treebasic <- tree(c~.-y,data = data.train)
plot(treebasic)
text(treebasic,pretty = 0)
cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
pruned <- prune.misclass(treebasic,best =
  cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
c_tree_pruned <- predict(pruned,newdata = data.test,
                         type = "class")
confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
accuracy1 <- confusion1$overall[1]

## boosted
# 10 fold CV
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
# grid for CV
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:20)*100, 
                        shrinkage = 0.05,
                        n.minobsinnode = 10)
# the minimum number of training set samples in
# a node to commence splitting (n.minobsinnode)
set.seed(527)
gbmFit2 <- train(c ~.-y, data = data.train, 
                 method = "gbm", 
                 trControl = trControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 distribution = "bernoulli")
plot(gbmFit2)
summary(gbmFit2)
c_boosted <- predict(gbmFit2,newdata=data.test)
confusion2 <- confusionMatrix(c_boosted,data.test$c)
accuracy2 <- confusion2$overall[1]

# random forests
tuneGrid <- expand.grid(.mtry = c(1: 5))
set.seed(123)
rf_mtry <- train(c~.-y,
                 data = data.train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE)
rf_mtry$finalModel
summary(rf_mtry)
plot(rf_mtry$finalModel)
plot(rf_mtry)
importance(rf_mtry$finalModel)
varImpPlot(rf_mtry$finalModel)

c_rf <- predict(rf_mtry$finalModel,data.test)
confusion3 <- confusionMatrix(c_rf,data.test$c)
accuracy3 <- confusion3$overall[1]

# control group 
r <- 50
set.seed(555)
accuracy.control <- matrix(NA,r,3)
colnames(accuracy.control) <- c("Pruned","Boosted","Random Forest")
for (i in 1:r) {
  data.train <- dgp2control(mu,sigma)
  data.test <- dgp2control(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy.control[i,1] <- confusion1$overall[1]
  gbmFit2 <- train(c ~.-y, data = data.train, 
                   method = "gbm", 
                   trControl = trControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid,
                   distribution = "bernoulli")
  c_boosted <- predict(gbmFit2,newdata=data.test)
  confusion2 <- confusionMatrix(c_boosted,data.test$c)
  accuracy.control[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy.control[i,3] <- confusion3$overall[1]
}
boxplot(accuracy.control,main="Accuracy control")
plot(accuracy.control[,2],accuracy.control[,3],xlab = "Accuracy of Boosting",
     ylab = "Accuracy of Random Forest", main = "Boosting v.s. RF")
abline(a=0, b=1, col='red')
covvar <- cor(data.train[,c(3:7)])
cov(data.train[,3],data.train$c)

# change to a non linear function 
r <- 50
set.seed(555)
accuracy <- matrix(NA,r,3)
colnames(accuracy) <- c("Pruned","Boosted","Random Forest")
for (i in 1:r) {
  data.train <- dgp2a(mu,sigma)
  data.test <- dgp2a(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy[i,1] <- confusion1$overall[1]
  gbmFit2 <- train(c ~.-y, data = data.train, 
                   method = "gbm", 
                   trControl = trControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid,
                   distribution = "bernoulli")
  c_boosted <- predict(gbmFit2,newdata=data.test)
  confusion2 <- confusionMatrix(c_boosted,data.test$c)
  accuracy[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy[i,3] <- confusion3$overall[1]
}
boxplot(accuracy,main="Accuracy")
plot(accuracy[,2],accuracy[,3],xlab = "Accuracy of Boosting",
     ylab = "Accuracy of Random Forest", main = "Boosting v.s. RF")
abline(a=0, b=1, col='red')
covvar <- cor(data.train[,c(3:7)])
cov(data.train[,3],data.train$c)



## high covariance case 

p <- 5
mu <- rep(0,p)
sigma2 <- matrix(c(2,1.8,0.1,0.1,0.1,
                  1.8,2.5,0.1,0.1,0.1,
                  0.1,0.1,3,0.1,0.1,
                  0.1,0.1,0.1,1.9,1.5,
                  0.1,0.1,0.1,1.5,2.3),nrow = p,byrow = T)

dgp2a <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-x[,1]+2*x[,2]+1.5*x[,1]*x[,2]+x[,3]^2+
    1/(abs(x[,4])+0.2)+log(abs(x[,5])+1)+0.5^x[,4]*x[,5]
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}
set.seed(527)
data.train <- dgp2a(mu,sigma2)
summary(data.train)
data.test <- dgp2a(mu,sigma2)
summary(data.test)
par(mfrow=c(2,3))
plot(data.train$X1,data.train$y,col=data.train$c)
plot(data.train$X2,data.train$y,col=data.train$c)
plot(data.train$X3,data.train$y,col=data.train$c)
plot(data.train$X4,data.train$y,col=data.train$c)
plot(data.train$X5,data.train$y,col=data.train$c)
hist(data.train$y)
par(mfrow=c(3,4))
plot(data.train$X1,data.train$X2,col=data.train$c)
plot(data.train$X1,data.train$X3,col=data.train$c)
plot(data.train$X1,data.train$X4,col=data.train$c)
plot(data.train$X1,data.train$X5,col=data.train$c)
plot(data.train$X2,data.train$X3,col=data.train$c)
plot(data.train$X2,data.train$X4,col=data.train$c)
plot(data.train$X2,data.train$X5,col=data.train$c)
plot(data.train$X3,data.train$X4,col=data.train$c)
plot(data.train$X3,data.train$X5,col=data.train$c)
plot(data.train$X4,data.train$X5,col=data.train$c)


trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
# grid for CV
gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:20)*100, 
                        shrinkage = 0.05,
                        n.minobsinnode = 10)

r <- 50
set.seed(555)
accuracy2 <- matrix(NA,r,3)
colnames(accuracy2) <- c("Pruned","Boosted","Random Forest")
for (i in 1:r) {
  data.train <- dgp2a(mu,sigma2)
  data.test <- dgp2a(mu,sigma2)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy2[i,1] <- confusion1$overall[1]
  gbmFit2 <- train(c ~.-y, data = data.train, 
                   method = "gbm", 
                   trControl = trControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid,
                   distribution = "bernoulli")
  c_boosted <- predict(gbmFit2,newdata=data.test)
  confusion2 <- confusionMatrix(c_boosted,data.test$c)
  accuracy2[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy2[i,3] <- confusion3$overall[1]
}
boxplot(accuracy2,main="Accuracy2")
plot(accuracy2[,2],accuracy2[,3],xlab = "Accuracy of Boosting",
     ylab = "Accuracy of Random Forest", main = "Boosting v.s. RF")
abline(a=0, b=1, col='red')
covvar <- cor(data.train[,c(3:7)])
cov(data.train[,3],data.train$c)


aaaa <- cbind(accuracy.control,accuracy,accuracy2)
boxplot(aaaa,main="Accuracy",col=c(
  "red","blue","darkgreen","red","blue","darkgreen",
  "red","blue","darkgreen"
),las=1,horizontal=T,names=c("Control","Control","Control",
"nonlinear","nonlinear","nonlinear",
"correlate","correlate","correlate"
),cex=1,lwd=2)
legend("topleft",legend = c("Pruned","Boosted","Random Forest"),
   col=c("red","blue","darkgreen"),pch = c(15,15,15))
## for highly correlated predictors, all trees perform
# better, since at large sample size, trees vary less
# and random forest performs better than boosting, 
# because r.f. has lower variance, 

dgp2combine <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-1.23*ifelse(x[,1]>=0.49,1.5,-1.5)+
    2*ifelse(x[,2]<=0.1,2,-2)-0.8*x[,3]+
    1.9*x[,4]+2*x[,5]+0.5*x[,4]*x[,5]+0.6*
    x[,1]*x[,2]+5*x[,2]*x[,3]*x[,4]+
    0.4*x[,1]*x[,5]
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}

r <- 50
set.seed(555)
accuracy3 <- matrix(NA,r,3)
colnames(accuracy3) <- c("Pruned","Boosted","Random Forest")
for (i in 1:r) {
  data.train <- dgp2combine(mu,sigma)
  data.test <- dgp2combine(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy3[i,1] <- confusion1$overall[1]
  gbmFit2 <- train(c ~.-y, data = data.train, 
                   method = "gbm", 
                   trControl = trControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid,
                   distribution = "bernoulli")
  c_boosted <- predict(gbmFit2,newdata=data.test)
  confusion2 <- confusionMatrix(c_boosted,data.test$c)
  accuracy3[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy3[i,3] <- confusion3$overall[1]
}
boxplot(accuracy3,main="Accuracy3")
plot(accuracy3[,2],accuracy3[,3],xlab = "Accuracy of Boosting",
     ylab = "Accuracy of Random Forest", main = "Boosting v.s. RF")
abline(a=0, b=1, col='red')

aaaa <- cbind(accuracy.control,accuracy,accuracy2,accuracy3)
boxplot(aaaa,main="Accuracy",col=c(
  "red","blue","darkgreen","red","blue","darkgreen",
  "red","blue","darkgreen","red","blue","darkgreen"
),las=1,horizontal=T,names=c("Control","Control","Control",
     "nonlinear","nonlinear","nonlinear",
"correlate","correlate","correlate","Combine","Combine","Combine"
),cex=1,lwd=2)
legend("topleft",legend = c("Pruned","Boosted","Random Forest"),
       col=c("red","blue","darkgreen"),pch = c(15,15,15))
## we can see that the combined case, boosted
# works at least not as bad as r.f., since it 
# could recognize spillover effect of some 
# weak predictors. Meanwhile, because boosting
# learns slowly, it could find the combination 
# effects finally, but for r.f. it has chances
# that not find such effect. 


r <- 50
set.seed(555)
accuracy.rf_bag <- matrix(NA,r,3)
colnames(accuracy.rf_bag) <- c("Pruned","Bagging","Random Forest")
for (i in 1:r) {
  data.train <- dgp2combine(mu,sigma)
  data.test <- dgp2combine(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy.rf_bag[i,1] <- confusion1$overall[1]
 bag <- randomForest(c~. -y, data = data.train,
                     mtry=5)
 predbag <- predict(bag,data.test)
 confusion2 <- confusionMatrix(predbag,data.test$c)
 accuracy.rf_bag[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy.rf_bag[i,3] <- confusion3$overall[1]
}
boxplot(accuracy.rf_bag,main="Accuracy",lwd=2)
plot(accuracy.rf_bag[,2],accuracy.rf_bag[,3],xlab = "Accuracy of Bagging",
     ylab = "Accuracy of Random Forest", main = "Bagging v.s. RF")
abline(a=0, b=1, col='red')
mean(accuracy.rf_bag[,3]==accuracy.rf_bag[,2])
mean(accuracy.rf_bag[,3]>accuracy.rf_bag[,2])
mean(accuracy.rf_bag[,3]<accuracy.rf_bag[,2])

#######
# control group 
r <- 50
set.seed(555)
accuracy.rf_bag2 <- matrix(NA,r,3)
colnames(accuracy.rf_bag2) <- c("Pruned","Bagging","Random Forest")
for (i in 1:r) {
  data.train <- dgp2a(mu,sigma)
  data.test <- dgp2a(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy.rf_bag2[i,1] <- confusion1$overall[1]
  bag <- randomForest(c~. -y, data = data.train,
                      mtry=5)
  predbag <- predict(bag,data.test)
  confusion2 <- confusionMatrix(predbag,data.test$c)
  accuracy.rf_bag2[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy.rf_bag2[i,3] <- confusion3$overall[1]
}
boxplot(accuracy.rf_bag2,main="Accuracy",lwd=2)
plot(accuracy.rf_bag2[,2],accuracy.rf_bag2[,3],xlab = "Accuracy of Bagging",
     ylab = "Accuracy of Random Forest", main = "Bagging v.s. RF")
abline(a=0, b=1, col='red')
mean(accuracy.rf_bag2[,3]==accuracy.rf_bag2[,2])
mean(accuracy.rf_bag2[,3]>accuracy.rf_bag2[,2])
mean(accuracy.rf_bag2[,3]<accuracy.rf_bag2[,2])

########################################
######################################3#
## with leading predictor x1
dgp2leading <- function(mu,sigma){
  error <- rnorm(n,0,0.5)
  x <- mvrnorm(n,mu,Sigma = sigma)
  y <-20*x[,1]+2*x[,2]+1.5*x[,1]*x[,2]+2*x[,3]^3+
    1/(abs(x[,4])+0.2)+log(abs(x[,5])+1)+0.5^x[,4]*x[,5]
  q <- quantile(y,c(0.2,0.8),names = F)
  c <- rep("high",n)
  low.j <- y<q[1] | y>q[2]
  c[low.j] <- "low"
  c <- as.factor(c)
  data <- data.frame(c,y,x)
  return(data)
}
set.seed(527)
data.train <- dgp2leading(mu,sigma)
summary(data.train)
data.test <- dgp2leading(mu,sigma)
summary(data.test)
plot(data.train$X1,data.train$c,col=data.train$c)
plot(data.train$X2,data.train$c,col=data.train$c)
plot(data.train$X3,data.train$c,col=data.train$c)
plot(data.train$X4,data.train$c,col=data.train$c)
plot(data.train$X5,data.train$c,col=data.train$c)
hist(data.train$y)
plot(data.train$X1,data.train$X2,col=data.train$c)
plot(data.train$X1,data.train$X3,col=data.train$c)
plot(data.train$X1,data.train$X4,col=data.train$c)
plot(data.train$X1,data.train$X5,col=data.train$c)
plot(data.train$X2,data.train$X3,col=data.train$c)
plot(data.train$X2,data.train$X4,col=data.train$c)
plot(data.train$X2,data.train$X5,col=data.train$c)
plot(data.train$X3,data.train$X4,col=data.train$c)
plot(data.train$X3,data.train$X5,col=data.train$c)
plot(data.train$X4,data.train$X5,col=data.train$c)

r <- 50
set.seed(555)
accuracy.rf_bag3 <- matrix(NA,r,3)
colnames(accuracy.rf_bag3) <- c("Pruned","Bagging","Random Forest")
for (i in 1:r) {
  data.train <- dgp2leading(mu,sigma)
  data.test <- dgp2leading(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy.rf_bag3[i,1] <- confusion1$overall[1]
  bag <- randomForest(c~. -y, data = data.train,
                      mtry=5)
  predbag <- predict(bag,data.test)
  confusion2 <- confusionMatrix(predbag,data.test$c)
  accuracy.rf_bag3[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy.rf_bag3[i,3] <- confusion3$overall[1]
}
boxplot(accuracy.rf_bag3,main="Accuracy",lwd=2)
plot(accuracy.rf_bag3[,2],accuracy.rf_bag3[,3],xlab = "Accuracy of Bagging",
     ylab = "Accuracy of Random Forest", main = "Bagging v.s. RF")
abline(a=0, b=1, col='red')
mean(accuracy.rf_bag3[,3]==accuracy.rf_bag3[,2])
mean(accuracy.rf_bag3[,3]>accuracy.rf_bag3[,2])
mean(accuracy.rf_bag3[,3]<accuracy.rf_bag3[,2])

########################################
# compare bagging and rf
brf <- cbind(accuracy.rf_bag,accuracy.rf_bag2,accuracy.rf_bag3)
boxplot(brf,names=c("combine","combine","combine",
                    "control","control","control",
                    "leading","leading","leading"),
        col=c("red","brown","darkgreen","red","brown","darkgreen",
              "red","brown","darkgreen"),lwd=2,
       las=1, horizontal=T,main="Bagging v.s. R.F.")
legend("topleft",legend = c("Pruned","Bagging","Random Forest"),
       col=c("red","brown","darkgreen"),pch = c(15,15,15))

######################################## 





## finally compare all with leading case

r <- 50
set.seed(555)
accuracy.all.leading <- matrix(NA,r,4)
colnames(accuracy.all.leading) <- c("Pruned","Boosted","Random Forest","Bagging")
gbmGrid2 <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:15)*150, 
                        shrinkage = 0.05,
                        n.minobsinnode = 10)
for (i in 1:r) {
  data.train <- dgp2leading(mu,sigma)
  data.test <- dgp2leading(mu,sigma)
  treebasic <- tree(c~.-y,data = data.train)
  cv.tessbasic <- cv.tree(treebasic,FUN = prune.misclass)
  pruned <- prune.misclass(treebasic,best =
                             cv.tessbasic$size[which.min(cv.tessbasic$dev)] )
  c_tree_pruned <- predict(pruned,newdata = data.test,
                           type = "class")
  confusion1 <- confusionMatrix(c_tree_pruned,data.test$c)
  accuracy.all.leading[i,1] <- confusion1$overall[1]
  gbmFit2 <- train(c ~.-y, data = data.train, 
                   method = "gbm", 
                   trControl = trControl, 
                   verbose = FALSE, 
                   tuneGrid = gbmGrid2,
                   distribution = "bernoulli")
  c_boosted <- predict(gbmFit2,newdata=data.test)
  confusion2 <- confusionMatrix(c_boosted,data.test$c)
  accuracy.all.leading[i,2] <- confusion2$overall[1]
  rf_mtry <- train(c~.-y,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE)
  c_rf <- predict(rf_mtry$finalModel,data.test)
  confusion3 <- confusionMatrix(c_rf,data.test$c)
  accuracy.all.leading[i,3] <- confusion3$overall[1]
  bag <- randomForest(c~. -y, data = data.train,
                      mtry=5)
  predbag <- predict(bag,data.test)
  confusion4 <- confusionMatrix(predbag,data.test$c)
  accuracy.all.leading[i,4] <- confusion4$overall[1]
}
boxplot(accuracy.all.leading,,lwd=2,main="Accuracy All Leading")
plot(accuracy.all.leading[,2],
     accuracy.all.leading[,4],xlab = "Accuracy of Boosting",
     ylab = "Accuracy of Bagging", main = "Boosting v.s. RF")
abline(a=0, b=1, col='red')
mean(accuracy.all.leading[,3]==
       accuracy.all.leading[,4])
mean(accuracy.all.leading[,3]>
       accuracy.all.leading[,4])
mean(accuracy.all.leading[,3]<
       accuracy.all.leading[,4])




