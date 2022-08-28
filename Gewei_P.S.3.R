rm(list=ls())
### (a)
#install.packages("mvtnorm")
library(mvtnorm)
###Random draws from a multivariate normal distribution
library(MASS)
###required to fit lda and qda commands below
set.seed(527)
n1 <- 300
n2 <- 500
sigma <- matrix(c(16,-2,-2,9),c(2,2),byrow = T)
mu1 <- c(-3,3)
mu2 <- c(5,5)
K_1 <- rep("Class 1",n1)
K_2 <- rep("Class 2",n2)
K <- c(K_1,K_2)

data.gen <- function(mu1,mu2,sigma){
x_class1 <- mvtnorm::rmvnorm(n1,mean = mu1,sigma = sigma)
x_class2 <- rmvnorm(n2,mu2,sigma)
total <- rbind(x_class1,x_class2)
K <- as.factor(K)
data <- data.frame(K,total[,1],total[,2])
colnames(data) <- c("K","x1","x2")
return(data)
}

sim <- data.gen(mu1,mu2,sigma)
sim.test <- data.gen(mu1,mu2,sigma)

colors <- c("coral2","cornflowerblue")
poin <- c(15,7)
plot(sim[,3],sim[,2],col=colors[sim$K],
     xlab = "x2",ylab = "x1",main = "Class 1&2",
     pch=poin[sim$K])
legend(10,-12,legend=c("Class 1","Class 2"),
       fill = c("coral2","cornflowerblue"))

plot(sim.test[,3],sim.test[,2],col=colors[sim.test$K],
     xlab = "x2",ylab = "x1",main = "Test Class 1&2",
     pch=poin[sim.test$K])
legend(10,-10,legend=c("Class 1","Class 2"),
       fill = c("coral2","cornflowerblue"))


## (b)
lda.m <- lda(K ~ x1 + x2,data = sim)
plot(lda.m)
lda.m
lda.prd <- predict(lda.m,sim)
lda.class <- lda.prd$class
lda.confusion <- table(lda.class,sim[,1])
lda.confusion

## try to change threshold
T <- 5/8 # T is the threshold to be class 2
lda.prdclass <- rep("Class 1",length(lda.prd$class))
lda.prdclass[lda.prd$posterior[,2]>T] <- "Class 2"
lda.prdclass <- as.factor(lda.prdclass)
lda.confusion.T <- table(lda.prdclass,sim[,1])
######

### qda case 
qda <- qda(K ~ x1 + x2,data = sim)
qda
qda.prd <- predict(qda,sim)
qda.class <- qda.prd$class
qda.confusion <- table(qda.class,sim[,1])
qda.confusion


## logistic case 
logistic <- glm(K ~ x1 + x2,data = sim,
                family=binomial)
logistic
plot(logistic)
log.prd <- predict(logistic,sim,type = "response")
plot(log.prd)
contrasts(sim[,1])  ## see the coding rule in R
log.prdclass <- rep("Class 1",length(log.prd))
log.prdclass[log.prd>0.5] <- "Class 2"
log.prdclass <- as.factor(log.prdclass)
log.confusion <- table(log.prdclass,sim[,1])
log.confusion

########### (c)
lda.prd.test <- predict(lda.m,sim.test)
lda.class.test <- lda.prd.test$class
lda.confusion.test <- table(lda.class.test,sim.test[,1])
lda.confusion.test
# function of TE PE need a 2x2 confusion table 
TE <- function(confusion){
  (confusion[2,1]+confusion[1,2])/
    (confusion[1,1]+confusion[2,2]+confusion[2,1]+confusion[1,2])
}
PE<- function(confusion){
  (confusion[2,1]+confusion[1,2])/
    (confusion[1,1]+confusion[2,2]+confusion[2,1]+confusion[1,2])
}

lda.TE <- TE(lda.confusion)
lda.PE <- PE(lda.confusion.test)

log.prd.test <- predict(logistic,
    newdata=sim.test,
                    type = "response")
log.prdclass.test <- rep("Class 1",length(log.prd))
log.prdclass.test[log.prd.test>0.5] <- "Class 2"
log.prdclass.test <- as.factor(log.prdclass.test)
log.confusion.test <- table(log.prdclass.test,sim.test[,1])
log.confusion.test

log.TE <- TE(log.confusion)
log.PE <- PE(log.confusion.test)
lda.TE
lda.PE
log.TE
log.PE
## (d)
# sensitivity = true predicted class 2/total true class 2
# specificity = true predicted class 1/total true class 1
# under 0.5 threshold
sen <- function(confusion){
  confusion[2,2]/sum(confusion[,2])
}
spe <- function(confusion){
  confusion[1,1]/sum(confusion[,1])
}

sen.lda <- sen(lda.confusion)
sen.lda.test <- sen(lda.confusion.test)
sen.lda
sen.lda.test
sen.log <- sen(log.confusion)
sen.log.test <- sen(log.confusion.test)
sen.log
sen.log.test

spe.lda <- spe(lda.confusion)
spe.lda.test <- spe(lda.confusion.test)
spe.lda
spe.lda.test
spe.log <- spe(log.confusion)
spe.log.test <- spe(log.confusion.test)
spe.log
spe.log.test

TE(lda.confusion.T)
sen(lda.confusion.T)
spe(lda.confusion.T)

##############################################
##############################################
##############################################
##############################################
############### exercise 2 ###################
##############################################
##############################################
##############################################
##############################################

## (a)
set.seed(527)
n1 <- 300
n2 <- 500
sigma <- matrix(c(16,-2,-2,9),c(2,2),byrow = T)
mu1 <- c(-3,3)
mu2 <- c(5,5)
K_1 <- rep("Class 1",n1)
K_2 <- rep("Class 2",n2)
K <- c(K_1,K_2)
R <- 100
results <- matrix(NA,100,8)
for (i in 1:R) {
  sim <- data.gen(mu1,mu2,sigma)
  sim.test <- data.gen(mu1,mu2,sigma)
  lda.m <- lda(K ~ x1 + x2,data = sim)#lda
  lda.prd <- predict(lda.m,sim)
  lda.class <- lda.prd$class
  lda.confusion <- table(lda.class,sim[,1])
  logistic <- glm(K ~ x1 + x2,data = sim,
                  family=binomial)#logistic
  log.prd <- predict(logistic,sim,type = "response")
  log.prdclass <- rep("Class 1",length(log.prd))
  log.prdclass[log.prd>0.5] <- "Class 2"
  log.prdclass <- as.factor(log.prdclass)
  log.confusion <- table(log.prdclass,sim[,1])
  lda.prd.test <- predict(lda.m,sim.test)
  lda.class.test <- lda.prd.test$class
  lda.confusion.test <- table(lda.class.test,sim.test[,1])
  log.prd.test <- predict(logistic,
                          newdata=sim.test,
                          type = "response")
  log.prdclass.test <- rep("Class 1",length(log.prd))
  log.prdclass.test[log.prd.test>0.5] <- "Class 2"
  log.prdclass.test <- as.factor(log.prdclass.test)
  log.confusion.test <- table(log.prdclass.test,sim.test[,1])
  log.TE <- TE(log.confusion)
  log.PE <- PE(log.confusion.test)
  lda.TE <- TE(lda.confusion)
  lda.PE <- PE(lda.confusion.test)
  results[i,1] <- lda.TE
  results[i,2] <- lda.PE
  results[i,3] <- log.TE
  results[i,4] <- log.PE
  results[i,5] <- sen(lda.confusion.test)
  results[i,6] <- spe(lda.confusion.test)
  results[i,7] <- sen(log.confusion.test)
  results[i,8] <- spe(log.confusion.test)
}
colnames(results) <- c("lda.te","lda.pe",
                       "log.te","log.pe",
                       "lda.sen","lda.spe",
                       "log.sen","log.spe")
targets <- c(mean(results[,1]),mean(results[,2]),
             mean(results[,3]),mean(results[,4]),
             mean(results[,5]),mean(results[,6]),
             mean(results[,7]),mean(results[,8]))
targets
boxplot(results[,1],results[,3])
boxplot(results[,2],results[,4],main="Prediction Error")
legend("bottom",legend=c("Left P.E.LDA","Right P.E.Logistic"))
boxplot(results[,5],results[,6],
       results[,7],results[,8],main="Prediction Error")
legend("bottomleft",
       legend=c("Top Left Sensitivity.LDA",
                "Bottom Left Specificity.LDA",
                "Top Right Sensitivity.Logistic",
                "Bottom Right Specificity.Logistic"))

#### (b)
# the objectives of lda and logit are both classifications
# LDA needs specific distribution, and Logistic needs
# some "between" data, thus, the bigger separation of 
# two distribution, the higher variance of Log, but LDA dosen't suffer

set.seed(527)
n1 <- 300
n2 <- 500
sigma <- matrix(c(16,-2,-2,9),c(2,2),byrow = T)
mu1 <- c(-3,3)
mu2 <- c(5,5)
K_1 <- rep("Class 1",n1)
K_2 <- rep("Class 2",n2)
K <- c(K_1,K_2)
R <- 5
results <- matrix(NA,5,2)
d <- seq(-15,15,0.5)
re2 <- c()

for (j in 1:length(d)) {
  mu2 <- c(d[j],d[j])
 for (i in 1:R) {
  sim <- data.gen(mu1,mu2,sigma)
  sim.test <- data.gen(mu1,mu2,sigma)
  lda.m <- lda(K ~ x1 + x2,data = sim)#lda
  logistic <- glm(K ~ x1 + x2,data = sim,
                  family=binomial)#logistic
  lda.prd.test <- predict(lda.m,sim.test)
  lda.class.test <- lda.prd.test$class
  lda.confusion.test <- table(lda.class.test,sim.test[,1])
  log.prd.test <- predict(logistic,
                          newdata=sim.test,
                          type = "response")
  log.prdclass.test <- rep("Class 1",length(log.prd))
  log.prdclass.test[log.prd.test>0.5] <- "Class 2"
  log.prdclass.test <- as.factor(log.prdclass.test)
  log.confusion.test <- table(log.prdclass.test,sim.test[,1])
  log.PE <- PE(log.confusion.test)
  lda.PE <- PE(lda.confusion.test)
  results[i,1] <- lda.PE
  results[i,2] <- log.PE
}
  targets <- c(mean(results[,1]),mean(results[,2]))
re2[j] <- (targets[1]-targets[2])/targets[1]
}
plot(d,re2,xlab = "Distance",
     ylab = "Relative difference of P.E.",
     main = "P.E Deviation of LDA and Logistic",
     pch=16)
lines(d[-c(1,2)],predict(lm(re2[-c(1,2)]~d[-c(1,2)]+I(d[-c(1,2)]^2))),
      col="red",lwd=2,lty=2)
mtext("(P.E.LDA - P.E.Logistic)/P.E.LDA", side = 3)
lines(d,rep(0,length(d)),col="blue",lwd=2)
summary(lm(re2~d+I(d^2)))
plot(d[21:41],re2[21:41],xlab = "Distance",
     ylab = "Relative difference of P.E.",
     main = "P.E Deviation of LDA and Logistic")
lines(d[20:41],rep(0,length(d[20:41])),col="red",
      lwd=2)

##############################################
# this also means, that the nearer two distributions
# the higher prob. that LDA plays worse than Log
results <- matrix(NA,R,2)
d <- seq(-15,15,0.1)
re2 <- c()

for (j in 1:length(d)) {
  mu2 <- c(d[j],d[j])
  for (i in 1:R) {
    sim <- data.gen(mu1,mu2,sigma)
    sim.test <- data.gen(mu1,mu2,sigma)
    lda.m <- lda(K ~ x1 + x2,data = sim)#lda
    lda.prd.test <- predict(lda.m,sim.test)
    lda.class.test <- lda.prd.test$class
    lda.confusion.test <- table(lda.class.test,sim.test[,1])
    logistic <- glm(K ~ x1 + x2,data = sim,
                    family=binomial)#logistic
    log.prd.test <- predict(logistic,
                            newdata=sim.test,
                            type = "response")
    log.prdclass.test <- rep("Class 1",length(log.prd))
    log.prdclass.test[log.prd.test>0.5] <- "Class 2"
    log.prdclass.test <- as.factor(log.prdclass.test)
    log.confusion.test <- table(log.prdclass.test,sim.test[,1])
    log.PE <- PE(log.confusion.test)
    lda.PE <- PE(lda.confusion.test)
    results[i,1] <- lda.PE
    results[i,2] <- log.PE
  }
  targets <- c(mean(results[,1]),mean(results[,2]))
  re2[j] <- (targets[1]-targets[2])/targets[1]
}
re20 <- na.omit(re2)
re20 <- re20[re20>-100]# remove the -inf 
mean(re20>=0)
plot(d,re2,xlab = "Distance",
     ylab = "Relative difference of P.E.",
     main = "P.E Deviation of LDA and Logistic")
#lines(d,predict(lm(re2~d+I(d^2))),col="red")
mtext("(P.E.LDA - P.E.Logistic)/P.E.LDA", side = 3)
lines(d,rep(0,length(d)),col="blue",lwd=2)
bz <- re2[re2>=0]
points(d[re2>=0],bz,col="red")

plot(d[70:230],re2[70:230],xlab = "Distance",
     ylab = "Relative difference of P.E.",
     main = "P.E Deviation of LDA and Logistic")
mtext("(P.E.LDA - P.E.Logistic)/P.E.LDA", side = 3)
lines(d,rep(0,length(d)),col="blue",lwd=2)
bz <- re2[re2>=0]
points(d[re2>=0],bz,col="red")

#colors <- c("coral2","cornflowerblue")
#poin <- c(15,7)
#plot(sim[,3],sim[,2],col=colors[sim$K],
#     xlab = "x2",ylab = "x1",main = "Class 1&2",
#     pch=poin[sim$K])
#legend(x="topleft",legend=c("Class 1","Class 2"),
#       fill = c("coral2","cornflowerblue"))

# Also the textbook says that for smaller n,
# LDA works better than logistic
############################################
# If one group¡¯s n is much less than the others
# then estimated logistic function may ¡°shift¡± to 
# the bigger sample class more, 
# which leads to higher variance


### (c)
set.seed(527)
n1 <- 300
n2 <- 500
sigma <- matrix(c(16,-2,-2,9),c(2,2),byrow = T)
mu1 <- c(-3,3)
mu2 <- c(5,5)
K_1 <- rep("Class 1",n1)
K_2 <- rep("Class 2",n2)
K <- c(K_1,K_2)

sim_c <- data.gen(mu1,mu2,sigma)
sim.test_c <- data.gen(mu1,mu2,sigma)

Threshold <- seq(0.05,0.95,0.05)
resultcurve <- matrix(NA,length(Threshold),6)
colnames(resultcurve) <- c("TE lda","1-SE lda",
  "1-SP lda","TE log","1-SE log","1-SP log")
lda.m <- lda(K ~ x1 + x2,data = sim_c)
logistic <- glm(K ~ x1 + x2,data = sim_c,
                family=binomial)
lda.prd <- predict(lda.m,sim.test_c)
log.prd <- predict(logistic,sim.test_c,type = "response")

for(i in 1:length(Threshold)) {
T <- Threshold[i]
lda.prdclass <- rep("Class 1",length(lda.prd$class))
lda.prdclass[lda.prd$posterior[,2]>T] <- "Class 2"
lda.prdclass <- as.factor(lda.prdclass)
lda.confusion.T <- table(lda.prdclass,sim.test_c[,1])

log.prdclass <- rep("Class 1",length(log.prd))
log.prdclass[log.prd>T] <- "Class 2"
log.prdclass <- as.factor(log.prdclass)
log.confusion <- table(log.prdclass,sim.test_c[,1])

  resultcurve[i,1] <- TE(lda.confusion.T)
  resultcurve[i,2] <- 1-sen(lda.confusion.T)
  resultcurve[i,3] <- 1-spe(lda.confusion.T)
  resultcurve[i,4] <- TE(log.confusion)
  resultcurve[i,5] <- 1-sen(log.confusion)
  resultcurve[i,6] <- 1-spe(log.confusion)
}

###### lda
plot(Threshold,resultcurve[,1],"l",col="black",
     ylab = "Error Rate",ylim = c(0,0.8),
     main = "LDA method")
lines(Threshold,resultcurve[,2],col="red",lty=2)
lines(Threshold,resultcurve[,3],col="blue",lty=3)
points(Threshold[match(min(resultcurve[,1]),
                       resultcurve[,1])],
       min(resultcurve[,1]),pch=18,col="green")
legend("topright",legend = c("Total Error Rate",
                  "Error Rate for True Class 2",
                "Error Rate for True Class 1"),
       lty = 1:3,col = c("black","red","blue"))

### log
plot(Threshold,resultcurve[,4],"l",col="black",
     ylab = "Error Rate",ylim = c(0,0.8),
     main = "Logisitc Model")
lines(Threshold,resultcurve[,5],col="red",lty=2)
lines(Threshold,resultcurve[,6],col="blue",lty=3)
points(Threshold[match(min(resultcurve[,4]),
                       resultcurve[,4])],
       min(resultcurve[,4]),pch=18,col="green")
legend("topright",legend = c("Total Error Rate",
                             "Error Rate for True Class 2",
                             "Error Rate for True Class 1"),
       lty = 1:3,col = c("black","red","blue"))

###################################################
### total 
plot(Threshold,resultcurve[,1],"l",col="black",
     ylab = "Error Rate",ylim = c(0,0.8),
     main = "TWO method")
grid()
lines(Threshold,resultcurve[,2],col="red",lty=2,lwd=1)
lines(Threshold,resultcurve[,3],col="blue",lty=3,lwd=1)
lines(Threshold,resultcurve[,4],col="burlywood1",lty=4,lwd=4)
lines(Threshold,resultcurve[,5],col="darkolivegreen",lty=5,lwd=2)
lines(Threshold,resultcurve[,6],col="darkorchid4",lty=6,lwd=2)
points(Threshold[match(min(resultcurve[,1]),
                       resultcurve[,1])],
       min(resultcurve[,1]),pch=18,col="brown1",cex=2)
points(Threshold[match(min(resultcurve[,4]),
                       resultcurve[,4])],
       min(resultcurve[,4]),pch=20,col="darkslategrey",cex=2)
legend("topright",legend = c("Total Error Rate (LDA)",
"Error Rate for Class 2 (LDA)","Error Rate for Class 1 (LDA)",
"Total Error Rate (LOG)",
"Error Rate for Class 2 (LOG)","Error Rate for Class 1 (LOG)",
"Lowest Test error of LDA","Lowest Test error of Logistic"),
       pch=c(NA,NA,NA,NA,NA,NA,18,20),
lwd=c(1,1,1,4,2,2,NA,NA),lty = c(1:6,NA,NA),col = c("black","red","blue","burlywood1",
                         "darkolivegreen","darkorchid4","brown1","darkslategrey"))






