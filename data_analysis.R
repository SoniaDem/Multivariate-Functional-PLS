##
## ========== DATA PREPARATION ==========
##

dd <- read.csv("~/MastersStuff/dataset.csv")
#normalising the data so one column addds up to one
d<-data.frame("ppm"=dd[,1],dd[,2:141]/sum(dd[,2:141]))
#ppm from 0 to 8
d<-subset(d,d$ppm>=0&d$ppm<=8)

#set bin size=0.16 ppm
bs=0.16
N = 50
for (i in 1:N){
  nam <- paste("X", i, sep = "")
  if (i<N)
  {assign(nam, colMeans(subset(d[,2:141],d$ppm >= min(d$ppm)+(i-1)*bs 
                               & d$ppm < min(d$ppm)+i*bs)))}
  else
  {assign(nam, colMeans(subset(d[,2:141],d$ppm>=min(d$ppm)+(i-1)*bs&d$ppm<max(d$ppm))))}
}

s=as.data.frame(matrix(0,nrow=140,ncol=N))
for (i in 1:N){
  s[,i]=as.data.frame(get(paste('X',i,sep="")))
}

s<-log(1000000*s+1)

# add a row of column names to the new data frame s to easily identify the rows
samples <- names(d)[-1]
s <- cbind(samples,s)

par(mfrow=c(1,1))
limit<-c(min(s[1,]),max(s[1,]))
plot(d[1,], type = "l", col = "black", xlab="ppm",ylab="intensity",main="subject 1")


##
## ========== GET THE CORRECT DATA TIME POINTS ==========
##

#subset the data to only take people before operation
numbers <- c(seq(1,18),21)
snam <- rep(NA, 19)

for (i in 1:19){
  nam <- paste("k", numbers[i], ".pre", sep = "")
  snam[i] <- nam
}

vars <- s$samples %in% snam
s.pre <- s[vars,]

#remove patient nr 15

s.pre <- subset(s.pre, s.pre$samples != 'k15.pre')
kt.pre <- s.pre

#subset the data to take people during and immediately after operation

numbers <- c(seq(1,18),21)[-15]
point <- c('.int', '.post', '.02')
s.vars <- rep(NA,3)

for(timep in point){
  snam <- rep(NA, 18)
  for (i in 1:18){
    nam <- paste("k", numbers[i], timep, sep = "")
    snam[i] <- nam
  }
  print(snam)
  vars <- s$samples %in% snam
  a <- s[vars,]
  assign(paste("kt", timep, sep = ""), a)
}

##
## ========== PREPARE DATA FOR FPLS ==========
##

library(fda)
library(RColorBrewer)

#set variables
bin <- seq(0,8, by=8/49)
nobs = length(bin)
knots <- seq(0,8, by = 8/49)
nknots = length(knots) - 1
norder = 4
nbasis = length(knots) + norder - 2

#create a bspline basis of 52 vectors in the range of 0 to 8
basisobj <- create.bspline.basis(c(min(bin), max(bin)), nbasis, norder, knots)
#plot the bsplie
# plot(basisobj, main = 'Basis Functions')

##
## ========== DEFINING THE SMOOTHED DATA ==========
##

lambda= 10^(-3)
fdParobj = fdPar(fdobj=basisobj, lambda=lambda)

xfdata.pre <- smooth.basis(bin, t(as.matrix(kt.pre[,-1])), fdParobj)$fd
xfdata.int <- smooth.basis(bin, t(as.matrix(kt.int[,-1])), fdParobj)$fd
xfdata.post <- smooth.basis(bin, t(as.matrix(kt.post[,-1])), fdParobj)$fd
xfdata.02 <- smooth.basis(bin, t(as.matrix(kt.02[,-1])), fdParobj)$fd

#plot the 4 time points in a 2,2 grid

par(mfrow=c(2,3))
# cat("\n")
plot(xfdata.pre, type="b", xlab='ppm')
title("Smoothed Pre-Transplant Data")
plot(xfdata.int, type="b", xlab='ppm')
title("Smoothed Intermediate Data")
plot(xfdata.post, type="b", xlab='ppm')
title("Smoothed Post-Transplant Data")
# plot(xfdata.02, type="b")
# title("2 Days Post-Transplant Data")

##
## ========== DEFINING THE RESPONSE VARIABLE ==========
##

# defining the function that will give a binary response variable given an outcome group
# fun feature: input must be a string

response <- function(ogroup) {
  c.red <- c(4,7,9,17,18)
  c.blue <- c(8,12,14,21)
  c.black <- c(1,2,3,5,6,10,11,13,16)
  
  if(ogroup == 'A'){
    success <-  c(c.black, c.red)
    acute.rejection = c.blue
  }
  
  if(ogroup == 'B'){
    success <-  c(c.black, c.blue)
    acute.rejection = c.red
  }
  
  if(ogroup == 'C'){
    success <-  c(c.black)
    acute.rejection = c(c.blue, c.red) 
  }
  
  if(ogroup != 'A' & ogroup != 'B' & ogroup != 'C'){
    print('error: please input one of the following outcome groups: A, B, C')
  }
  
  response <- rep(0, 18)
  
  a <- rep(NA, length(success))
  n=1
  
  for (i in success){
    a[n]<- paste("k", i, ".pre", sep = "")
    n=n+1
  }  
  
  j = which(kt.pre$samples %in% a )
  response[j] = 1
  
  return(as.numeric(response))
}

##
## ========== PREPARATION FOR FPLS ==========
##

library(ftsa)
library(fda.usc)
library(caret)
library(ROCR)
library(pROC)

#tranform an fd object to an fdata object
fda.kt.pre <- fdata(xfdata.pre)
fda.kt.int <- fdata(xfdata.int)
fda.kt.post <- fdata(xfdata.post)
fda.kt.02 <- fdata(xfdata.02)

#classification function: given the fpls model will output the classification

class <- function(pls_mod){
  fitted <- pls_mod$fitted.values
  
  alpha <- - inprod.fdata(pls_mod$beta.est, mean(pls_mod$fdataobj))
  
  len <- length(pls_mod$fdataobj)
  Phi <- fitted + rep(alpha, len)
  
  Phi_O <- rep(NA, len)
  for (i in 1:len){
    if (Phi[i] < mean(Phi)-var(Phi)){Phi_O[i]=0}
    if (Phi[i] >= mean(Phi)-var(Phi)){Phi_O[i]=1}
  }
  
  return(Phi_O)
}

##
## ========== FUNCTIONS TO MEASURE MODELS ==========
##

auc_score <- function(pred, truth, ...) {
  predob = prediction(pred, truth)
  area <- auc(truth, pred)
  return(area)
}

#function to plot the ROC plot
rocplot <- function(pred, truth, ...) {
  predob = prediction(pred, truth)
  perf = performance(predob, "tpr", "fpr")
  plot(perf, ...)
  area <- auc(truth, pred)
  area <- format(round(area, 4), nsmall = 4)
  text(x=0.8, y=0.1, labels = paste("AUC =", area))
  
  # the reference x=y line
  segments(x0=0, y0=0, x1=1, y1=1, col="gray", lty=2)
}

#function that outputs accuracy given a confusion matrix
acc <- function(table){
  return(( (table[1,1] + table[2,2]) / sum(table)))
}

##
## ========== FPLS ON DIFFERENT TIME POINTS ==========
##

# set the response variable
y <- response('C')
yy <- as.factor(y)

#removing some patients for two time points as they were not present
y.post <- y[-c(13,16)]
yy.post <- as.factor(y.post)
y.02 <- y[-c(1,16)]
yy.02 <- as.factor(y.02)

# PRE

pls_model_pre <- fregre.pls(fda.kt.pre, y , l = 4, lambda = 0)
res.pre <- round(pls_model_pre$fitted.values)
tab.pre <- confusionMatrix(as.factor(res.pre), yy)$table
acc.pre <- acc(tab.pre)
rocplot(as.numeric(res.pre), y)

# INT

pls_model_int <- fregre.pls(fda.kt.int, y , l = 4, lambda = 0)
res.int <- round(pls_model_int$fitted.values)
tab.int <- confusionMatrix(as.factor(res.int), yy)$table
acc.int <- acc(tab.int)
rocplot(as.numeric(res.int), y)

# POST

pls_model_post <- fregre.pls(fda.kt.post, y.post , l = 4, lambda = 0)
res.post <- round(pls_model_post$fitted.values)
tab.post <- confusionMatrix(as.factor(res.post), yy.post)$table
acc.post <- acc(tab.post)
rocplot(as.numeric(res.post), y.post)

# 2 DAYS

pls_model_02 <- fregre.pls(fda.kt.02, y.02 , l = 4, lambda = 0)
res.02 <- class(pls_model_02)
tab.02 <- confusionMatrix(as.factor(res.02), yy.02)$table
acc.02 <- acc(tab.02)
rocplot(as.numeric(res.02), y.02)

print(acc.pre, acc.int, acc.post, acc.02)
##
## ========== FPLS - GENERATING ACCURACY PLOTS ==========
##

for (a in c('A', 'B', 'C')){
  
  y <- response(a)
  print(y)
  yy <- as.factor(y)
  accuracy.pre <- rep(NA,4)
  auc.pre <- rep(NA,4)
  
  for (i in 2:5){
    pls_model <- fregre.pls(fda.kt.pre, y , l = i, lambda = 0)
    resp<- class(pls_model)
    table <- confusionMatrix(as.factor(resp), yy)$table
    accuracy.pre[i-1] <- acc(table)
    auc.pre[i-1] <- auc_score(as.numeric(resp), yy)
  }
  
  assign(paste('accuracy.pre', a, sep = ""), accuracy.pre)
  assign(paste('auc.pre', a, sep = ""), auc.pre)  
}

par(mfrow=c(1,3))
plot(accuracy.preA, main='Outcome A Accuracy', xlab='Number of PLS Components', ylab = 'Accuracy', type='l', col='royalblue3', ylim=c(0.6,1))
plot(accuracy.preB, main='Outcome B Accuracy', xlab='Number of PLS Components', ylab = 'Accuracy', type='l', col='royalblue3', ylim=c(0.6,1))
plot(accuracy.preC, main='Outcome C Accuracy', xlab='Number of PLS Components', ylab = 'Accuracy', type='l', col='royalblue3', ylim=c(0.6,1))

##
## ========== MULTIVARIATE MODEL: ATTEMPT USING PC COMPONENTS ==========
##
y <- response('C')
yy <- as.factor(y)

fdatapls.pre <- fdata2pls(fda.kt.pre, y, ncomp = 2, lambda = 0)
T.pre <- fdatapls.pre$x

fdatapls.int <- fdata2pls(fda.kt.int, y, ncomp = 2, lambda = 0)
T.int <- fdatapls.int$x

#linking two X matrices together


lm.mod <- glm(T.int[,1] ~ T.pre)
lm.mod2 <- glm(T.int[,2] ~ T.pre)
H_t <-  cbind(lm.mod$residuals, lm.mod2$residuals)

tildeX <- as.matrix(cbind(T.pre, H_t))
fin.mod.residuals <- glm(y ~ tildeX, family='binomial', control = list(maxit = 50))

outcome.fin.res <- as.factor(round(fin.mod.residuals$fitted.values))
fin.res.tab <- confusionMatrix(outcome.fin.res, yy)$table
acc.mul1 <- acc(fin.res.tab) 
rocplot(as.numeric(outcome.fin.res), y)

##
## ========== MULTIVARIATE MODEL WITH 3 TIME POINTS==========
##

#kt.post has two people missing 13 & 17

fdatapls1.pre <- fdata2pls(fda.kt.pre[-c(13,16)], as.numeric(y.post), ncomp = 2, lambda = 0)
T.pre1 <- fdatapls1.pre$x

fdatapls1.int <- fdata2pls(fda.kt.int[-c(13,16)], as.numeric(y.post), ncomp = 1, lambda = 0)
T.int <- fdatapls1.int$x

fdatapls.post <- fdata2pls(fda.kt.post, as.numeric(y.post), ncomp = 1, lambda = 0)
T.post <- fdatapls.post$x


#linking two X matrices together

lm.mod <- glm(T.int[,1] ~ T.pre1)
#lm.mod2 <- glm(T.int[,2] ~ T.pre)
H_t <-  cbind(lm.mod$residuals) #, lm.mod2$residuals)

lm.mod.post <- lm(T.post ~ T.pre1)

H_t2 <- lm.mod.post$residuals

postbyint <- lm(H_t2 ~ H_t)
Htt <- postbyint$residuals

#T.post <- append(T.post, NA, after=12)
#T.post <- append(T.post, NA, after=15)
#lm.mod.post <- lm(T.post ~ T.pre1)

#H_t2 <- lm.mod.post$residuals
#H_t2 <- append(H_t2 , NA, after=12)
#H_t2<- append(H_t2 , NA, after=15)

#postbyint <- lm(H_t2 ~ H_t)
#Htt <- postbyint$residuals
#Htt <- append(Htt , mean(Htt), after=12)
#Htt<- append(Htt , mean(Htt), after=15)

tildeX3 <- cbind(T.pre1, H_t, Htt)
fin.mod3.residuals <- glm(y.post ~ tildeX3, family='binomial', control = list(maxit = 50), na.action = na.pass)

outcome.fin.res3 <- as.factor(round(fin.mod3.residuals$fitted.values))
fin.res.tab3 <- confusionMatrix(outcome.fin.res3, as.factor(y.post))$table
rocplot(as.numeric(outcome.fin.res3), as.numeric(y.post))
acc.mul3 = acc(fin.res.tab3)

##
## ========== FPCA ==========
##

# set the response variable
y <- response('A')
yy <- as.factor(y)

#removing some patients for two time points as they were not present
y.post <- y[-c(13,16)]
yy.post <- as.factor(y.post)
y.02 <- y[-c(1,16)]
yy.02 <- as.factor(y.02)

#the follwing atasets can be considered input into FPCA: xfdata.pre, xfdata.int, xfdata.post, xfdata.02
mod.pca5 <- pca.fd(xfdata.pre, nharm = 5, harmfdPar=fdPar(xfdata.pre), centerfns = TRUE)
mod_scores <- mod.pca5$scores
sc.data <- as.data.frame(mod_scores)
sc.data <-  cbind(sc.data, response)
sc.data$response <- as.factor(sc.data$response)

#logistic regression
mod_l <- glm(response ~ .  , data = sc.data, family='binomial')
preds.l = predict(mod_l, newdata=sc.data, type = "response")
preds.l <- as.factor(round(preds.l))
cm.glm = confusionMatrix(preds.l, sc.data[,length(sc.data)])
print(cm.glm$table)
rocplot(as.numeric(preds.l), as.numeric(y)) 

#random forest
mod_rf <- randomForest(response ~ ., data = sc.data)
cm.rf <- mod_rf$confusion
rocplot(as.numeric(mod_rf$predicted), as.numeric(y)) 

