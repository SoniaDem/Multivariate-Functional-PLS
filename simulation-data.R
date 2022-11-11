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
library(fda.usc)

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

lambda= 10^(-20)
fdParobj = fdPar(fdobj=basisobj, lambda=lambda)

testing = kt.pre[,-1]

xfdata.pre <- smooth.basis(bin, t(as.matrix(testing)), fdParobj)$fd

xfdata.int <- smooth.basis(bin, t(as.matrix(kt.int[,-1])), fdParobj)$fd

xfdata.post <- smooth.basis(bin, t(as.matrix(kt.post[,-1])), fdParobj)$fd


##
## ========== SIMULATION DATA: PCA COMPONENTS ==========
##

# 1st TIME POIINT

# create a PCA decomposition usign 3 components
pca_pre <- pca.fd(xfdata.pre, nharm = 5, harmfdPar=fdPar(xfdata.pre),
                  centerfns = TRUE)

# isolate the harmonics and change their class to match the PLS approach and 
# make it easier to multiply with the basis matrices 
hr_pre <- pca_pre$harmonics
new_hr_pre <- fdata(hr_pre)
sc_pre <- pca_pre$scores


# 2nd TIME POINT

# create a PCA decomposition usign 3 components
pca_int <- pca.fd(xfdata.int, nharm = 5, harmfdPar=fdPar(xfdata.int),
                  centerfns = TRUE)

# isolate the harmonics and change their class to match the PLS approach and 
# make it easier to multiply with the basis matrices 
hr_int <- pca_int$harmonics
new_hr_int <- fdata(hr_int)
sc_int <- pca_int$scores

# 3rd TIME POINT

# create a PCA decomposition usign 3 components
pca_post<- pca.fd(xfdata.post, nharm = 5, harmfdPar=fdPar(xfdata.post),
                  centerfns = TRUE)

# isolate the harmonics and change their class to match the PLS approach and 
# make it easier to multiply with the basis matrices 
hr_post <- pca_post$harmonics
new_hr_post <- fdata(hr_post)
sc_post <- pca_post$scores

##
## ========== CALCULATE THE DISTRIBUTION OF THE PCA SCORES ======= 
##

c.black <- c(1,2,3,5,6,10,11,13,15)
c.black.post <- c(1,2,3,5,6,10,11,15)
c.red <- c(4,7,8,9,12,14,16,17,18)
c.red.post <- c(4,7,8,9,12,13,15,16)

#1st TIME POINT
mean_pre_0 <- rep(NA, 5)
mean_pre_1 <- rep(NA, 5)
var_pre_0 <- rep(NA, 5)
var_pre_1 <- rep(NA, 5)

for ( i in 1:5){
  mean_pre_0[i] <- mean(sc_pre[c.red,i])
  mean_pre_1[i] <- mean(sc_pre[c.black,i]) 
  var_pre_0[i] <- var(sc_pre[c.red,i])
  var_pre_1[i] <- var(sc_pre[c.black,i]) 
}

#2nd TIME POINT

mean_int_0 <- rep(NA, 5)
mean_int_1 <- rep(NA, 5)
var_int_0 <- rep(NA, 5)
var_int_1 <- rep(NA, 5)

for ( i in 1:5){
  mean_int_0[i] <- mean(sc_int[c.red,i])
  mean_int_1[i] <- mean(sc_int[c.black,i])  
  var_int_0[i] <- var(sc_int[c.red,i])
  var_int_1[i] <- var(sc_int[c.black,i]) 
}

#3rd TIME POINT

mean_post_0 <- rep(NA, 5)
mean_post_1 <- rep(NA, 5)
var_post_0 <- rep(NA, 5)
var_post_1 <- rep(NA, 5)

for ( i in 1:5){
  mean_post_0[i] <- mean(sc_post[c.red.post,i])
  mean_post_1[i] <- mean(sc_post[c.black,i])  
  var_post_0[i] <- var(sc_post[c.red.post,i])
  var_post_1[i] <- var(sc_post[c.black,i]) 
}


#Adjust the variables a bit
mean_pre_0 <- mean_pre_0*10
mean_pre_1 <- mean_pre_1*10
var_pre_0 <- var_pre_0*5
var_pre_1 <- var_pre_1*5

mean_int_0 <- mean_int_0*10
mean_int_1 <- mean_int_1*10
var_int_0 <- var_int_0*5
var_int_1 <- var_int_1*5

mean_post_0 <- mean_post_0*10
mean_post_1 <- mean_post_1*10
var_post_0 <- var_post_0*5
var_post_1 <- var_post_1*5


##
## =========== MULTIVARIATE NORMAL DISTRIBUTION ============
##
# were missing patient 13 and 16 in the third time point y.post <- y[-c(13,16)]

#for group 0

for (i in 1:5){
  var11_0 <- var_pre_0[i] 
  var22_0 <- var_int_0[i]
  var33_0 <- var_post_0[i]
  cov12_0 <- cov(sc_pre[c.black,i], sc_int[c.black,i])
  cov13_0 <- cov(sc_pre[c.black,i][-8], sc_post[c.black.post,i])
  cov23_0 <- cov(sc_int[c.black,i][-8], sc_post[c.black.post,i])
  
  n0 <- 50
  mmu_0 <- matrix(c(mean_pre_0[i], mean_int_0[i], mean_post_0[i]), ncol=1)
  mcov_0 <- matrix(c(var11_0, cov12_0, cov13_0, 
                     cov12_0, var22_0, cov23_0,
                     cov13_0, cov23_0, var33_0), ncol=3)
  bivar <- mvrnorm(n = n0, mu = mmu_0, Sigma = mcov_0)  # Random sample from bivariate normal distribution
  name <- paste('zeroscores',i, sep = "", collapse = NULL)
  assign(name, bivar)
}

#for group 1

for (i in 1:5){
  var11_1 <- var_pre_1[i] 
  var22_1 <- var_int_1[i]
  var33_1 <- var_post_1[i]
  cov12_1 <- cov(sc_pre[c.red,i], sc_int[c.red,i])
  cov13_1 <- cov(sc_pre[c.red,i][-7], sc_post[c.red.post,i])
  cov23_1 <- cov(sc_int[c.red,i][-7], sc_post[c.red.post,i])
  
  n1 <- 50
  mmu_1 <- matrix(c(mean_pre_1[i], mean_int_1[i], mean_post_1[i]), ncol=1)
  mcov_1 <- matrix(c(var11_1, cov12_1, cov13_1, 
                     cov12_1, var22_1, cov23_1,
                     cov13_1, cov23_1, var33_1), ncol=3)
  bivar <- mvrnorm(n = n1, mu = mmu_1, Sigma = mcov_1)  # Random sample from bivariate normal distribution
  name <- paste('onescores',i, sep = "", collapse = NULL)
  assign(name, bivar)
}

## ========= TURN THE VECTORS INTO MATRICES FOR LATEER ==============
library(wakefield)

##attempt to make a mixed vector of scores 1


random_vec <- r_sample_logical(100, prob = NULL, name = "Logical")
true_values <- which(random_vec == TRUE)

mixedscores1 <- zeroscores1
mixedscores1[true_values] = onescores1[true_values]

averagescore = (2*zeroscores1+onescores1)/3

scorelist_zero <- cbind(zeroscores1[,1]*0.05, zeroscores2[,1], zeroscores3[,1], zeroscores4[,1], zeroscores5[,1],
                        zeroscores1[,2]*0.05, zeroscores2[,2], zeroscores3[,2], zeroscores4[,2], zeroscores5[,2],
                        zeroscores1[,3]*0.05, zeroscores2[,3], zeroscores3[,3], zeroscores4[,3], zeroscores5[,3])

scorelist_one <- cbind(zeroscores1[,1]*0.05, zeroscores2[,1], zeroscores3[,1],  zeroscores4[,1],  zeroscores5[,1],
                       zeroscores1[,2]*0.05, zeroscores2[,2], zeroscores3[,2],  zeroscores4[,2],  zeroscores5[,2],
                       zeroscores1[,3]*0.05, zeroscores2[,3], zeroscores3[,3],  zeroscores4[,3],  zeroscores5[,3])

for (i in 1:5){
  T0_pre <- replicate(52, scorelist_zero[,i])
  T0_int <- replicate(52, scorelist_zero[,i+5])
  T0_post <- replicate(52, scorelist_zero[,i+10])
  T1_pre <-  replicate(52, scorelist_one[,i])
  T1_int <-  replicate(52, scorelist_one[,i+5])
  T1_post <-  replicate(52, scorelist_one[,i+10])
  
  name_pre <- paste('Tpre',i, sep = "", collapse = NULL)
  name_int <- paste('Tint',i, sep = "", collapse = NULL)
  name_post <- paste('Tpost',i, sep = "", collapse = NULL)
  
  assign(name_pre, rbind(T0_pre,T1_pre))
  assign(name_int, rbind(T0_int,T1_int))
  assign(name_post, rbind(T0_post,T1_post))
}  


## ============ USE THE MATRICES TO GENERATE DATA ========

fhr_pre <- fdata(hr_pre)
fhr_int <- fdata(hr_int)
fhr_post <- fdata(hr_post)

Tlist_pre <- cbind(Tpre1, Tpre2, Tpre3, Tpre4, Tpre5)
Tlist_int <- cbind(Tint1, Tint2, Tint3, Tint4, Tint5)
Tlist_post <- cbind(Tpost1, Tpost2, Tpost3, Tpost4, Tpost5)

#add noise

noisepre <- rnorm(ncol(Tlist_pre)*nrow(Tlist_pre))
noiseint <- rnorm(ncol(Tlist_int)*nrow(Tlist_int))
noisepost <- rnorm(ncol(Tlist_post)*nrow(Tlist_post))

scale <- function(noise){(noise-min(noise))/(max(noise)-min(noise))*(1.2e-5-2e-6)+2e-6}
w.noisepre <- scale(noisepre)
w.noiseint <- scale(noiseint)
w.noisepost <- scale(noisepost)

npre_mat <- matrix(NA, nrow=100, ncol=260)
for (i in 1:100){
  for (k in 1:260)
    npre_mat[i,k] = w.noisepre[(i-1)*260+k]
}

nint_mat <- matrix(NA, nrow=100, ncol=260)
for (i in 1:100){
  for (k in 1:260)
    nint_mat[i,k] = w.noiseint[(i-1)*260+k]
}

npost_mat <- matrix(NA, nrow=100, ncol=260)
for (i in 1:100){
  for (k in 1:260)
    npost_mat[i,k] = w.noisepost[(i-1)*260+k]
}

Tlist_pre = Tlist_pre + npre_mat
Tlist_int = Tlist_int + nint_mat
Tlist_post = Tlist_post + nint_mat

pre_sim <- (fhr_pre$data[1,])*t(Tpre1)
int_sim <- (fhr_int$data[1,])*t(Tint1)
post_sim <- (fhr_post$data[1,])*t(Tpost1)

for (i in 2:5){
  pre_sim = pre_sim + (fhr_pre$data[i,])*t(Tlist_pre[,i:(i+51)])
  int_sim = int_sim + (fhr_int$data[i,])*t(Tlist_int[,i:(i+51)])
  post_sim = post_sim + (fhr_post$data[i,])*t(Tlist_post[,i:(i+51)])
}

data_mean_int <- c(0,colMeans(kt.int[,-1]),0)
data_mean_pre <- c(0, colMeans(kt.pre[,-1]),0)
data_mean_post <- c(0, colMeans(kt.post[,-1]),0)
#
pre_sim = pre_sim + data_mean_pre
int_sim = int_sim + data_mean_int
post_sim = post_sim + data_mean_post


## ============ CREATE SPLINES FOR THE DATA ========

#set variables
bin <- seq(0,8, by=8/51)
nobs = length(bin)
knots <- seq(0,8, by = 8/51)
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

lambda= 10^(-15)
fdParobj = fdPar(fdobj=basisobj, lambda=lambda)

sim_xfdata.pre <- smooth.basis(bin, pre_sim, fdParobj)$fd
sim_xfdata.int <- smooth.basis(bin, int_sim, fdParobj)$fd
sim_xfdata.post <- smooth.basis(bin, post_sim, fdParobj)$fd

#comparing real to simulated
par(mfrow=c(2,3))
plot(sim_xfdata.pre, type="b")
title("Simulated Pre-Transplant Data")
plot(sim_xfdata.int, type="b")
title("Simulated Intermediate Data")
plot(sim_xfdata.post, type="b")
title("Simulated Post-Transplant Data")
plot(xfdata.pre, type="b")
title("Real Pre-Transplant Data")
plot(xfdata.int, type="b")
title("Real Intermediate Data")
plot(xfdata.post, type="b")
title("Real Post-Transplant Data")

par(mfrow=c(1,3))
plot((xfdata.pre), type="b", xlab='ppm')
title("Smoothed Pre-Transplant Data")
plot(xfdata.int, type="b", xlab='ppm')
title("Smoothed Intermediate Data")
plot(xfdata.post, type="b", xlab='ppm')
title("Smoothed Post_Transplant Data")

#comparing real to simulated
par(mfrow=c(2,3))
plot(sim_xfdata.pre[1:50], type="b")
title("Simulated Pre-Transplant Data 0")
plot(sim_xfdata.int[1:50], type="b")
title("Simulated Int-Transplant Data 0")
plot(sim_xfdata.post[1:50], type="b")
title("Simulated Post-Transplant Data 0")
plot(sim_xfdata.pre[51:100], type="b")
title("Simulated Pre-Transplant Data 1")
plot(sim_xfdata.int[51:100], type="b")
title("Simulated Int-Transplant Data 1")
plot(sim_xfdata.post[51:100], type="b")
title("Simulated Post-Transplant Data 1")


## ========== PLS MODELS ON THE DATA =========
library(caret)
library(ROCR)
library(pROC)


sim.fda.pre <- fdata(sim_xfdata.pre)
sim.fda.int <- fdata(sim_xfdata.int)
sim.fda.post <- fdata(sim_xfdata.post)

#response variable
y = c(rep(0, 50), rep(1, 50))
yy = as.factor(y)

# PRE
#max 7 components
pls_model_pre <- fregre.pls(sim.fda.pre, y , l = 3, lambda = 0)
res.pre <- round(pls_model_pre$fitted.values)
auc.pre <- auc(y, res.pre)
rocplot(as.numeric(res.pre), y)

# INT

pls_model_int <- fregre.pls(sim.fda.int, y , l = 3, lambda = 0)
res.int <- round(pls_model_int$fitted.values)
auc.int <- auc(y, res.int)
rocplot(as.numeric(res.int), y)

# POST

pls_model_post <- fregre.pls(sim.fda.post, y , l = 3, lambda = 0)
res.post <- round(pls_model_post$fitted.values)
auc.post <- auc(y, res.post)
rocplot(as.numeric(res.post), y)

## MULTIVARIATE MODEL

fdatapls.pre <- fdata2pls(sim.fda.pre, y, ncomp = 2, lambda = 0)
T.pre <- fdatapls.pre$x

fdatapls.int <- fdata2pls(sim.fda.int, y, ncomp = 1, lambda = 0)
T.int <- fdatapls.int$x

fdatapls.post <- fdata2pls(sim.fda.post, y, ncomp = 1, lambda = 0)
T.post <- fdatapls.post$x

#linking two X matrices together

lm.mod <- glm(T.int[,1] ~ T.pre)
#lm.mod2 <- glm(T.int[,2] ~ T.pre)
H_t <-  cbind(lm.mod$residuals) #, lm.mod2$residuals)

lm.mod.post <- lm(T.post ~ T.pre)

H_t2 <- lm.mod.post$residuals

postbyint <- lm(H_t2 ~ H_t)
Htt <- postbyint$residuals

tildeX3 <- cbind(T.pre, H_t, Htt)
fin.mod3.residuals <- glm(y ~ tildeX3, family='binomial', control = list(maxit = 50), na.action = na.pass)

outcome.fin.res3 <- as.factor(round(fin.mod3.residuals$fitted.values))
fin.res.tab3 <- confusionMatrix(outcome.fin.res3, as.factor(y))$table
rocplot(as.numeric(outcome.fin.res3), as.numeric(y))
auc.mixed <- auc(y, as.numeric(outcome.fin.res3))


auc_pre_list[j] = auc.pre
auc_int_list[j] = auc.int
auc_mixed_list[j] = auc.mixed


### ===== 500 TIMES ====

auc_pre_list = rep(NA, 500)
auc_int_list = rep(NA, 500)
auc_post_list = rep(NA, 500)
auc_mixed_list = rep(NA, 500)
auc_mixed2_list = rep(NA, 500)

for(j in 1:500){
  
  #for group 0
  
  for (i in 1:5){
    var11_0 <- var_pre_0[i] 
    var22_0 <- var_int_0[i]
    var33_0 <- var_post_0[i]
    cov12_0 <- cov(sc_pre[c.black,i], sc_int[c.black,i])
    cov13_0 <- cov(sc_pre[c.black,i][-8], sc_post[c.black.post,i])
    cov23_0 <- cov(sc_int[c.black,i][-8], sc_post[c.black.post,i])
    
    n0 <- 50
    mmu_0 <- matrix(c(mean_pre_0[i], mean_int_0[i], mean_post_0[i]), ncol=1)
    mcov_0 <- matrix(c(var11_0, cov12_0, cov13_0, 
                       cov12_0, var22_0, cov23_0,
                       cov13_0, cov23_0, var33_0), ncol=3)
    bivar <- mvrnorm(n = n0, mu = mmu_0, Sigma = mcov_0)  # Random sample from bivariate normal distribution
    name <- paste('zeroscores',i, sep = "", collapse = NULL)
    assign(name, bivar)
  }
  
  #for group 1
  
  for (i in 1:5){
    var11_1 <- var_pre_1[i] 
    var22_1 <- var_int_1[i]
    var33_1 <- var_post_1[i]
    cov12_1 <- cov(sc_pre[c.red,i], sc_int[c.red,i])
    cov13_1 <- cov(sc_pre[c.red,i][-7], sc_post[c.red.post,i])
    cov23_1 <- cov(sc_int[c.red,i][-7], sc_post[c.red.post,i])
    
    n1 <- 50
    mmu_1 <- matrix(c(mean_pre_1[i], mean_int_1[i], mean_post_1[i]), ncol=1)
    mcov_1 <- matrix(c(var11_1, cov12_1, cov13_1, 
                       cov12_1, var22_1, cov23_1,
                       cov13_1, cov23_1, var33_1), ncol=3)
    bivar <- mvrnorm(n = n1, mu = mmu_1, Sigma = mcov_1)  # Random sample from bivariate normal distribution
    name <- paste('onescores',i, sep = "", collapse = NULL)
    assign(name, bivar)
  }
  
  ## ========= TURN THE VECTORS INTO MATRICES FOR LATEER ==============
  library(wakefield)
  
  ##attempt to make a mixed vector of scores 1
  
  
  random_vec <- r_sample_logical(100, prob = NULL, name = "Logical")
  true_values <- which(random_vec == TRUE)
  
  mixedscores1 <- zeroscores1
  mixedscores1[true_values] = onescores1[true_values]
  
  averagescore = (2*zeroscores1+onescores1)/3
  
  scorelist_zero <- cbind(zeroscores1[,1]*0.05, zeroscores2[,1], zeroscores3[,1], zeroscores4[,1], zeroscores5[,1],
                          zeroscores1[,2]*0.05, zeroscores2[,2], zeroscores3[,2], zeroscores4[,2], zeroscores5[,2],
                          zeroscores1[,3]*0.05, zeroscores2[,3], zeroscores3[,3], zeroscores4[,3], zeroscores5[,3])
  
  scorelist_one <- cbind(onescores1[,1]*0.05, onescores2[,1], onescores3[,1], onescores4[,1], onescores5[,1],
                         onescores1[,2]*0.05, onescores2[,2], onescores3[,2], onescores4[,2], onescores5[,2],
                         onescores1[,3]*0.05, onescores2[,3], onescores3[,3], onescores4[,3], onescores5[,3])
  
  for (i in 1:5){
    T0_pre <- replicate(52, scorelist_zero[,i])
    T0_int <- replicate(52, scorelist_zero[,i+5])
    T0_post <- replicate(52, scorelist_zero[,i+10])
    T1_pre <-  replicate(52, scorelist_one[,i])
    T1_int <-  replicate(52, scorelist_one[,i+5])
    T1_post <-  replicate(52, scorelist_one[,i+10])
    
    name_pre <- paste('Tpre',i, sep = "", collapse = NULL)
    name_int <- paste('Tint',i, sep = "", collapse = NULL)
    name_post <- paste('Tpost',i, sep = "", collapse = NULL)
    
    assign(name_pre, rbind(T0_pre,T1_pre))
    assign(name_int, rbind(T0_int,T1_int))
    assign(name_post, rbind(T0_post,T1_post))
  }  
  
  
  ## ============ USE THE MATRICES TO GENERATE DATA ========
  
  fhr_pre <- fdata(hr_pre)
  fhr_int <- fdata(hr_int)
  fhr_post <- fdata(hr_post)
  
  Tlist_pre <- cbind(Tpre1, Tpre2, Tpre3, Tpre4, Tpre5)
  Tlist_int <- cbind(Tint1, Tint2, Tint3, Tint4, Tint5)
  Tlist_post <- cbind(Tpost1, Tpost2, Tpost3, Tpost4, Tpost5)
  
  #add noise
  
  noisepre <- rnorm(ncol(Tlist_pre)*nrow(Tlist_pre))
  noiseint <- rnorm(ncol(Tlist_int)*nrow(Tlist_int))
  noisepost <- rnorm(ncol(Tlist_post)*nrow(Tlist_post))
  
  scale <- function(noise){(noise-min(noise))/(max(noise)-min(noise))*(1.2e-5-2e-6)+2e-6}
  w.noisepre <- scale(noisepre)
  w.noiseint <- scale(noiseint)
  w.noisepost <- scale(noisepost)
  
  npre_mat <- matrix(NA, nrow=100, ncol=260)
  for (i in 1:100){
    for (k in 1:260)
      npre_mat[i,k] = w.noisepre[(i-1)*260+k]
  }
  
  nint_mat <- matrix(NA, nrow=100, ncol=260)
  for (i in 1:100){
    for (k in 1:260)
      nint_mat[i,k] = w.noiseint[(i-1)*260+k]
  }
  
  npost_mat <- matrix(NA, nrow=100, ncol=260)
  for (i in 1:100){
    for (k in 1:260)
      npost_mat[i,k] = w.noisepost[(i-1)*260+k]
  }
  
  Tlist_pre = Tlist_pre + npre_mat
  Tlist_int = Tlist_int + nint_mat
  Tlist_post = Tlist_post + nint_mat
  
  pre_sim <- (fhr_pre$data[1,])*t(Tpre1)
  int_sim <- (fhr_int$data[1,])*t(Tint1)
  post_sim <- (fhr_post$data[1,])*t(Tpost1)
  
  for (i in 2:5){
    pre_sim = pre_sim + (fhr_pre$data[i,])*t(Tlist_pre[,i:(i+51)])
    int_sim = int_sim + (fhr_int$data[i,])*t(Tlist_int[,i:(i+51)])
    post_sim = post_sim + (fhr_post$data[i,])*t(Tlist_post[,i:(i+51)])
  }
  
  data_mean_int <- c(0,colMeans(kt.int[,-1]),0)
  data_mean_pre <- c(0, colMeans(kt.pre[,-1]),0)
  data_mean_post <- c(0, colMeans(kt.post[,-1]),0)
  #
  pre_sim = pre_sim + data_mean_pre
  int_sim = int_sim + data_mean_int
  post_sim = post_sim + data_mean_post
  
  
  ## ============ CREATE SPLINES FOR THE DATA ========
  
  #set variables
  bin <- seq(0,8, by=8/51)
  nobs = length(bin)
  knots <- seq(0,8, by = 8/51)
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
  
  lambda= 10^(-15)
  fdParobj = fdPar(fdobj=basisobj, lambda=lambda)
  
  sim_xfdata.pre <- smooth.basis(bin, pre_sim, fdParobj)$fd
  sim_xfdata.int <- smooth.basis(bin, int_sim, fdParobj)$fd
  sim_xfdata.post <- smooth.basis(bin, post_sim, fdParobj)$fd
  
  ## ========== PLS MODELS ON THE DATA =========
  
  sim.fda.pre <- fdata(sim_xfdata.pre)
  sim.fda.int <- fdata(sim_xfdata.int)
  sim.fda.post <- fdata(sim_xfdata.post)
  
  #response variable
  y = c(rep(0, 50), rep(1, 50))
  yy = as.factor(y)
  
  # PRE
  #max 7 components
  pls_model_pre <- fregre.pls(sim.fda.pre, y , l = 4, lambda = 0)
  res.pre <- round(pls_model_pre$fitted.values)
  auc.pre <- auc(y, res.pre)
  # rocplot(as.numeric(res.pre), y)
  
  # INT
  
  pls_model_int <- fregre.pls(sim.fda.int, y , l = 4, lambda = 0)
  res.int <- round(pls_model_int$fitted.values)
  auc.int <- auc(y, res.int)
  #rocplot(as.numeric(res.int), y)
  
  # POST
  
  pls_model_post <- fregre.pls(sim.fda.post, y , l = 4, lambda = 0)
  res.post <- round(pls_model_post$fitted.values)
  auc.post <- auc(y, res.post)
  # rocplot(as.numeric(res.post), y)
  
  ## MULTIVARIATE MODEL
  
  fdatapls.pre <- fdata2pls(sim.fda.pre, y, ncomp = 2, lambda = 0)
  T.pre <- fdatapls.pre$x
  
  fdatapls.int <- fdata2pls(sim.fda.int, y, ncomp = 2, lambda = 0)
  T.int <- fdatapls.int$x
  
  #linking two X matrices together
  
  lm.mod <- glm(T.int[,1] ~ T.pre)
  lm.mod2 <- glm(T.int[,2] ~ T.pre)
  H_t <-  cbind(lm.mod$residuals, lm.mod2$residuals)
  
  tildeX <- as.matrix(cbind(T.pre, H_t))
  fin.mod.residuals <- glm(y ~ tildeX, family='binomial', control = list(maxit = 50))
  
  outcome.fin.res <- round(fin.mod.residuals$fitted.values)
  auc.mixed2 <- auc(y, outcome.fin.res)
  
  
  ## MULTIVARIATE MODEL
  
  fdatapls.pre <- fdata2pls(sim.fda.pre, y, ncomp = 2, lambda = 0)
  T.pre <- fdatapls.pre$x
  
  fdatapls.int <- fdata2pls(sim.fda.int, y, ncomp = 1, lambda = 0)
  T.int <- fdatapls.int$x
  
  fdatapls.post <- fdata2pls(sim.fda.post, y, ncomp = 1, lambda = 0)
  T.post <- fdatapls.post$x
  
  #linking two X matrices together
  
  lm.mod <- glm(T.int[,1] ~ T.pre)
  #lm.mod2 <- glm(T.int[,2] ~ T.pre)
  H_t <-  cbind(lm.mod$residuals) #, lm.mod2$residuals)
  
  lm.mod.post <- lm(T.post ~ T.pre)
  
  H_t2 <- lm.mod.post$residuals
  
  postbyint <- lm(H_t2 ~ H_t)
  Htt <- postbyint$residuals
  
  tildeX3 <- cbind(T.pre, H_t, Htt)
  fin.mod3.residuals <- glm(y ~ tildeX3, family='binomial', control = list(maxit = 50), na.action = na.pass)
  
  outcome.fin.res3 <- as.factor(round(fin.mod3.residuals$fitted.values))
  fin.res.tab3 <- confusionMatrix(outcome.fin.res3, as.factor(y))$table
  #rocplot(as.numeric(outcome.fin.res3), as.numeric(y))
  auc.mixed <- auc(y, as.numeric(outcome.fin.res3))
  
  
  auc_pre_list[j] = auc.pre
  auc_int_list[j] = auc.int
  auc_post_list[j] = auc.post
  auc_mixed2_list[j] = auc.mixed2
  auc_mixed_list[j] = auc.mixed
}

par(mfrow=c(1,4))
mean(auc_pre_list)
var(auc_pre_list)
hist(auc_pre_list)

mean(auc_int_list)
var(auc_int_list)
hist(auc_int_list)

mean(auc_post_list)
var(auc_post_list)
hist(auc_post_list)

mean(auc_mixed_list)
var(auc_mixed_list)
hist(auc_mixed_list)

boxplot(auc_pre_list, auc_int_list, auc_post_list, auc_mixed2_list, auc_mixed_list,
        main = "Multiple boxplots for comparision",
        at = c(1,2,3,4,5),
        names = c("ozone", "normal", "temp", "normal", 'poop'),
        las = 2,
        col = c("orange","red"),
        border = "brown")


#set up dataset


simulation_1correct_1small = cbind(auc_pre_list, auc_int_list, auc_post_list, auc_mixed2_list, auc_mixed_list)

sum(((auc_pre_list < auc_mixed_list) & (auc_int_list< auc_mixed_list) & (auc_post_list< auc_mixed_list))==TRUE)

sum(((auc_pre_list < auc_mixed2_list) & (auc_int_list< auc_mixed2_list) & (auc_post_list< auc_mixed2_list))==TRUE)
##
## --------------- BOXPLOTS ----------------
##


library(broman)
brocolors("crayons")


par(mfrow=c(3,1))

boxplot(data_first_score_small[,1], data_first_score_small[,2], data_first_score_small[,3], data_first_score_small[,4],
        data_1score_wrong[,1], data_1score_wrong[,2], data_1score_wrong[,3], data_1score_wrong[,4],
        data_15score_wrong_1smol[,1], data_15score_wrong_1smol[,2], data_15score_wrong_1smol[,3], data_15score_wrong_1smol[,4],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:4, 6:9, 11:14),
        names = rep(c("PRE", "INT", "POST", "MUL"),3),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc"),
        horizontal = TRUE,
        border = c( rep('brown',4),  rep('#066F0A',4), rep('midnight blue',4)))

abline(v=c(0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), col="#E0E0E0")

boxplot(data_first_score_small[,1], data_first_score_small[,2], data_first_score_small[,3], data_first_score_small[,4],
        data_1score_wrong[,1], data_1score_wrong[,2], data_1score_wrong[,3], data_1score_wrong[,4],
        data_15score_wrong_1smol[,1], data_15score_wrong_1smol[,2], data_15score_wrong_1smol[,3], data_15score_wrong_1smol[,4],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:4, 6:9, 11:14),
        names = rep(c("PRE", "INT", "POST", "MUL"),3),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc"),
        horizontal = TRUE,
        border = c( rep('brown',4),  rep('#066F0A',4), rep('midnight blue',4)), add = TRUE, ann = FALSE)

legend("bottomleft", inset=.02, legend=c("Sim A", "Sim B", "Sim C"),
       fill=c("#faa76c", "#14994B", "#6699cc"), cex=1)
#########################

boxplot(data_first_score_small[,1], data_first_score_small[,2], data_first_score_small[,3], data_first_score_small[,4],
        data_1score_small[,1], data_1score_small[,2], data_1score_small[,3], data_1score_small[,4],
        data_15score_wrong_1smol[,1], data_15score_wrong_1smol[,2], data_15score_wrong_1smol[,3], data_15score_wrong_1smol[,4],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:4, 6:9, 11:14),
        names = rep(c("PRE", "INT", "POST", "MUL"),3),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc"),
        horizontal = TRUE,
        border = c( rep('brown',4),  rep('#066F0A',4), rep('midnight blue',4)))

abline(v=c(0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), col="#E0E0E0")

boxplot(data_first_score_small[,1], data_first_score_small[,2], data_first_score_small[,3], data_first_score_small[,4],
        data_1score_small[,1], data_1score_small[,2], data_1score_small[,3], data_1score_small[,4],
        data_15score_wrong_1smol[,1], data_15score_wrong_1smol[,2], data_15score_wrong_1smol[,3], data_15score_wrong_1smol[,4],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:4, 6:9, 11:14),
        names = rep(c("pre", "int", "post", "mul"),3),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc"),
        horizontal = TRUE,
        border = c( rep('brown',4),  rep('#066F0A',4), rep('midnight blue',4)), add = TRUE, ann = FALSE)

legend("bottomleft", inset=.02, legend=c("Sim A", "Sim B", "Sim C"),
       fill=c("#faa76c", "#14994B", "#6699cc"), cex=1)


##########################################

#########################
SIM_A = simulation_1correct_1small
SIM_B = simulation_1wrong_1small
SIM_C = simulation_13wrong_1small
SIM_D = simulation_135wrong_1small

boxplot(SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        SIM_B[,1], SIM_B[,2], SIM_B[,3], SIM_B[,4], SIM_B[,5],
        SIM_C[,1], SIM_C[,2], SIM_C[,3], SIM_C[,4], SIM_C[,5],
        SIM_D[,1], SIM_D[,2], SIM_D[,3], SIM_D[,4], SIM_D[,5],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:5, 7:11, 13:17, 19:23),
        names = rep(c("PRE", "INT", "POST","2 MUL", "3 MUL"),4),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc", "#6699cc",
                "#ADA1E0", "#ADA1E0", "#ADA1E0", "#7A6CEF", "#7A6CEF"),
        horizontal = TRUE,
        border = c( rep('brown',5),  rep('#066F0A',5), rep('midnight blue',5), rep("#2E066D", 5)))

abline(v=c(0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), col="#E0E0E0")

boxplot(SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        SIM_B[,1], SIM_B[,2], SIM_B[,3], SIM_B[,4], SIM_B[,5],
        SIM_C[,1], SIM_C[,2], SIM_C[,3], SIM_C[,4], SIM_C[,5],
        SIM_D[,1], SIM_D[,2], SIM_D[,3], SIM_D[,4], SIM_D[,5],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:5, 7:11, 13:17, 19:23),
        names = rep(c("PRE", "INT", "POST","2 MUL", "3 MUL"),4),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc", "#6699cc",
                "#ADA1E0", "#ADA1E0", "#ADA1E0", "#7A6CEF", "#7A6CEF"),
        horizontal = TRUE,
        border = c( rep('brown',5),  rep('#066F0A',5), rep('midnight blue',5), rep("#2E066D", 5)), add = TRUE, ann = FALSE)

legend("bottomleft", inset=.02, legend=c("Sim A", "Sim B", "Sim C", "Sim D"),
       fill=c("#faa76c", "#14994B", "#6699cc", "#7A6CEF"), cex=1)


##########################################

#########################
SIM_A = simulation_1correct_1small
SIM_B = simulation_1wrong_1small
SIM_C = simulation_13wrong_1small
SIM_D = simulation_135wrong_1small

boxplot(SIM_D[,1], SIM_D[,2], SIM_D[,3], SIM_D[,4], SIM_D[,5],
        SIM_C[,1], SIM_C[,2], SIM_C[,3], SIM_C[,4], SIM_C[,5],
        SIM_B[,1], SIM_B[,2], SIM_B[,3], SIM_B[,4], SIM_B[,5],
        SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:5, 7:11, 13:17, 19:23),
        names = rep(c("pre", "int", "post","mul. (2)", "mul. (3)"),4),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc", "#6699cc",
                "#ADA1E0", "#ADA1E0", "#ADA1E0", "#7A6CEF", "#7A6CEF"),
        horizontal = TRUE,
        border = c( rep('brown',5),  rep('#066F0A',5), rep('midnight blue',5), rep("#2E066D", 5)))

abline(v=seq(0.5,1, by=0.025), col="#E0E0E0")

boxplot(SIM_D[,1], SIM_D[,2], SIM_D[,3], SIM_D[,4], SIM_D[,5],
        SIM_C[,1], SIM_C[,2], SIM_C[,3], SIM_C[,4], SIM_C[,5],
        SIM_B[,1], SIM_B[,2], SIM_B[,3], SIM_B[,4], SIM_B[,5],
        SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        main = "Multiple Boxplots for Simulated Data",
        at = c(1:5, 7:11, 13:17, 19:23),
        names = rep(c("pre", "int", "post","mul. (2)", "mul. (3)"),4),
        las = 2,
        col = c("#fce883","#fce883","#fce883", "#faa76c", "#faa76c",
                "#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B",
                "#c5d0e6" , "#c5d0e6" ,"#c5d0e6" , "#6699cc", "#6699cc",
                "#ADA1E0", "#ADA1E0", "#ADA1E0", "#7A6CEF", "#7A6CEF"),
        horizontal = TRUE,
        border = c( rep('brown',5),  rep('#066F0A',5), rep('midnight blue',5), rep("#2E066D", 5)), add = TRUE, ann = FALSE)

legend("topleft", inset=.02, legend=c("Scenario A", "Scenario B", "Scenario C", "Scenario D"),
       fill=c("#7A6CEF",  "#6699cc", "#14994B","#faa76c" ), cex=1)

#### ===============================


boxplot(SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        main = "AUC distribution for Simulated Data",
        at = c(1:5),
        names = c("pre", "int", "post","mul. (2)", "mul. (3)"),
        las = 2,
        col = c("#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B"),
        horizontal = TRUE,
        border = c(rep('#066F0A',5)))

abline(v=seq(0.5,1, by=0.025), col="#E0E0E0")


boxplot(SIM_A[,1], SIM_A[,2], SIM_A[,3], SIM_A[,4], SIM_A[,5],
        main = "AUC distribution for Simulated Data",
        at = c(1:5),
        names = c("pre", "int", "post","mul. (2)", "mul. (3)"),
        las = 2,
        col = c("#AFDD91"   ,"#AFDD91","#AFDD91"  , "#14994B", "#14994B"),
        horizontal = TRUE,
        border = c(rep('#066F0A',5)), add=TRUE, ann=FALSE)


