#SNR for functional logistic regression


### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))


library(FADPclust)
library(fda)
library(tidyverse)
library(tictoc)
library(randomForest)
library(precrec)
library(ROCit)
library(caret)

library(LocationSmoothedBART)

library(PRROC)
library(pROC)

### Clear Environment ###

##################################################################

badk <- read_csv("ORGAN DATA/Bad_Kidney_2Dfunc_Sept8.csv")
goodk <- read_csv("ORGAN DATA/Good_Kidney_2Dfunc_Sept8.csv")

colSums(is.na(goodk))
goodk[is.na(goodk)] <- mean(goodk$sphericity, na.rm=TRUE)

var_names = names(goodk)[-(1:2)]

#Only want some of the functional predictors to keep
var_names[c(1:10, 12, 13, 15, 17, 18)]

name_vec = c(1:10, 12, 13, 15, 17, 18)

var_names2 = var_names[name_vec]
num_predictors = length(var_names2)

#Data cleaning
#colSums(is.na(goodk))
#goodk[is.na(goodk)] <- mean(goodk$sphericity1, na.rm=TRUE)
#colSums(is.na(badk))


set_good <- split(goodk,goodk$name)
set_bad <- split(badk,badk$name)


#Extracts the coefficients at specif intervals across all funcitons. This case gives 11 values
#per shape feature

grid = seq(.05, from=0, to=1)
numPoinst=length(grid)+2
bb <- create.bspline.basis(rangeval = c(0,1), breaks = grid, nbasis = numPoinst)
plot(bb)

data_set <- NULL

num_predictors = length(var_names2)



for (j in 1:length(set_good)) {
  data_temp <- NULL

  for (i in name_vec ) {
    f <- Data2fd(set_good[[j]][[i+2]], basisobj = bb)
    data_temp <- c(data_temp, f$coefs[,1])
    #print(paste(j, '---', i))
  }
  data_set <- rbind(data_set, data_temp)
}

for (j in 1:length(set_bad)) {
  data_temp <- NULL

  for (i in name_vec ) {
    f <- Data2fd(set_bad[[j]][[i+2]], basisobj = bb)
    data_temp <- c(data_temp, f$coefs[,1])
  }
  data_set <- rbind(data_set, data_temp)

}

data_set
#Logistic Regression
outcomes = c(rep(0,260), rep(1, 52))



#create Train and Test sets
df = as_tibble(data_set)

outcomes = c(rep(0,260), rep(1, 52))
outcomes_class = c(rep("Acceptable",260), rep("Unacceptable", 52))
true_y = outcomes

# Testing out FBART
num_intvls = dim(data_set)[2]/num_predictors

points = data_set

n_point =  dim(data_set)[2]/num_predictors



#ALL SET UP
tic()
set.seed(1)

flog_num_predictors = 12

tempbasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)
grid <- seq(0,1,1/(length(points[1,1:n_point])-1))

templist      = vector("list", flog_num_predictors+1)
templist[[1]] = rep(1,  nrow(points))

for (i in 0:(flog_num_predictors-1) ) {
  tempSmooth= smooth.basis(grid, t(points[,(n_point*i+1):(n_point*i+n_point)]), tempbasis)
  templist[[i+2]] = tempSmooth$fd


}



conbasis = create.constant.basis(c(0,1))
betabasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)

betalist = vector("list",flog_num_predictors+1)
betalist[[1]] = conbasis

for (i in 2:(flog_num_predictors+1) ) {
  betalist[[i]] = betabasis

}


fRegressList2 = fRegress(log(abs(true_y-1e-1)/(1-abs(true_y-1e-1))), templist, betalist)
yhat <- exp(fRegressList2$yhatfdobj)/(1+exp(fRegressList2$yhatfdobj))


plot(yhat)

###############
#get the SNR from logistic functional regression

true_y

mean(true_y)/ sd(true_y)


var(yhat)/var(true_y-yhat)






# Compute predicted probabilities
predicted_probs <- yhat

# Calculate the signal strength (variance of the logits)
logits <- log(predicted_probs / (1 - predicted_probs))

signal_variance <- var(logits)

# Calculate average noise (variance of binomial distribution for each predicted probability)
noise_variance <- mean(predicted_probs * (1 - predicted_probs))


# Calculate SNR
snr <- signal_variance / noise_variance

# Print SNR
print(paste("Estimated SNR:", snr))


