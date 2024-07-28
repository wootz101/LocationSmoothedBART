# Plots for labeled kidney shape statistics
# Getting functional Data from Kidney data-set for FBART
# 10 REPLICATES OF CV

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





#######################################################################################
#######################################################################################

num_intvls


a_int=0.5
b_int=2

ntree = 100
ndpost = 200

binaryOffset = 25
rho_fp = NULL
rho_int = NULL
sigma_int = 2



nskip = 100
kk= 2
power = 2
base = 0.5


##################################
#Replicate holders


AUC_cv_r = NULL
PR_cv_r = NULL
F1_cv_r = NULL

Sens_cv_r = NULL
Spec_cv_r = NULL

Sens_y_cv_r = NULL
Spec_y_cv_r = NULL

youdens_cv_r = NULL

#Time these models
flog_time_all=list()
fbart_time_all=list()
bart_time_all=list()
dart_time_all=list()
rf_time_all=list()

# Begin Replicates
for(cv_r in 1:2){
  #######################################################
  #Begin the 5 fold cross validation


  numfold=10
  aucFLR <- 1:numfold
  betas_fold <- list()
  temp_train_fold <- list()
  temp_test_fold <- list()

  varProb <- NULL
  kfold <- cbind(true_y, 1:numfold)[,2]


  # For logistic Regression
  #ds2 <- as.data.frame(data_set)

  set.seed(cv_r)
  shuffled_data = sample(1:nrow(df))
  ds2_shuf <-  df[shuffled_data,]
  true_y_shuf <- true_y[shuffled_data]


  # The truth
  all_test_y <- NULL
  all_train_y <- NULL

  #functional Regresion
  all_test_y_obj <- NULL
  all_train_y_obj <- NULL

  #FuncBART
  all_test_y_fb <- NULL
  all_train_y_fb <- NULL

  #BART
  all_test_y_BART <- NULL
  all_train_y_BART <- NULL

  #DART
  all_test_y_DART <- NULL
  all_train_y_DART <- NULL

  #RF
  all_test_y_RF <- NULL
  all_train_y_RF <- NULL






  for (k in 1:numfold) {


    ################################################
    #Logistic Regression

    train <- ds2_shuf[kfold != k,]
    test <- ds2_shuf[kfold == k,]

    y_train = true_y_shuf[kfold != k]
    y_test = true_y_shuf[kfold == k]



    # Random Shuffled

    #sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.8,0.2))
    #train  <- df[sample, ]
    #test   <- df[!sample, ]
    #y_train = true_y[sample]
    #y_test = true_y[!sample]
    ######################################################################
    #FunctionalBART

    tic()
    set.seed(1)
    m0 = p_lsBART(num_predictors = num_predictors,
                  num_intvls =num_intvls,
                  x.train = as.matrix(train) ,
                  y.train = y_train,
                  x.test =  as.matrix(test),
                  #rho=rho,
                  #rho_fp=rho,
                  rho_int=rho_int,
                  sigma_int = sigma_int,
                  a_int=a_int,
                  b_int=b_int,
                  k = kk,
                  power = power,
                  base = base,
                  binaryOffset = binaryOffset,
                  sparse = TRUE,
                  ntree= ntree,
                  ndpost=ndpost,
                  nskip = nskip,
                  dart_fp = TRUE,
                  dart_int = TRUE)

    all_test_y_fb = c(all_test_y_fb, m0$prob.test.mean)
    all_train_y_fb = c(all_train_y_fb, m0$prob.train.mean)
    print("______________________________________________________")
    print(paste("Finished FuncBART --", cv_r, "--",k))

    time1 =toc()
    fbart_time_all[(cv_r-1)*10+k] = str_split(time1$callback_msg, ' ')[[1]][1]


    ###########################################################################



    #ALL SET UP
    tic()
    set.seed(1)

    flog_num_predictors = 8

    tempbasis = create.bspline.basis(rangeval = c(0,1), nbasis=8)
    grid <- seq(0,1,1/(length(train[1,1:n_point])-1))

    templist      = vector("list", flog_num_predictors+1)
    templist[[1]] = rep(1,  nrow(train))

    for (i in 0:(flog_num_predictors-1) ) {
      tempSmooth= smooth.basis(grid, t(train[,(n_point*i+1):(n_point*i+n_point)]), tempbasis)
      templist[[i+2]] = tempSmooth$fd


    }

    #t(train[,(n_point*0+1):(n_point*0+n_point)])

    conbasis = create.constant.basis(c(0,1))
    betabasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)

    betalist = vector("list",flog_num_predictors+1)
    betalist[[1]] = conbasis

    for (i in 2:(flog_num_predictors+1) ) {
      betalist[[i]] = betabasis

    }


    #outcomes_train = outcomes_shuf[kfold != k]
    #outcomes_test = outcomes_shuf[kfold == k]
    #y_train = cc[sample]*5+rnorm(sum(sample), mean=0, sd=2)


    #fRegressList2 = fRegress( y_train , templist, betalist)
    fRegressList2 = fRegress(log(abs(y_train-1e-1)/(1-abs(y_train-1e-1))), templist, betalist)
    yhat <- exp(fRegressList2$yhatfdobj)/(1+exp(fRegressList2$yhatfdobj))

    #plot(fRegressList2$yhatfdobj)


    #fRegressList2$betaestlist[[2]]
    #plot(fRegressList2$betaestlist[[2]])
    #plot(fRegressList2$betaestlist[[4]])

    # Test set

    #Create fucntional data for test set
    testlist      = vector("list", num_predictors+1)
    testlist[[1]] = rep(1,  nrow(test))
    grid <- seq(0,1,1/(length(test[1,1:n_point])-1))

    for (i in 0:(flog_num_predictors-1)) {
      testSmooth= smooth.basis(grid, t(test[  ,(n_point*i+1):(n_point*i+n_point) ]), tempbasis)
      testlist[[i+2]] = testSmooth$fd

    }

    blist <- fRegressList2$betaestlist

    #Use the coefficients from the training set to predict from unseeen
    #test set. Call it b_coefs

    b_coefs <- NULL
    b_coefs <- rbind(b_coefs, rep(blist[[1]]$fd$coefs, length(y_test)))
    for (i in 2:(flog_num_predictors+1)) {
      b_coefs <- rbind(b_coefs, inprod(blist[[i]]$fd, testlist[[i]]))
    }

    #Answer Prediction
    yhat_test_fReg = colSums(b_coefs)
    yhat_test <- exp(yhat_test_fReg)/(1+exp(yhat_test_fReg))

    all_test_y_obj <- c(all_test_y_obj, yhat_test)
    all_train_y_obj <- c(all_train_y_obj, yhat)

    time1 =toc()
    flog_time_all[(cv_r-1)*10+k] = str_split(time1$callback_msg, ' ')[[1]][1]





    print("______________________________________________________")
    print(paste("Finished Log Reg --", cv_r, "--", k))

    ###########################################################################

    ########################################################################
    #BART
    tic()
    set.seed(1)
    m2 = pbart(x.train = as.matrix(train) , y.train = y_train, ntree = ntree,
               x.test =  as.matrix(test),  sparse = FALSE, binaryOffset = NULL,
               ndpost=ndpost, nskip = nskip)

    all_test_y_BART = c(all_test_y_BART, m2$prob.test.mean)
    all_train_y_BART = c(all_train_y_BART, m2$prob.train.mean)
    print("______________________________________________________")
    print(paste("Finished  BART --",cv_r, "--", k))

    time1 =toc()
    bart_time_all[(cv_r-1)*10+k] = str_split(time1$callback_msg, ' ')[[1]][1]

    ########################################################################
    #DART
    tic()
    set.seed(1)
    m3 = pbart(x.train = as.matrix(train) , y.train = y_train, ntree = ntree,
               x.test =  as.matrix(test),  sparse = TRUE, binaryOffset = NULL,
               ndpost=ndpost, nskip = nskip)

    all_test_y_DART = c(all_test_y_DART, m3$prob.test.mean)
    all_train_y_DART = c(all_train_y_DART, m3$prob.train.mean)
    print("______________________________________________________")
    print(paste("Finished  DART --",cv_r, "--", k))

    time1 =toc()
    dart_time_all[(cv_r-1)*10+k] = str_split(time1$callback_msg, ' ')[[1]][1]

    ########################################################################
    #Testing using Radnom Forest
    tic()
    set.seed(1)

    m5 <- randomForest(y_train ~ ., data=train, importance=TRUE, proximity=TRUE, mtry=ndpost, ntree=ntree)
    p5 = predict(m5, test)

    all_test_y_RF = c(all_test_y_RF, p5)
    all_train_y_RF = c(all_train_y_RF, m5$predicted)

    #all_test_y_RF = c(all_test_y_RF, m3$prob.test.mean+.1)
    #all_train_y_RF = c(all_train_y_RF,  m3$prob.train.mean+0.1)


    print("______________________________________________________")

    print(paste("Finished  RF --", cv_r, "--",k))
    time1 =toc()
    rf_time_all[(cv_r-1)*10+k] = str_split(time1$callback_msg, ' ')[[1]][1]


    # THE TRUTH
    all_test_y <- c(all_test_y, y_test)
    all_train_y <- c(all_train_y, y_train)

  }


  CV1 = cbind(all_test_y_obj, all_test_y_fb, all_test_y_BART,
              all_test_y_DART, all_test_y_RF, all_test_y )


  #msmdat2 <- mmdata(CV1[,1:4],
  #                 c(CV1[,5]), modnames = c("LogReg", "funcBART", "BART", "DART" ))

  # msmdat2 <- mmdata(CV1[ ,1:5],
  # c(CV1[,6]), modnames = c("LogReg", "funcBART", "BART", "DART", "RF" ))

  #mscurves <- evalmod(msmdat2)

  #autoplot(mscurves)
  #aa <- auc(mscurves)
  # Results
  #resultsDF = data.frame(matrix(round(aa$aucs, digits=3), ncol=5))
  # names(resultsDF) = c("LogReg","funcBART", "BART", "DART", "RF"  )
  #resultsDF



  AUC_cv_r = data.frame(rbind(AUC_cv_r , c(roc.curve( CV1[CV1[ ,6]==1 ,1], CV1[ CV1[ ,6]==0 ,1])$auc,
                                           roc.curve( CV1[CV1[ ,6]==1 ,2], CV1[ CV1[ ,6]==0 ,2])$auc,
                                           roc.curve( CV1[CV1[ ,6]==1 ,3], CV1[ CV1[ ,6]==0 ,3])$auc,
                                           roc.curve( CV1[CV1[ ,6]==1 ,4], CV1[ CV1[ ,6]==0 ,4])$auc,
                                           roc.curve( CV1[CV1[ ,6]==1 ,5], CV1[ CV1[ ,6]==0 ,5])$auc )) )

  PR_cv_r = data.frame(rbind(PR_cv_r,c(pr.curve( CV1[CV1[ ,6]==1 ,1], CV1[ CV1[ ,6]==0 ,1])$auc.integral,
                                       pr.curve( CV1[CV1[ ,6]==1 ,2], CV1[ CV1[ ,6]==0 ,2])$auc.integral,
                                       pr.curve( CV1[CV1[ ,6]==1 ,3], CV1[ CV1[ ,6]==0 ,3])$auc.integral,
                                       pr.curve( CV1[CV1[ ,6]==1 ,4], CV1[ CV1[ ,6]==0 ,4])$auc.integral,
                                       pr.curve( CV1[CV1[ ,6]==1 ,5], CV1[ CV1[ ,6]==0 ,5])$auc.integral )) )

  cm_flog = confusionMatrix( as.factor(as.integer(as.logical(CV1[,1]>0.5))),
                             as.factor(CV1[,6]), mode = "everything", positive="1")

  cm_fb = confusionMatrix( as.factor(as.integer(as.logical(CV1[,2]>0.5))),
                           as.factor(CV1[,6]), mode = "everything", positive="1")

  cm_bart = confusionMatrix( as.factor(as.integer(as.logical(CV1[,3]>0.5))),
                             as.factor(CV1[,6]), mode = "everything", positive="1")

  cm_dart = confusionMatrix( as.factor(as.integer(as.logical(CV1[,4]>0.5))),
                             as.factor(CV1[,6]), mode = "everything", positive="1")

  cm_rf = confusionMatrix( as.factor(as.integer(as.logical(CV1[,5]>0.5))),
                           as.factor(CV1[,6]), mode = "everything", positive="1")


  rocobj <- roc(CV1[,6], CV1[,1])
  cm_flog_y = coords(rocobj, "best")

  rocobj <- roc(CV1[,6], CV1[,2])
  cm_fb_y = coords(rocobj, "best")

  rocobj <- roc(CV1[,6], CV1[,3])
  cm_bart_y = coords(rocobj, "best")

  rocobj <- roc(CV1[,6], CV1[,4])
  cm_dart_y = coords(rocobj, "best")

  rocobj <- roc(CV1[,6], CV1[,5])
  cm_rf_y = coords(rocobj, "best")





  F1_cv_r = data.frame(rbind(F1_cv_r, c(cm_flog$byClass[7], cm_fb$byClass[7],
                                        cm_bart$byClass[7], cm_dart$byClass[7],
                                        cm_rf$byClass[7])) )

  Sens_cv_r  = data.frame( rbind(Sens_cv_r, c(cm_flog$byClass[1], cm_fb$byClass[1],
                                              cm_bart$byClass[1], cm_dart$byClass[1],
                                              cm_rf$byClass[1])) )

  Spec_cv_r = data.frame( rbind(Spec_cv_r, c(cm_flog$byClass[2], cm_fb$byClass[2],
                                             cm_bart$byClass[2], cm_dart$byClass[2],
                                             cm_rf$byClass[2])) )

  Sens_y_cv_r  = data.frame( rbind(Sens_y_cv_r, c(cm_flog_y$sensitivity[1], cm_fb_y$sensitivity[1],
                                                  cm_bart_y$sensitivity[1], cm_dart_y$sensitivity[1],
                                                  cm_rf_y$sensitivity[1])) )

  Spec_y_cv_r = data.frame( rbind(Spec_y_cv_r, c(cm_flog_y$specificity[1], cm_fb_y$specificity[1],
                                                 cm_bart_y$specificity[1], cm_dart_y$specificity[1],
                                                 cm_rf_y$specificity[1])) )

  youdens_cv_r = data.frame( rbind(youdens_cv_r, c(cm_flog_y$threshold[1], cm_fb_y$threshold[1],
                                                   cm_bart_y$threshold[1], cm_dart_y$threshold[1],
                                                   cm_rf_y$threshold[1])) )

  names(AUC_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(PR_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )

  names(F1_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(Sens_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(Spec_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(Sens_y_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(Spec_y_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )
  names(youdens_cv_r) = c("LogReg","lsBART", "BART", "DART", "RF"  )

}


AUC_cv_r
PR_cv_r
F1_cv_r
Sens_cv_r
Spec_cv_r
Sens_y_cv_r
Spec_y_cv_r
youdens_cv_r



round(colMeans(AUC_cv_r), 3)
round(colMeans(PR_cv_r), 3)
round(colMeans(F1_cv_r), 3)
round(colMeans(Sens_cv_r), 3)
round(colMeans(Spec_cv_r), 3)


median(AUC_cv_r[,1])
median(AUC_cv_r[,2])
median(AUC_cv_r[,3])
median(AUC_cv_r[,4])
median(AUC_cv_r[,5])


round(sd(AUC_cv_r[,1]), 3)
round(sd(AUC_cv_r[,2]), 3)
round(sd(AUC_cv_r[,3]), 3)
round(sd(AUC_cv_r[,4]), 3)
round(sd(AUC_cv_r[,5]), 3)


round(sd(PR_cv_r[,1]), 3)
round(sd(PR_cv_r[,2]), 3)
round(sd(PR_cv_r[,3]), 3)
round(sd(PR_cv_r[,4]), 3)
round(sd(PR_cv_r[,5]), 3)


round(sd(F1_cv_r[,1]), 3)
round(sd(F1_cv_r[,2]), 3)
round(sd(F1_cv_r[,3]), 3)
round(sd(F1_cv_r[,4]), 3)
round(sd(F1_cv_r[,5]), 3)

round(sd(Sens_cv_r[,1]), 3)
round(sd(Sens_cv_r[,2]), 3)
round(sd(Sens_cv_r[,3]), 3)
round(sd(Sens_cv_r[,4]), 3)
round(sd(Sens_cv_r[,5]), 3)

round(sd(Spec_cv_r[,1]), 3)
round(sd(Spec_cv_r[,2]), 3)
round(sd(Spec_cv_r[,3]), 3)
round(sd(Spec_cv_r[,4]), 3)
round(sd(Spec_cv_r[,5]), 3)


round(colMeans(Sens_y_cv_r), 3)
round(colMeans(Spec_y_cv_r), 3)
round(colMeans(youdens_cv_r), 3)


colMeans(AUC_cv_r)
colMeans(PR_cv_r)
colMeans(F1_cv_r)
colMeans(Sens_cv_r)
colMeans(Spec_cv_r)

flog_time_all[[10]]

c(flog_time_all[[1]],
  flog_time_all[[2]],
  flog_time_all[[3]],
  flog_time_all[[4]],
  flog_time_all[[5]],
  flog_time_all[[6]],
  flog_time_all[[7]],
  flog_time_all[[8]],
  flog_time_all[[9]],
  flog_time_all[[10]]
)
as.numeric(flog_time_all[[10]])

mean( as.numeric(c(flog_time_all[[1]],
                   flog_time_all[[2]],
                   flog_time_all[[3]],
                   flog_time_all[[4]],
                   flog_time_all[[5]],
                   flog_time_all[[6]],
                   flog_time_all[[7]],
                   flog_time_all[[8]],
                   flog_time_all[[9]],
                   flog_time_all[[10]]
)) )

mean( as.numeric(c(fbart_time_all[[1]],
                   fbart_time_all[[2]],
                   fbart_time_all[[3]],
                   fbart_time_all[[4]],
                   fbart_time_all[[5]],
                   fbart_time_all[[6]],
                   fbart_time_all[[7]],
                   fbart_time_all[[8]],
                   fbart_time_all[[9]],
                   fbart_time_all[[10]]
)) )

mean( as.numeric(c(bart_time_all[[1]],
                   bart_time_all[[2]],
                   bart_time_all[[3]],
                   bart_time_all[[4]],
                   bart_time_all[[5]],
                   bart_time_all[[6]],
                   bart_time_all[[7]],
                   bart_time_all[[8]],
                   bart_time_all[[9]],
                   bart_time_all[[10]]
)) )

mean( as.numeric(c(dart_time_all[[1]],
                   dart_time_all[[2]],
                   dart_time_all[[3]],
                   dart_time_all[[4]],
                   dart_time_all[[5]],
                   dart_time_all[[6]],
                   dart_time_all[[7]],
                   dart_time_all[[8]],
                   dart_time_all[[9]],
                   dart_time_all[[10]]
)) )

mean( as.numeric(c(rf_time_all[[1]],
                   rf_time_all[[2]],
                   rf_time_all[[3]],
                   rf_time_all[[4]],
                   rf_time_all[[5]],
                   rf_time_all[[6]],
                   rf_time_all[[7]],
                   rf_time_all[[8]],
                   rf_time_all[[9]],
                   rf_time_all[[10]]
)) )




sd( as.numeric(c(flog_time_all[[1]],
                 flog_time_all[[2]],
                 flog_time_all[[3]],
                 flog_time_all[[4]],
                 flog_time_all[[5]],
                 flog_time_all[[6]],
                 flog_time_all[[7]],
                 flog_time_all[[8]],
                 flog_time_all[[9]],
                 flog_time_all[[10]]
)) )

sd( as.numeric(c(fbart_time_all[[1]],
                 fbart_time_all[[2]],
                 fbart_time_all[[3]],
                 fbart_time_all[[4]],
                 fbart_time_all[[5]],
                 fbart_time_all[[6]],
                 fbart_time_all[[7]],
                 fbart_time_all[[8]],
                 fbart_time_all[[9]],
                 fbart_time_all[[10]]
)) )

sd( as.numeric(c(bart_time_all[[1]],
                 bart_time_all[[2]],
                 bart_time_all[[3]],
                 bart_time_all[[4]],
                 bart_time_all[[5]],
                 bart_time_all[[6]],
                 bart_time_all[[7]],
                 bart_time_all[[8]],
                 bart_time_all[[9]],
                 bart_time_all[[10]]
)) )

sd( as.numeric(c(dart_time_all[[1]],
                 dart_time_all[[2]],
                 dart_time_all[[3]],
                 dart_time_all[[4]],
                 dart_time_all[[5]],
                 dart_time_all[[6]],
                 dart_time_all[[7]],
                 dart_time_all[[8]],
                 dart_time_all[[9]],
                 dart_time_all[[10]]
)) )

sd( as.numeric(c(rf_time_all[[1]],
                 rf_time_all[[2]],
                 rf_time_all[[3]],
                 rf_time_all[[4]],
                 rf_time_all[[5]],
                 rf_time_all[[6]],
                 rf_time_all[[7]],
                 rf_time_all[[8]],
                 rf_time_all[[9]],
                 rf_time_all[[10]]
)) )






fbart_time_all
bart_time_all
dart_time_all
rf_time_all

mean(as.numeric(flog_time_all[[10]]))
mean(as.numeric(fbart_time_all[[10]]))
mean(as.numeric(bart_time_all[[10]]))
mean(as.numeric(dart_time_all[[10]]))
mean(as.numeric(rf_time_all[[10]]))

get_long_table <- function(table, cv_r){

  return_table = as.data.frame(cbind(
    c(table[,1],table[,2],table[,3],table[,4],table[,5]),
    c(rep("LogReg",cv_r), rep("lsBART", cv_r), rep("BART", cv_r), rep("DART", cv_r), rep("RF", cv_r) )
  ))

  return_table$V1=as.numeric(return_table$V1)
  names(return_table) = c("value", "model")

  return(return_table)
}




CV1




par(mfrow=c(1,1))
barplot(colMeans(m0$varprob_fp),  main="lsBART functional predictors prob.",
        space=0, names.arg  = paste0(var_names2), las=2)


par(mfrow=c(2,4))
for (i in 1:num_predictors) {
  barplot(colMeans(m0$varprob_times[[i]]), main=paste("", var_names2[i]),
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ) , space=0)
}


library(pROC)
rocobj <- roc(CV1[,6], CV1[,2])
best_threshold = coords(rocobj, "best")
best_threshold

AUC_cv_r
colMeans(AUC_cv_r)

ggplot(get_long_table(AUC_cv_r, cv_r), aes(x = model, value, colour=model))+
  geom_boxplot()+
  labs(y = "AUC Value",
       title = "ROC AUC BoxPlot")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))




colMeans(PR_cv_r)

ggplot(get_long_table(PR_cv_r, cv_r), aes(x = model, value, colour=model))+
  geom_boxplot()+
  labs(y = " PR AUC Value",
       title = "Precision-Recall AUC BoxPlot")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

ggplot(get_long_table(F1_cv_r, cv_r), aes(x = model, value, colour=model))+
  geom_boxplot()+
  labs(y = "F1 Score",
       title = "F1 BoxPlot")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16))



ggplot(get_long_table(Sens_cv_r, cv_r), aes(x = model, value, colour=model))+
  geom_boxplot()+
  labs(y = "Sensitivity",
       title = "Sensitivity BoxPlot")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16))



ggplot(get_long_table(Spec_cv_r, cv_r), aes(x = model, value, colour=model))+
  geom_boxplot()+
  labs(y = "Specificity Value",
       title = "Specificity BoxPlot")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16))






