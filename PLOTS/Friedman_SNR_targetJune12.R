####################################################

# Changing SNR

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
library(FADPclust)
library(fda)
library(tidyverse)
library(randomForest)
library(BART)
library(tictoc)
library(LocationSmoothedBART)
library(treeshap)
library(glmnet)
library(nnet)
library(caret)

#Generate a 1000 for each class
set.seed(34)

X_1 = rnorm(10, sd=1)
X_2 = rnorm(10, sd=1)
X_3 = rnorm(10, sd=1)
X_4 = rnorm(10, sd=1)
X_5 = rnorm(10, sd=1)


c1_s_data = NULL
c2_s_data = NULL
c3_s_data = NULL
c4_s_data = NULL
c5_s_data = NULL

for(i in 1:2000){
  c1_s_data = cbind(c1_s_data, X_1 + rnorm(10,sd=1) )
  c2_s_data = cbind(c2_s_data, X_2 + rnorm(10,sd=1) )
  c3_s_data = cbind(c3_s_data, X_3 + rnorm(10,sd=1) )
  c4_s_data = cbind(c4_s_data, X_4 + rnorm(10,sd=1) )
  c5_s_data = cbind(c5_s_data, X_5 + rnorm(10,sd=1) )
}
sim_data_1 = cbind(c1_s_data, c2_s_data, c3_s_data, c4_s_data, c5_s_data)



####################
# Simulated Data
bb = create.bspline.basis(rangeval = c(0,1), norder=8)

############################################
# CREATE THE POINTS

d1 = Data2fd(c1_s_data, basis=bb)
d1_points = t(d1$coefs)%*%t(eval.basis(seq(by =0.01,from=0,to=1), bb ))

d2 = Data2fd(c2_s_data, basis=bb)
d2_points = t(d2$coefs)%*%t(eval.basis(seq(by =0.01,from=0,to=1), bb ))

d3 = Data2fd(c3_s_data, basis=bb)
d3_points = t(d3$coefs)%*%t(eval.basis(seq(by =0.01,from=0,to=1), bb ))

d4 = Data2fd(c4_s_data, basis=bb)
d4_points = t(d4$coefs)%*%t(eval.basis(seq(by =0.01,from=0,to=1), bb ))

d5 = Data2fd(c5_s_data, basis=bb)
d5_points = t(d5$coefs)%*%t(eval.basis(seq(by =0.01,from=0,to=1), bb ))




d1_points = d1_points+rnorm( dim(d1_points)[2], sd=1 )
d2_points = d2_points+rnorm( dim(d2_points)[2], sd=1 )
d3_points = d3_points+rnorm( dim(d3_points)[2], sd=1)
d4_points = d4_points+rnorm( dim(d4_points)[2], sd=1 )
d5_points = d5_points+rnorm( dim(d5_points)[2], sd=1)


num_intvls = dim(d1_points)[2]


# CREATE THE Y VARIABLE

#Friedman MARS paper sim data

#The Time Points we want to infer
t1 = 1:10
t2 = 21:30
t3 = 41:50
t4 = 61:70
t5 = 81:90

friedmanFunc <- function(x){

  res = 10*sin(pi*x[1]*x[2]) + 20*(x[3]-0.5)^2 + 10*x[4] + 5*x[5]
  return(res)
}


friedmanFunc2 <- function(points, ts){

  res = 10*rowSums(sin(pi* points[[1]][, ts[[1]]] * points[[2]][, ts[[2]] ]) ) +
    20*rowSums((points[[3]][, ts[[3]]  ] -0.5)^2) +
    10*rowSums((points[[4]][, ts[[4]]  ]))+
    5*rowSums((points[[5]][, ts[[5]]  ]))
  return(res)
}


points1 = list(d1_points, d2_points, d3_points, d4_points, d5_points)
ts1 = list(t1,t2,t3,t4,t5)



####################################################################################################
####################################################################################################
tempbasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)
grid = seq(.01, from=0, to=1)
numPoinst=length(grid)+2
bb <- create.bspline.basis(rangeval = c(0,1), breaks = grid, nbasis = numPoinst)


flog_test_all=NULL
fbart_test_all=NULL
bart_test_all=NULL
dart_test_all=NULL
rf_test_all=NULL
lasso_test_all=NULL
ridge_test_all=NULL
nn_test_all=NULL


flog_time_all=NULL
fbart_time_all=NULL
bart_time_all=NULL
dart_time_all=NULL
rf_time_all=NULL
lasso_time_all=NULL
ridge_time_all=NULL
nn_time_all=NULL


flog_SNR=NULL
fbart_SNR=NULL
bart_SNR=NULL
dart_SNR=NULL
rf_SNR=NULL

num_predictors = 15
points = cbind(d1_points,
                 d2_points,
                 d3_points,
                 d4_points,
                 d5_points)

n_point = length(grid)

for(i in 1:(num_predictors-5) ){
      d_rand_points = matrix(rnorm(dim(d1_points)[1]*dim(d1_points)[2], sd=5),
                             ncol = dim(d1_points)[2],
                             nrow = dim(d1_points)[1])
      points = cbind(points, d_rand_points)
}




#snr_target = c( 1, 5, 10, 100)

#snr_target = c( 1, 2.5, 5, 7.5, 10, 12)

#snr_target = c( 1, 10, 100, 500, 1000)

#snr_target = c(0.5, 1, 2.5, 5, 10, 25, 50, 100)

#snr_target = c(0.5, 1, 5, 10, 50, 100)

snr_target = c( 1, 5, 10, 50, 100)

for (num_p in c(1:length(snr_target)) ) {
set.seed(1)
  #create new noisy predictor

  #the truth
  true_signal = friedmanFunc2(points1, ts1)
  # Add Gaussian noise

  signal_variance <- var(true_signal)
  noise_variance <- signal_variance / snr_target[num_p]

  true_y = true_signal + rnorm(length(true_signal),
                               sd = sqrt(noise_variance))




  #create Train and Test sets
  df= data.frame(points)
  df = as_tibble(df)

  #use 50% of dataset as training set and 30% as test set
  set.seed(num_p)
  sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.5,0.5))
  train  <- df[sample, ]
  test   <- df[!sample, ]

  y_train = true_y[sample]
  y_test = true_y[!sample]

  ###################################################################
  set.seed(314)
  tic()
  lasso_model <- glmnet(as.matrix(train), y_train, alpha = 1)
  lasso_model$beta

  cv_lasso <- cv.glmnet(as.matrix(train), y_train, alpha = 1)

  # Best lambda value
  best_lambda <- cv_lasso$lambda.min
  final_lasso_model <- glmnet(as.matrix(train), y_train, alpha = 1, lambda = best_lambda)


  predictions_lasso <- predict(final_lasso_model, newx = as.matrix(test))

  rmse_lasso = sqrt(mean( (y_test- predictions_lasso)^2 ))
  lasso_test_all[num_p]= rmse_lasso

  print(paste("LASSO --", num_p))

  time1 =toc()
  lasso_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ###################################################################


  ###################################################################
  set.seed(314)
  tic()
  ridge_model <- glmnet(as.matrix(train), y_train, alpha = 0)
  ridge_model$beta

  cv_ridge <- cv.glmnet(as.matrix(train), y_train, alpha = 0)

  # Best lambda value
  best_lambda <- cv_ridge$lambda.min
  final_ridge_model <- glmnet(as.matrix(train), y_train, alpha = 0, lambda = best_lambda)


  predictions_ridge <- predict(final_ridge_model, newx = as.matrix(test))

  rmse_ridge = sqrt(mean( (y_test- predictions_ridge)^2 ))
  ridge_test_all[num_p]= rmse_ridge

  print(paste("Ridge --", num_p))

  time1 =toc()
  ridge_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ###################################################################

  ###################################################################
  tic()
  # Normalize the predictor variables
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }



  # Normalize the data
  train_norm <- as.data.frame(lapply(train, normalize))

  # Perform PCA
  pca <- preProcess(train_norm, method = "pca", pcaComp = 50)  # Reduce to 50 components
  train_pca <- predict(pca, train_norm)

  # Combine the PCA components with the response variable
  train_pca <- data.frame(train_pca, y_train)

  # Train the neural network model
  set.seed(314)  # For reproducibility
  nn_model <- nnet(y_train ~ ., data = train_pca, size = 5, maxit = 1500,
                   linout = TRUE, trace = FALSE)

  # Make predictions
  test_norm <- as.data.frame(lapply(test, normalize))
  test_pca <- predict(pca, test_norm)
  predictions_nn <- predict(nn_model, test_pca)


  rmse_nn <- sqrt(mean((predictions_nn - y_test)^2))
  rmse_nn
  nn_test_all[num_p]= rmse_nn

  print(paste("NN --", num_p))

  time1 =toc()
  nn_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ###################################################################
  ###################################################################





  tic()

  if(num_p<100){
    #ALL SET UP

    tempbasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)
    grid <- seq(0,1,1/(length(train[1,1:n_point])-1))

    templist      = vector("list", num_predictors+1)
    templist[[1]] = rep(1,  nrow(train))

    for (i in 0:(num_predictors-1) ) {
      tempSmooth= smooth.basis(grid, t(train[,(n_point*i+1):(n_point*i+n_point)]), tempbasis)
      templist[[i+2]] = tempSmooth$fd
      print(i+2)

    }

    #t(train[,(n_point*0+1):(n_point*0+n_point)])

    conbasis = create.constant.basis(c(0,1))
    betabasis = create.bspline.basis(rangeval = c(0,1), nbasis=4)

    betalist = vector("list",num_predictors+1)
    betalist[[1]] = conbasis

    for (i in 2:(num_predictors+1) ) {
      betalist[[i]] = betabasis

    }

    #Logistic Regression

    fRegressList2 = fRegress(y_train , templist, betalist)

    # Test set

    #Create fucntional data for test set
    testlist      = vector("list", num_predictors+1)
    testlist[[1]] = rep(1,  nrow(test))
    grid <- seq(0,1,1/(length(test[1,1:n_point])-1))

    for (i in 0:(num_predictors-1)) {
      testSmooth= smooth.basis(grid, t(test[  ,(n_point*i+1):(n_point*i+n_point) ]), tempbasis)
      testlist[[i+2]] = testSmooth$fd

    }

    blist <- fRegressList2$betaestlist

    #Use the coefficients from the training set to predict from unseeen
    #test set. Call it b_coefs

    b_coefs <- NULL
    b_coefs <- rbind(b_coefs, rep(blist[[1]]$fd$coefs, length(y_test)))
    for (i in 2:(num_predictors+1)) {
      b_coefs <- rbind(b_coefs, inprod(blist[[i]]$fd, testlist[[i]]))
    }

    #Answer Prediction
    yhat_test_fReg = colSums(b_coefs)

    rmse_FLOG = sqrt(mean( (y_test-yhat_test_fReg)^2 ))
    flog_test_all[num_p]= rmse_FLOG

    flog_SNR[num_p]= var(y_test)/var((y_test-yhat_test_fReg))

    print(paste("FLOG --", num_p))
  }
  else{    flog_test_all[num_p]= NaN
  print(paste("FLOG --", num_p)) }

  time1 =toc()
  flog_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ###########################################################
  ########################################################################
  #BART
  set.seed(314)

  tic()
  m2 = wbart(x.train = as.matrix(train) , y.train = y_train, ntree=200,
             x.test =  as.matrix(test),  sparse = FALSE,  ndpost=100, nskip = 100)


  rmse_BART = sqrt(mean( (y_test- m2$yhat.test.mean)^2 ))
  bart_test_all[num_p]= rmse_BART

  print(paste("BART --", num_p))

  time1 =toc()

  bart_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]

  ########################################################################
  #DART
  set.seed(314)
  tic()
  m3 = wbart(x.train = as.matrix(train) , y.train = y_train, ntree=200,
             x.test =  as.matrix(test),  sparse = TRUE,  ndpost=100, nskip = 100)


  rmse_DART = sqrt(mean( (y_test- m3$yhat.test.mean)^2 ))
  dart_test_all[num_p]= rmse_DART

  dart_SNR[num_p]= var(y_test)/var((y_test-m3$yhat.test.mean))

  print(paste("DART --", num_p))

  time1 =toc()
  dart_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ########################################################################
  #Testing using Radnom Forest
  set.seed(314)
  tic()
  m5 <- randomForest(y_train ~ ., data=train, importance=TRUE,
                     proximity=TRUE, ntree=200, mtry=100)

  p5 = predict(m5, test)

  rmse_RF = sqrt(mean( (y_test- p5)^2 ))
  rf_test_all[num_p]= rmse_RF
  rf_SNR[num_p]= var(y_test)/var((y_test-p5))


  print(paste("Random Forest --", num_p))

  time1 =toc()

  rf_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]

  ###########################################################
  #Functional Bart
  set.seed(314)
  tic()
  m0 = w_lsBART(num_predictors = num_predictors,
                num_intvls =num_intvls,
                x.train = as.matrix(train),
                y.train = y_train,
                x.test =  as.matrix(test),
                sparse = TRUE,
                ndpost = 100,
                nskip = 100,
                ntree=200,
                dart_fp = TRUE,
                dart_int = TRUE)


  ###########################################
  #### RMSE (SEE HOW GOOD IT IS)
  rmse_FBART = sqrt(mean( (y_test-m0$yhat.test.mean)^2 ))
  fbart_test_all[num_p]= rmse_FBART

  fbart_SNR[num_p]= var(y_test)/var((y_test-m0$yhat.test.mean))

  print(paste("lsBART --", num_p))

  time1 =toc()
  fbart_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]

}


# Create a dataset
data <- data.frame(
  x = snr_target,
  y1 = flog_test_all,
  y2 = fbart_test_all,
  y3 = bart_test_all,
  y4 = dart_test_all,
  y5 = rf_test_all,
  y6 = lasso_test_all,
  y7 = ridge_test_all,
  y8 = nn_test_all

)




ggplot(data, aes(x = x)) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  geom_line(aes(y = y6, color = "LASSO"), size = 1.5) +
  geom_point(aes(y = y6, color = "LASSO"), size = 3) +
  geom_line(aes(y = y7, color = "Ridge"), size = 1.5) +
  geom_point(aes(y = y7, color = "Ridge"), size = 3) +
  geom_line(aes(y = y8, color = "NN"), size = 1.5) +
  geom_point(aes(y = y8, color = "NN"), size = 3) +
  labs(x = "SNR", y = "RMSE",
       title = "Test Set Root Mean Squared Error") +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1)
  )




ggplot(data, aes(x = x)) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  geom_line(aes(y = y6, color = "LASSO"), size = 1.5) +
  geom_point(aes(y = y6, color = "LASSO"), size = 3) +
  geom_line(aes(y = y7, color = "Ridge"), size = 1.5) +
  geom_point(aes(y = y7, color = "Ridge"), size = 3) +
  geom_line(aes(y = y8, color = "NN"), size = 1.5) +
  geom_point(aes(y = y8, color = "NN"), size = 3) +
 # scale_x_log10() +  # Add this line to set the x-axis to log scale
  scale_x_continuous(trans='log10') +
  labs(x = "SNR", y = "RMSE",
       title = "Test Set Root Mean Squared Error") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 0.1)
  )





names(data) = c("x", "F. Reg.","lsBART", "BART", "DART", "RF", "LASSO", "Ridge", "NN" )
data
colMeans(data)










































ggplot(data, aes(x = log(x))) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  geom_line(aes(y = y6, color = "LASSO"), size = 1.5) +
  geom_point(aes(y = y6, color = "LASSO"), size = 3) +
  geom_line(aes(y = y7, color = "Ridge"), size = 1.5) +
  geom_point(aes(y = y7, color = "Ridge"), size = 3) +
  geom_line(aes(y = y8, color = "NN"), size = 1.5) +
  geom_point(aes(y = y8, color = "NN"), size = 3) +
  labs(x = "log scale SNR", y = "RMSE",
       title = "Test Set Root Mean Squared Error") +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1)
  )





























# Create a dataset for plots
data2 <- data.frame(
  x = noise_sig,
  y1 = as.numeric(flog_time_all),
  y2 = as.numeric(fbart_time_all),
  y3 = as.numeric(bart_time_all),
  y4 = as.numeric(dart_time_all),
  y5 = as.numeric(rf_time_all)

)

# Create a dataset for plots
data_SNR <- data.frame(
  x = noise_sig,
  y1 = as.numeric(flog_SNR),
  y2 = as.numeric(fbart_SNR),
  y3 = as.numeric(bart_SNR),
  y4 = as.numeric(dart_SNR),
  y5 = as.numeric(rf_SNR)

)

# Create a dataset for plots
data_SNR_RMSE <- data.frame(
  x1 = as.numeric(flog_SNR),
  x2 = as.numeric(fbart_SNR),
  x3 = as.numeric(bart_SNR),
  x4 = as.numeric(dart_SNR),
  x5 = as.numeric(rf_SNR),
  y1 = flog_test_all,
  y2 = fbart_test_all,
  y3 = bart_test_all,
  y4 = dart_test_all,
  y5 = rf_test_all
)


data_SNR_RMSE





# Transform the data to long format for both x and y values
data_long <- data_SNR_RMSE %>%
  pivot_longer(cols = starts_with("x"), names_to = "model_x", values_to = "x") %>%
  pivot_longer(cols = starts_with("y"), names_to = "model_y", values_to = "y") %>%
  filter(substr(model_x, 2, 2) == substr(model_y, 2, 2)) %>%
  select(x, y, model = model_x)

data_long <- data_SNR_RMSE %>%
  pivot_longer(cols = starts_with("x"), names_to = "model_x", values_to = "x") %>%
  pivot_longer(cols = starts_with("y"), names_to = "model_y", values_to = "y") %>%
  filter(substr(model_x, 2, 2) == substr(model_y, 2, 2)) %>%
  select(x, y, model = model_x) %>%
  mutate(model = recode(model,
                        "x1" = "Func. Reg.",
                        "x2" = "lsBART",
                        "x3" = "BART",
                        "x4" = "DART",
                        "x5" = "RF"))




ggplot(data_long, aes(x = x, y = y, color = model, group = model)) +
  geom_line(size = 1.5, alpha = 0.75) +
  geom_point(size = 3) +
  labs(title = "Line plot of SNR vs RMSE for different models",
       x = "SNR",
       y = "RMSE") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 0.1))








# Create a line plot with points and different colors for each line
ggplot(data, aes(x = x)) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  labs(x = "Sigma", y = "RMSE",
       title = "Test Set Root Mean Squared Error") +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1)
  )

# Create a line plot with points and different colors for each line
ggplot(data_SNR, aes(x = x)) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  labs(x = "Sigma", y = "SNR",
       title = "SNR") +
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1)
  )


# Create a line plot with points and different colors for each line
ggplot(data2, aes(x = x)) +
  geom_line(aes(y = y1, color = "F. Reg."), size = 1.5) +
  geom_point(aes(y = y1, color = "F. Reg."), size = 3) +
  geom_line(aes(y = y2, color = "lsBART"), size = 1.5) +
  geom_point(aes(y = y2, color = "lsBART"), size = 3) +
  geom_line(aes(y = y3, color = "BART"), size = 1.5) +
  geom_point(aes(y = y3, color = "BART"), size = 3) +
  geom_line(aes(y = y4, color = "DART"), size = 1.5) +
  geom_point(aes(y = y4, color = "DART"), size = 3) +
  geom_line(aes(y = y5, color = "RF"), size = 1.5) +
  geom_point(aes(y = y5, color = "RF"), size = 3) +
  labs(x = "Sigma ", y = "Seconds",
       title = "Time vs Number of Predictors")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1))

