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
  c1_s_data = cbind(c1_s_data, X_1 + rnorm(10,sd=0.5) )
  c2_s_data = cbind(c2_s_data, X_2 + rnorm(10,sd=0.5) )
  c3_s_data = cbind(c3_s_data, X_3 + rnorm(10,sd=0.5) )
  c4_s_data = cbind(c4_s_data, X_4 + rnorm(10,sd=0.5) )
  c5_s_data = cbind(c5_s_data, X_5 + rnorm(10,sd=0.5) )
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




d1_points = d1_points+rnorm( dim(d1_points)[2], sd=0.05 )
d2_points = d2_points+rnorm( dim(d2_points)[2], sd=0.05 )
d3_points = d3_points+rnorm( dim(d3_points)[2], sd=0.05)
d4_points = d4_points+rnorm( dim(d4_points)[2], sd=0.05 )
d5_points = d5_points+rnorm( dim(d5_points)[2], sd=0.05)


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

#the truth
true_y = friedmanFunc2(points1, ts1)
# Add Gaussian noise
true_y = true_y + rnorm(length(true_y))


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

flog_time_all=NULL
fbart_time_all=NULL
bart_time_all=NULL
dart_time_all=NULL
rf_time_all=NULL


for (num_p in c(1:1) ) {
  num_predictors = 5*num_p
  points = cbind(d1_points,
                 d2_points,
                 d3_points,
                 d4_points,
                 d5_points)

  n_point = length(grid)

  if(num_p>1){
    for(i in 1:(num_predictors-5) ){
      d_rand_points = matrix(rnorm(dim(d1_points)[1]*dim(d1_points)[2], sd=5),
                             ncol = dim(d1_points)[2],
                             nrow = dim(d1_points)[1])
      points = cbind(points, d_rand_points)
    }
  }



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
    print(paste("FLOG --", num_p))
  }
  else{    flog_test_all[num_p]= NaN
  print(paste("FLOG --", num_p)) }

  time1 =toc()
  flog_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ###########################################################
  ########################################################################
  #BART

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
  tic()
  m3 = wbart(x.train = as.matrix(train) , y.train = y_train, ntree=200,
             x.test =  as.matrix(test),  sparse = TRUE,  ndpost=100, nskip = 100)


  rmse_DART = sqrt(mean( (y_test- m3$yhat.test.mean)^2 ))
  dart_test_all[num_p]= rmse_DART
  print(paste("DART --", num_p))

  time1 =toc()
  dart_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]


  ########################################################################
  #Testing using Radnom Forest
  tic()

  m5 <- randomForest(y_train ~ ., data=train, importance=TRUE,
                     proximity=TRUE, ntree=200, mtry=100)

  p5 = predict(m5, test)

  rmse_RF = sqrt(mean( (y_test- p5)^2 ))
  rf_test_all[num_p]= rmse_RF

  print(paste("Random Forest --", num_p))

  time1 =toc()

  rf_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]

  ###########################################################
  #Functional Bart

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

  print(paste("lsBART --", num_p))

  time1 =toc()
  fbart_time_all[num_p] = str_split(time1$callback_msg, ' ')[[1]][1]

}


# Create a dataset
data <- data.frame(
  x = 5*c(1:10),
  y1 = flog_test_all,
  y2 = fbart_test_all,
  y3 = bart_test_all,
  y4 = dart_test_all,
  y5 = rf_test_all

)


# Create a dataset for plots
data2 <- data.frame(
  x = 5*c(1:10),
  y1 = as.numeric(flog_time_all),
  y2 = as.numeric(fbart_time_all),
  y3 = as.numeric(bart_time_all),
  y4 = as.numeric(dart_time_all),
  y5 = as.numeric(rf_time_all)

)



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
  labs(x = "Number of Predictor Functions ", y = "RMSE",
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
  labs(x = "Number of Predictor Functions ", y = "Seconds",
       title = "Time vs Number of Predictors")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=0.1))






###########################################
#See the functional Predictor Probabilities and Location Probabilities
par(mfrow=c(1,1))
barplot(colMeans(m0$varprob_fp),  main="Functional Predictor Probabilities",
        space=0, names.arg  = paste0("X_", 1:num_predictors))


barplot(colMeans(m0$varprob_times[[1]]), main="X_1 points lsBART",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[1]]), lwd=1, col="green", lty=2)

barplot(colMeans(m0$varprob_times[[2]]), main="X_2 points lsBART",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[2]]), lwd=1, col="green", lty=2)

barplot(colMeans(m0$varprob_times[[3]]), main="X_3 points lsBART",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[3]]), lwd=1, col="green", lty=2)

barplot(colMeans(m0$varprob_times[[4]]), main="X_4 points lsBART",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[4]]), lwd=1, col="green", lty=2)

barplot(colMeans(m0$varprob_times[[5]]), main="X_5 points lsBART",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[5]]), lwd=1, col="green", lty=2)


#############################
# See the true location in predictor functions
# and the location Probabilities
dev.off()
par(cex.lab = .01,  # Axis labels size
    cex.axis = 1.5, # Axis tick labels size
    cex.main = 1.5, # Main title size
    cex.sub = .01)  # Subtitle size

par(mfrow=c(2,1))
x


plot(Data2fd(c2_s_data, basis=bb ))
abline(v = (ts1[[2]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[2]]), main="X_2 Location Probability",
        names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 )) , space=0)
abline(v = (ts1[[2]]), lwd=1, col="green", lty=2)


plot(Data2fd(c3_s_data, basis=bb ))
abline(v = (ts1[[3]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[3]]), main="X_3 Location Probability",
        names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ))  , space=0)
abline(v = (ts1[[3]]), lwd=1, col="green", lty=2)


plot(Data2fd(c4_s_data, basis=bb ))
abline(v = (ts1[[4]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[4]]), main="X_4 Location Probability",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[4]]), lwd=1, col="green", lty=2)


plot(Data2fd(c5_s_data, basis=bb ))
abline(v = (ts1[[5]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[5]]), main="X_5 Location Probability",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[5]]), lwd=1, col="green", lty=2)

par(mfrow=c(1,1))

##########################################################################
# Getting Plots


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
  c1_s_data = cbind(c1_s_data, X_1 + rnorm(10,sd=0.5) )
  c2_s_data = cbind(c2_s_data, X_2 + rnorm(10,sd=0.5) )
  c3_s_data = cbind(c3_s_data, X_3 + rnorm(10,sd=0.5) )
  c4_s_data = cbind(c4_s_data, X_4 + rnorm(10,sd=0.5) )
  c5_s_data = cbind(c5_s_data, X_5 + rnorm(10,sd=0.5) )
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




d1_points = d1_points+rnorm( dim(d1_points)[2], sd=0.5 )
d2_points = d2_points+rnorm( dim(d2_points)[2], sd=0.5 )
d3_points = d3_points+rnorm( dim(d3_points)[2], sd=0.5)
d4_points = d4_points+rnorm( dim(d4_points)[2], sd=0.5 )
d5_points = d5_points+rnorm( dim(d5_points)[2], sd=0.5)


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

#the truth
true_y = friedmanFunc2(points1, ts1)
# Add Gaussian noise
true_y = true_y + rnorm(length(true_y))


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

flog_time_all=NULL
fbart_time_all=NULL
bart_time_all=NULL
dart_time_all=NULL
rf_time_all=NULL
par(mfrow=c(1,1))
plot(Data2fd(c1_s_data, basis=bb ), ylab = "X_1(t)")
abline(v = (ts1[[1]])/dim(d1_points)[2], lwd=1.5, col="green", lty=2)

plot(Data2fd(c2_s_data, basis=bb ), ylab = "X_2(t)")
abline(v = (ts1[[2]])/dim(d1_points)[2], lwd=1.5, col="green", lty=2)

plot(Data2fd(c3_s_data, basis=bb ), ylab = "X_3(t)")
abline(v = (ts1[[3]])/dim(d1_points)[2], lwd=1.5, col="green", lty=2)

plot(Data2fd(c4_s_data, basis=bb ), ylab = "X_4(t)")
abline(v = (ts1[[4]])/dim(d1_points)[2], lwd=1.5, col="green", lty=2)

plot(Data2fd(c5_s_data, basis=bb ), ylab = "X_5(t)")
abline(v = (ts1[[5]])/dim(d1_points)[2], lwd=1.5, col="green", lty=2)





