# LocationSmoothedBART
Code to run lsBART

## Overview

LocationSmoothedBART (lsBART) is an innovative R package that extends the Bayesian Additive Regression Trees (BART) model. It is specifically designed for scalar on function regression and is particularly suited for medical image analysis in radiation treatment planning. This package emphasizes identifying significant regions within predictor functions, ensuring accurate and reliable contour delineation.

## Features

- Extends the traditional BART model for more precise medical image analysis.
- Utilizes shape features for contour quality assurance, independent of the imaging platform.
- Incorporates Shapley values for detailed error analysis in contour delineation.

## Installation

``` r
# To install the latest version from GitHub:
devtools::install_github("wootz101/LocationSmoothedBART")
```

## Usage

Load the package and apply lsBART to your data:

``` r
library(fda)
library(tidyverse)
library(LocationSmoothedBART)



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

partialfriedmanFunc3 <- function(points, ts){
  
  res = 10*rowSums(sin(pi* points[[1]][, ts[[1]]] ) )+
    20*rowSums((points[[2]][, ts[[2]] ] -0.5)^2)
  return(res)
}


points1 = list(d1_points, d2_points)
ts1 = list(t1,t2)

#the truth
true_y = partialfriedmanFunc3(points1, ts1)
# Add Gaussian noise
true_y = true_y + rnorm(length(true_y))


num_predictors = 5
points = cbind(d1_points,
               d2_points,
               d3_points,
               d4_points,
               d5_points)

#create Train and Test sets
df = data.frame(points)
df = as_tibble(df)

#use 50% of dataset as training set and 30% as test set
set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.5,0.5))
train  <- df[sample, ]
test   <- df[!sample, ]

y_train = true_y[sample]
y_test = true_y[!sample]



## Little snippet example
m0 = w_lsBART(num_predictors = num_predictors,
              num_intvls =num_intvls,
              x.train = as.matrix(train),
              y.train = y_train,
              x.test =  as.matrix(test),
              sparse = TRUE,
              ndpost = 100,
              nskip = 100,
              ntree = 200,
              dart_fp = TRUE,
              dart_int = TRUE)




###########################################
#See the functional Predictor Probabilities and Location Probabilities
par(mfrow=c(1,1))
barplot(colMeans(m0$varprob_fp),  main="Functional Predictor Probabilities",
        space=0, names.arg  = paste0("X_", 1:num_predictors))


#See Location probabilities
par(mfrow=c(2,1))
plot(Data2fd(c1_s_data, basis=bb ))
abline(v = (ts1[[1]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[1]]), main="X_1 Location Probability",
        names.arg  = as.character(1:num_intvls) , space=0)
abline(v = (ts1[[1]]), lwd=1, col="green", lty=2)



plot(Data2fd(c2_s_data, basis=bb ))
abline(v = (ts1[[2]])/dim(d1_points)[2], lwd=1, col="green", lty=2)
barplot(colMeans(m0$varprob_times[[2]]), main="X_2 Location Probability",
        names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 )) , space=0)
abline(v = (ts1[[2]]), lwd=1, col="green", lty=2)

```

## Documentation

For detailed documentation and vignettes on how to use LocationSmoothedBART, visit our [documentation page](https://[your-github-username].github.io/LocationSmoothedBART/).

## Contributing

Contributions to the LocationSmoothedBART package are welcome. Please refer to our contributing guidelines for more information.

## Citation

If you use LocationSmoothedBART in your research, please cite:

```
-- Not Published Yet --
```

## Contact

Zachary Wooten
ztw5@rice.edu
(https://wootz101.github.io/Website2/)
