\name{LocationSmoothedBART-package}
\alias{LocationSmoothedBART-package}
\alias{LocationSmoothedBART}
\docType{package}
\title{
  Location Smoothed Bayesian Additive Regression Trees (lsBART)
}
\description{
  LocationSmoothedBART (lsBART) is an advanced R package that extends the Bayesian Additive Regression Trees (BART) model. It is specifically designed for scalar on function regression, making it particularly suited for medical image analysis in radiation treatment planning. This package focuses on identifying significant regions within predictor functions, ensuring accurate and reliable contour delineation.
}
\details{
  The lsBART package provides tools to apply a Bayesian tree-based model for scalar on function regression. It extends the traditional BART model by incorporating multi-level sparsity and spatial smoothing techniques. A key feature is the integration of Shapley values for detailed inference on binary response variables.

  The most important functions in this package include:

  \describe{
    \item{\code{w_lsBART}}{The core function to fit the location smoothed BART model to your data for continuous outcome variables.}
    \item{\code{p_lsBART}}{The probit location smoothed BART model for classification and binary outcomes.}
    \item{\code{splitTrees}}{Splits the tree string into individual trees for further processing.}
    \item{\code{getBARTTree}}{Extracts the structure of a specific tree from the BART model.}
    \item{\code{renameSplitVar}}{Renames the split variables in the tree data frame based on the dataset column names.}
    \item{\code{BART.unify}}{Unifies the BART model into a format suitable for treeshap calculations.}
    \item{\code{getBARTprediction}}{Traverses the BART tree and returns the prediction for a given observation.}
  }

  \subsection{w_lsBART}{
  Fits a Bayesian additive regression trees (BART) model for scalar-on-function regression with location smoothing.

  \emph{Input Arguments:}
  \describe{
    \item{\code{num_predictors}}{Number of functional predictors.}
    \item{\code{num_intvls}}{Number of points or intervals of each function.}
    \item{\code{x.train}}{Training data matrix of predictors.}
    \item{\code{y.train}}{Response vector for the training data.}
    \item{\code{x.test}}{Test data matrix of predictors (default is an empty matrix).}
    \item{\code{sparse}}{Logical flag indicating whether to use sparse BART.}
    \item{\code{theta}}{Sparsity parameter for BART.}
    \item{\code{omega}}{Sparsity parameter for BART.}
    \item{\code{a}}{Hyperparameter for the Dirichlet process prior.}
    \item{\code{b}}{Hyperparameter for the Dirichlet process prior.}
    \item{\code{augment}}{Logical flag indicating whether to augment the BART model.}
    \item{\code{rho}}{Proportion of predictors to be considered for splitting at each node (default is NULL).}
    \item{\code{dart_fp}}{Logical flag indicating whether to use DART for functional predictors.}
    \item{\code{theta_fp}}{Sparsity parameter for functional predictors.}
    \item{\code{omega_fp}}{Sparsity parameter for functional predictors.}
    \item{\code{a_fp}}{Hyperparameter for the Dirichlet process prior for functional predictors.}
    \item{\code{b_fp}}{Hyperparameter for the Dirichlet process prior for functional predictors.}
    \item{\code{augment_fp}}{Logical flag indicating whether to augment the DART model for functional predictors.}
    \item{\code{rho_fp}}{Proportion of functional predictors to be considered for splitting at each node (default is NULL).}
    \item{\code{dart_int}}{Logical flag indicating whether to use DART for interactions.}
    \item{\code{theta_int}}{Sparsity parameter for interactions.}
    \item{\code{omega_int}}{Sparsity parameter for interactions.}
    \item{\code{a_int}}{Hyperparameter for the Dirichlet process prior for interactions.}
    \item{\code{b_int}}{Hyperparameter for the Dirichlet process prior for interactions.}
    \item{\code{augment_int}}{Logical flag indicating whether to augment the DART model for interactions.}
    \item{\code{rho_int}}{Proportion of interaction terms to be considered for splitting at each node (default is NULL).}
    \item{\code{sigma_int}}{Standard deviation of the noise for the response variable.}
    \item{\code{xinfo}}{Matrix containing the cut points for the predictors.}
    \item{\code{usequants}}{Logical flag indicating whether to use quantiles for cut points.}
    \item{\code{cont}}{Logical flag indicating whether to treat all predictors as continuous.}
    \item{\code{rm.const}}{Logical flag indicating whether to remove constant predictors.}
    \item{\code{sigest}}{Initial estimate of the error standard deviation.}
    \item{\code{sigdf}}{Degrees of freedom for the inverse-chi-squared prior on the error variance.}
    \item{\code{sigquant}}{Quantile for the inverse-chi-squared prior on the error variance.}
    \item{\code{k}}{Prior parameter for the terminal node parameter.}
    \item{\code{power}}{Power parameter for the base tree prior.}
    \item{\code{base}}{Base parameter for the tree prior.}
    \item{\code{sigmaf}}{Standard deviation of the prior for the terminal node parameter.}
    \item{\code{lambda}}{Scale parameter for the inverse-chi-squared prior on the error variance.}
    \item{\code{fmean}}{Mean of the response variable.}
    \item{\code{w}}{Weights for the observations.}
    \item{\code{ntree}}{Number of trees in the ensemble.}
    \item{\code{numcut}}{Number of cut points for each predictor.}
    \item{\code{ndpost}}{Number of posterior samples to draw.}
    \item{\code{nskip}}{Number of MCMC iterations to skip before starting to collect posterior samples.}
    \item{\code{keepevery}}{Thinning interval for keeping MCMC samples.}
    \item{\code{nkeeptrain}}{Number of posterior samples to keep for the training data.}
    \item{\code{nkeeptest}}{Number of posterior samples to keep for the test data.}
    \item{\code{nkeeptestmean}}{Number of posterior mean samples to keep for the test data.}
    \item{\code{nkeeptreedraws}}{Number of posterior tree draws to keep.}
    \item{\code{printevery}}{Interval for printing progress information.}
    \item{\code{transposed}}{Logical flag indicating whether the input matrices are transposed.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{sigma}}{Posterior samples of the error standard deviation.}
    \item{\code{yhat.train.mean}}{Posterior mean of the predicted values for the training data.}
    \item{\code{yhat.train}}{Posterior samples of the predicted values for the training data.}
    \item{\code{yhat.test.mean}}{Posterior mean of the predicted values for the test data.}
    \item{\code{yhat.test}}{Posterior samples of the predicted values for the test data.}
    \item{\code{varcount}}{Posterior samples of the count of how often each predictor is used in the splits.}
    \item{\code{varcount.mean}}{Mean count of how often each predictor is used in the splits.}
    \item{\code{varprob}}{Posterior samples of the probability of each predictor being used in the splits.}
    \item{\code{varprob.mean}}{Mean probability of each predictor being used in the splits.}
    \item{\code{varcount_fp}}{Posterior samples of the count of how often each functional predictor is used in the splits.}
    \item{\code{varprob_fp}}{Posterior samples of the probability of each functional predictor being used in the splits.}
    \item{\code{varcount_int}}{Posterior samples of the count of how often each interaction term is used in the splits.}
    \item{\code{varprob_int}}{Posterior samples of the probability of each interaction term being used in the splits.}
    \item{\code{varprob_times}}{List of matrices showing the probability of each time point being used in the splits for each functional predictor.}
    \item{\code{treedraws}}{List containing the tree draws, including the cutpoints and trees as character strings.}
    \item{\code{proc.time}}{The processing time for the function call.}
    \item{\code{mu}}{The mean of the response variable used for centering.}
  }
  }

    \subsection{p_lsBART}{
  Fits a probit Bayesian additive regression trees (BART) model for binary outcomes with location smoothing. Similar input and outputs as w_lsBART so here we only display inputs and outputs unique to p_lsBART.

  \emph{Input Arguments:}
  \describe{
    \item{\code{binaryOffset}}{Offset for the binary response variable.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{prob.train}}{Posterior samples of the predicted probabilities for the training data.}
    \item{\code{prob.train.mean}}{Posterior mean of the predicted probabilities for the training data.}
    \item{\code{prob.test}}{Posterior samples of the predicted probabilities for the test data.}
    \item{\code{prob.test.mean}}{Posterior mean of the predicted probabilities for the test data.}
  }
  }

  \subsection{splitTrees}{
  Splits the tree string into individual trees for further processing.

  \emph{Input Arguments:}
  \describe{
    \item{\code{s}}{A string containing all trees from the BART model.}
    \item{\code{num_trees}}{An integer representing the number of trees.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{tree_vector}}{A vector of strings, each representing a single tree.}
  }
  }

  \subsection{getBARTTree}{
  Extracts the structure of a specific tree from the BART model.

  \emph{Input Arguments:}
  \describe{
    \item{\code{bartModel}}{A BART model object containing tree draws and cut points.}
    \item{\code{k}}{An integer representing the tree number to extract (default is 1).}
    \item{\code{labelVar}}{A logical flag indicating whether to label split variables (default is FALSE).}
    \item{\code{ntree}}{An integer representing the total number of trees.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{tree_df}}{A data frame representing the structure of the specified tree.}
  }
  }

  \subsection{renameSplitVar}{
  Renames the split variables in the tree data frame based on the dataset column names.

  \emph{Input Arguments:}
  \describe{
    \item{\code{tree_df}}{A data frame containing the tree structure with split variables.}
    \item{\code{dataset}}{A data frame containing the dataset used to train the model.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{tree_df}}{A data frame with split variable IDs replaced by their respective column names.}
  }
  }

  \subsection{BART.unify}{
  Unifies the BART model into a format suitable for treeshap calculations.

  \emph{Input Arguments:}
  \describe{
    \item{\code{rf_model}}{A BART model object.}
    \item{\code{data}}{A data frame containing the dataset.}
    \item{\code{numtree}}{An integer representing the number of trees.}
    \item{\code{ndpost}}{An integer representing the number of posterior draws.}
    \item{\code{ndpost_draw}}{An integer representing the specific posterior draw to unify.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{ret}}{A list containing a unified model and data, formatted for use with the TreeShap package.}
  }
  }

  \subsection{getBARTprediction}{
  Traverses the BART tree and returns the prediction for a given observation.

  \emph{Input Arguments:}
  \describe{
    \item{\code{tree_df}}{A data frame representing the structure of a tree.}
    \item{\code{observation}}{A named vector or data frame row containing the observation to predict.}
  }

  \emph{Outputs:}
  \describe{
    \item{\code{prediction}}{A numeric value representing the predicted response for the given observation.}
  }
  }


  For detailed examples and tutorials, please refer to the vignettes and the documentation available on the package's GitHub page.
}
\author{
  Maintainer: Zachary Wooten <ztw5@rice.edu>
}
\references{
  Wooten, Z. et al. (2024). "Location smoothed BART: A method for interpretable and robust quality assurance of organ contours in radiotherapy treatment planning" (Under Review)

  This work is built on the BART R package by Sparapani, Spanbauer, and McCulloch. For additional background information and related methodologies, refer to the BART R package documentation and associated research papers.

  Sparapani R, Spanbauer C, McCulloch R (2021). “Nonparametric Machine Learning and Efficient Computation with Bayesian
  Additive Regression Trees: The BART R Package.” _Journal of Statistical Software_, *97*(1), 1-66.
  doi:10.18637/jss.v097.i01 <https://doi.org/10.18637/jss.v097.i01>.

  For the implementing the TreeShap algorithm within the BART framework we built on the TreeShap package by Komisarczyk, Kozminski, Maksymiuk, and Biecek. For additional background information and related methodologies, refer to the TreeShap package documentation and associated research papers.

  Komisarczyk K, Kozminski P, Maksymiuk S, Biecek P (2023). _treeshap: Fast SHAP values computation for tree ensemble
  models_. R package version 0.1.1.9001, <https://github.com/ModelOriented/treeshap>.


}
\keyword{ package }
\seealso{
  \code{\link{BART}}, \code{\link{fda}}, \code{\link{tidyverse}}
}
\examples{
  \dontrun{
     # Example of fitting the lsBART model to simulated data

     # Load necessary libraries
     library(fda)
     library(tidyverse)
     library(LocationSmoothedBART)

     # Generate simulated data
     set.seed(34)
     X_1 <- rnorm(10, sd=1)
     X_2 <- rnorm(10, sd=1)
     c1_s_data <- replicate(2000, X_1 + rnorm(10, sd=0.5))
     c2_s_data <- replicate(2000, X_2 + rnorm(10, sd=0.5))

     bb <- create.bspline.basis(rangeval = c(0, 1), norder=8)
     d1 <- Data2fd(c1_s_data, basis=bb)
     d1_points <- t(d1$coefs) %*% t(eval.basis(seq(0, 1, by=0.01), bb))
     d2 <- Data2fd(c2_s_data, basis=bb)
     d2_points <- t(d2$coefs) %*% t(eval.basis(seq(0, 1, by=0.01), bb))

     # Combine points data
     points <- cbind(d1_points, d2_points)
     df <- as_tibble(data.frame(points))

     # Split into train and test sets
     set.seed(1)
     sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.5,0.5))
     train <- df[sample, ]
     test <- df[!sample, ]

     # Generate response variable
     t1 <- 1:10
     t2 <- 21:30
     true_y <- 10 * rowSums(sin(pi * d1_points[, t1])) +
               20 * rowSums((d2_points[, t2] - 0.5)^2) + rnorm(nrow(d1_points))
     y_train <- true_y[sample]
     y_test <- true_y[!sample]

     # Fit the model
     m0 <- w_lsBART(num_predictors = 2, num_intvls = ncol(d1_points),
                    x.train = as.matrix(train), y.train = y_train,
                    x.test = as.matrix(test), sparse = TRUE,
                    ndpost = 100, nskip = 100, ntree = 200,
                    dart_fp = TRUE, dart_int = TRUE)


     # Calculate RMSE
     predictions <- colMeans(m0$yhat.test)
     rmse <- sqrt(mean((y_test - predictions)^2))
     print(paste("RMSE:", rmse))

     # Plot predictions vs y_test
      plot(y_test, predictions, main="Predictions vs Actual",
        xlab="Actual y_test", ylab="Predictions", pch=19, col="blue")
        abline(0, 1, col="red", lwd=2)


     # Plot functional predictor probabilities
     barplot(colMeans(m0$varprob_fp), main="Functional Predictor Probabilities",
             space=0, names.arg = paste0("X_", 1:2))

     # Plot location probabilities for each predictor
     par(mfrow=c(2,1), mar=c(4, 4, 2, 1))
     plot(Data2fd(c1_s_data, basis=bb))
     abline(v = (t1) / ncol(d1_points), lwd=1, col="green", lty=2)
     barplot(colMeans(m0$varprob_times[[1]]), main="X_1 Location Probability",
             names.arg = as.character(1:ncol(d1_points)), space=0)
     abline(v = t1, lwd=1, col="green", lty=2)

     plot(Data2fd(c2_s_data, basis=bb))
     abline(v = (t2) / ncol(d2_points), lwd=1, col="green", lty=2)
     barplot(colMeans(m0$varprob_times[[2]]), main="X_2 Location Probability",
             names.arg = as.character(1:ncol(d2_points)), space=0)
     abline(v = t2, lwd=1, col="green", lty=2)


  }
}
