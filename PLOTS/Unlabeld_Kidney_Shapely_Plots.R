# Shapely Value plots for Unlabeled Kidney feature importance

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))
library(FADPclust)
library(tictoc)
library(randomForest)
library(fda)
library(precrec)
library(tidyverse)
library(ROCit)
library(randomForest)
library(LocationSmoothedBART)
library(treeshap)

source("TreeShap_Source.R")

####################################################################################################################################
####################################################################################################################################

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


#Flag from 36 kidney contours that are unlabeled
var_names[c(1:10, 12, 13, 15, 17, 18)]
name_vec = c(1:10, 12, 13, 15, 17, 18)
var_names2 = var_names[name_vec]
num_predictors = length(var_names2)

unlabeledk <- read_csv("/Users/Wootz/Desktop/ORGAN DATA/Unlabeled_Kidney_2Dfunc_Sept24.csv")


set_unlab <- split(unlabeledk, unlabeledk$name)

factor_unlab = factor(unlabeledk$name)
#Gives the plan name
levels(factor_unlab)

#Only want some of the functional predictors to keep

grid = seq(.01, from=0, to=1)
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

num_intvls = dim(data_set)[2]/num_predictors
data_set

#create Train and Test sets
df = as_tibble(data_set)

data_set_unlab <- NULL

num_predictors = length(var_names2)



for (j in 1:length(set_unlab)) {
  data_temp <- NULL

  for (i in name_vec ) {
    f <- Data2fd(set_unlab[[j]][[i+2]], basisobj = bb)
    data_temp <- c(data_temp, f$coefs[,1])
    #print(paste(j, '---', i))
  }
  data_set_unlab <- rbind(data_set_unlab, data_temp)
}




data_set_unlab
########################################################
# set up data
par(mfrow=c(1,1))
df_unlab = as_tibble(data_set_unlab)

plot( as.numeric(df_unlab[19, 1:num_intvls]) )

outcomes = c(rep(0,260), rep(1, 52))
outcomes_class = c(rep("Acceptable",260), rep("Unacceptable", 52))
true_y = outcomes

# Testing out FBART



num_intvls = dim(data_set)[2]/num_predictors
dim(train)

train  <- as.matrix(df)
test   <- as.matrix(df_unlab)
y_train = true_y


#best case
binaryOffset = NULL
rho_int = 1
rho_fp = 20
sigma_int = 2

nskip = 100
kk = 2
power = 2
base = 0.5

ntree = 50
ndpost = 100


set.seed(1)
m00 = p_lsBART(num_predictors = num_predictors,
                                    num_intvls =num_intvls,
                                    x.train = train ,
                                    y.train = y_train,
                                    x.test =  test,
                                    rho_int=rho_int,
                                    rho_fp = rho_fp,
                                    sigma_int=sigma_int,
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

# Accuracy of hits

# Best threshold for training data
library(pROC)
rocobj <- roc(outcomes, m00$prob.train.mean)
rocobj
best_threshold = coords(rocobj, "best")
best_threshold


table(m00$prob.train.mean>0.5,
      outcomes)

table(m00$prob.train.mean>best_threshold$threshold[1],
      outcomes)

#Traditional threshold value
confusionMatrix( as.factor( as.integer(m00$prob.train.mean>0.5) ),
                 as.factor(outcomes), mode = "everything", positive="1")

#Youden's index
confusionMatrix( as.factor( as.integer(m00$prob.train.mean>best_threshold$threshold[1]) ),
                 as.factor(outcomes), mode = "everything", positive="1")


roc.curve( m00$prob.train.mean[outcomes==1],  m00$prob.train.mean[outcomes==0])
pr.curve(m00$prob.train.mean[outcomes==1],  m00$prob.train.mean[outcomes==0] )



##################################################################################
# Interpret results
prediction_results = tibble(names = levels(factor_unlab),
                            prob = m00$prob.test.mean)
#get credible intervals
library(bayestestR)



upper_ci=NULL
lower_ci=NULL
stdev=NULL
for (i  in 1:dim(m00$prob.test)[2]) {
  ci1 = ci(m00$prob.test[,i], method = "HDI", ci=0.95)
  upper_ci[i] = ci1$CI_high
  lower_ci[i] = ci1$CI_low
  stdev[i] = sd(m00$prob.test[,i])
}

df_test_probs = tibble(index = 1:dim(m00$prob.test)[2],
                       index_order = order(m00$prob.test.mean),
                       mean_prob_sort = sort(m00$prob.test.mean),
                       mean_prob = m00$prob.test.mean,
                       upper_ci = upper_ci,
                       lower_ci =lower_ci,
                       stdev = stdev)



par(mfrow=c(1,1))
barplot(colMeans(m00$varprob_fp),  main="lsBART functional predictors prob.",
        space=0, names.arg  = paste0(var_names2), las=2)


par(mfrow=c(2,4))
for (i in 1:num_predictors) {
  barplot(colMeans(m00$varprob_times[[i]]), main=paste("", var_names2[i]),
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ) , space=0)
}
par(mfrow=c(1,1))



################################################################################################
#Understanding prediction results


prediction_results$names[ df_test_probs$mean_prob>0.5]
which(df_test_probs$mean_prob>0.5)

prediction_results$names[ df_test_probs$mean_prob+df_test_probs$stdev>0.5]
which(df_test_probs$mean_prob+df_test_probs$stdev>0.5)

obs_ids = which(df_test_probs$mean_prob+df_test_probs$stdev>0.5)
obs_ids_names = prediction_results$names[ df_test_probs$mean_prob+df_test_probs$stdev>0.5]




#############
# Now we have the truth on errors in the plans
N=0
Y=1

unlab_truth = c(N,
                Y,
                N,
                Y,
                Y,
                N,
                Y,
                Y,
                N,
                Y,
                N,
                Y,
                N,
                N,
                N,
                N,
                N,
                Y,
                N,
                N,
                N,
                Y,
                Y,
                N,
                N,
                Y,
                N,
                Y,
                N,
                N,
                N,
                Y,
                N,
                N,
                Y,
                Y )
table(unlab_truth)



table(df_test_probs$mean_prob+df_test_probs$stdev>0.5,
      unlab_truth)

confusionMatrix( as.factor( as.integer(df_test_probs$mean_prob+df_test_probs$stdev>0.5) ),
                 as.factor(unlab_truth), mode = "everything", positive="1")



ggplot(df_test_probs, aes(x=index, y=mean_prob, color= factor(unlab_truth) ))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_prob - stdev, ymax=mean_prob+stdev), width=1,
                position=position_dodge(0.05))+
  geom_hline(yintercept=c(0.5, best_threshold$threshold[1]),
             linetype='dashed', color=c( 'grey', 'black'))+
  ylim(-0.05, 1)+
  ggtitle("Prediction Probabilities with St. Dev")+
  xlab("Plan Index")+
  ylab("Probability")+
  theme_classic()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16))+
  scale_color_manual(values=c("green", "red"))



################################################################################################################################################################################################
#Figure 6

ggplot(df_test_probs, aes(x=index, y=mean_prob_sort, color= factor( unlab_truth[index_order] ) ))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean_prob_sort - stdev[index_order], ymax=mean_prob_sort+stdev[index_order]),
                width=.9, size=1.1,
                position=position_dodge(0.5))+
  geom_hline(yintercept=c(0.5, best_threshold$threshold[1]),
             linetype='dashed', color=c( 'grey', 'black'))+
  ylim(-0.15, 1)+
  ggtitle("Ordered Probabilities with St. Dev.")+
  xlab("Ordered Plan Index")+
  ylab("Probability")+
  theme_classic()+
  theme(plot.title = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),  # Legend title
        legend.text = element_text(size = 18),   # Legend text
        axis.text=element_text(size=20))+
  labs(color = "Errors")+
  scale_color_manual(values=c("green", "red"))





# Performance result
roc.curve( df_test_probs$mean_prob[unlab_truth==1],  df_test_probs$mean_prob[unlab_truth==0] )
pr.curve( df_test_probs$mean_prob[unlab_truth==1],  df_test_probs$mean_prob[unlab_truth==0] )


# Results using Youden's index
confusionMatrix( as.factor( as.integer(df_test_probs$mean_prob> best_threshold$threshold[1]) ),
                 as.factor(unlab_truth), mode = "everything", positive="1")




rocobj <- roc( unlab_truth, df_test_probs$mean_prob+df_test_probs$stdev)
co1 = coords(rocobj, "best")
co1

co1$threshold
confusionMatrix( as.factor( as.integer(df_test_probs$mean_prob+df_test_probs$stdev> co1$threshold) ),
                 as.factor(unlab_truth), mode = "everything", positive="1")




################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################


#Understanding prediction results

obs_ids = which(df_test_probs$mean_prob+df_test_probs$stdev>0.5)
obs_ids_names = prediction_results$names[ df_test_probs$mean_prob+df_test_probs$stdev>0.5]

obs_ids = which(df_test_probs$mean_prob>best_threshold$threshold)
obs_ids_names = prediction_results$names[ df_test_probs$mean_prob>best_threshold$threshold]

var_names2[12]= "eccentricity"


################################################################################################################################################################################################
################################################################################################################################################################################################
#Shapely Value plots

shap_predictor_plot = function(df_shapely_total){
  ggplot(df_shapely_total, aes(x=name, y=shapely, fill = shapely>0)) +
    geom_bar(stat="identity") +
    theme_bw()+
    coord_flip()+
    ggtitle("Predictor Importance")+
    xlab("")+
    ylab("Shapely Value")+
    scale_fill_manual(values = c("darkblue", "skyblue"),
                      name = "Sign",
                      breaks = c(TRUE, FALSE),
                      labels = c("Positive", "Negative"))+
    theme(axis.text.x = element_text(angle = 0, hjust = 1))+
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.position = "none",
          axis.text=element_text(size=14))
}


shap_loc_mean_plot= function(fp_num, test_obs, shap_name ){

  i=fp_num

  plot_data= test[test_obs,][(num_intvls*(i-1)+1):(num_intvls*i)]

  mean(plot_data)
  min(plot_data)
  plot_points_imp = colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)]

  #good transformation for visualization
  plot_points_imp2 =mean(plot_data)*plot_points_imp/max(plot_points_imp)



  # See the mean good kidney:

  f0 = Data2fd( t(as.matrix(df[1:260 ,(num_intvls*(i-1)+1):(num_intvls*i)])))

  mean_area = mean.fd(f0)
  pd = Data2fd(plot_data)
  pd_mean= mean.fd(f0)
  df_shape_loc11 = tibble( percentage = c(0, pd$basis$params, 1),
                           line = pd$coefs[c(-1, -length(pd$coefs))] ,
                           mean_line = pd_mean$coefs[c(-1, -length(pd$coefs))],
                           shapely = plot_points_imp
  )

  df_shape_loc11

  multiplier <- max(df_shape_loc11$line)/max(df_shape_loc11$shapely)
  ggplot(df_shape_loc11, aes(x = percentage)) +
    geom_bar(aes(y = shapely*multiplier,  fill = shapely>0), stat = "identity", alpha = 0.7) +
    geom_line(aes(y = mean_line), color = "green", linetype = 2, linewidth = 2) +
    geom_line(aes(y = line), color = "red", linewidth = 1.5) +
    scale_fill_manual(values = c("darkblue", "skyblue"),
                      name = "Sign",
                      breaks = c(TRUE, FALSE),
                      labels = c("Positive", "Negative"))+
    ggtitle("")+
    ylab(" ")+
    xlab("Slice Percentage")+
    scale_y_continuous(name = shap_name,
                       sec.axis = sec_axis(~./multiplier, name=""))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 1))+
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 16),
          legend.position = "none",
          axis.text=element_text(size=18))
}




#######################################################################################################################################
#2
j=2
#Gettting shapely values
obs_ids_names[j]
obs_id = obs_ids[j]
test_obs = obs_ids[j]
obs1= data.frame(df_unlab)[obs_id,]


unif_list2=NULL
tic()
for (i in 1:ndpost) {
  bart_unify = BART.unify(m00, data.frame(df), ntree, ndpost, i )
  treeshap1 <- treeshap(bart_unify,  obs1 , verbose = F, interactions = F)
  unif_list2 = rbind( unif_list2, treeshap1$shaps)
  print(i)
}
toc()


par(mfrow=c(1,1))
barplot(colMeans(unif_list2), main = paste("Shapely Value Average Obs: ", obs_ids_names[j]))


#Per function
par(mfrow=c(2,4))
shape_spahely_total=c()
for (i in 1:num_predictors) {
  barplot(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)],
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ),
          main = paste("Obs-", obs_id, var_names2[i]),
          ylim=c(min(colMeans(unif_list2)), max(colMeans(unif_list2))))

  shape_spahely_total[i] = sum(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)])
}
par(mfrow=c(1,1))
df_shapely_total = tibble(name = var_names2, shapely = shape_spahely_total)

#Figure 7.a
shap_predictor_plot(df_shapely_total)

#Figure 7.c
shap_loc_mean_plot(3, test_obs, shap_name="Mean Radius")


### Get all the individual function positions
shap_loc_mean_plot(1, test_obs, shap_name="Area" )
shap_loc_mean_plot(2, test_obs, shap_name="Perimeter" )
shap_loc_mean_plot(3, test_obs, shap_name="Mean Radius")
shap_loc_mean_plot(4, test_obs, shap_name="SD Radius")
shap_loc_mean_plot(5, test_obs, shap_name="Min Radius")
shap_loc_mean_plot(6, test_obs, shap_name="Max Radius")
shap_loc_mean_plot(7, test_obs, shap_name="Centroid Size" )
shap_loc_mean_plot(8, test_obs, shap_name="Spphericity" )
shap_loc_mean_plot(9, test_obs, shap_name="Convexity" )
shap_loc_mean_plot(10, test_obs, shap_name="Roundness" )
shap_loc_mean_plot(11, test_obs, shap_name="Circularity" )
shap_loc_mean_plot(12, test_obs, shap_name="Eccentricity" )
shap_loc_mean_plot(13, test_obs, shap_name="Elongation" )
shap_loc_mean_plot(14, test_obs, shap_name="Rectangularity" )
shap_loc_mean_plot(15, test_obs, shap_name="Solidity" )

#######################################################################################################################################
#######################################################################################################################################

j=7
#Gettting shapely values
obs_ids_names[j]
obs_id = obs_ids[j]
test_obs = obs_ids[j]
obs1= data.frame(df_unlab)[obs_id,]


unif_list2=NULL
tic()
for (i in 1:ndpost) {
  bart_unify = BART.unify(m00, data.frame(df), ntree, ndpost, i )
  treeshap1 <- treeshap(bart_unify,  obs1 , verbose = F, interactions = F)
  unif_list2 = rbind( unif_list2, treeshap1$shaps)
  print(i)
}
toc()


par(mfrow=c(1,1))
barplot(colMeans(unif_list2), main = paste("Shapely Value Average Obs: ", obs_ids_names[j]))
summary(colMeans(unif_list2))

#Per function
par(mfrow=c(2,4))
shape_spahely_total=c()
for (i in 1:num_predictors) {
  barplot(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)],
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ),
          main = paste("Obs-", obs_id, var_names2[i]),
          ylim=c(min(colMeans(unif_list2)), max(colMeans(unif_list2))))

  shape_spahely_total[i] = sum(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)])
}


df_shapely_total = tibble(name = var_names2, shapely = shape_spahely_total)

#Figure 8.a)
shap_predictor_plot(df_shapely_total)

#Figure 8.c)
shap_loc_mean_plot(3, test_obs, shap_name="Mean Radius")


### Get all the individual function positions
shap_loc_mean_plot(1, test_obs, shap_name="Area" )
shap_loc_mean_plot(2, test_obs, shap_name="Perimeter" )
shap_loc_mean_plot(3, test_obs, shap_name="Mean Radius")
shap_loc_mean_plot(4, test_obs, shap_name="SD Radius")
shap_loc_mean_plot(5, test_obs, shap_name="Min Radius")
shap_loc_mean_plot(6, test_obs, shap_name="Max Radius")
shap_loc_mean_plot(7, test_obs, shap_name="Centroid Size" )
shap_loc_mean_plot(8, test_obs, shap_name="Sphericity" )
shap_loc_mean_plot(9, test_obs, shap_name="Convexity" )
shap_loc_mean_plot(10, test_obs, shap_name="Roundness" )
shap_loc_mean_plot(11, test_obs, shap_name="Circularity" )
shap_loc_mean_plot(12, test_obs, shap_name="Eccentricity" )
shap_loc_mean_plot(13, test_obs, shap_name="Elongation" )
shap_loc_mean_plot(14, test_obs, shap_name="Rectangularity" )
shap_loc_mean_plot(15, test_obs, shap_name="Solidity" )

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

j=9
#Gettting shapely values
obs_ids_names[j]
obs_id = obs_ids[j]
test_obs = obs_ids[j]
obs1= data.frame(df_unlab)[obs_id,]


unif_list2=NULL
tic()
for (i in 1:ndpost) {
  bart_unify = BART.unify(m00, data.frame(df), ntree, ndpost, i )
  treeshap1 <- treeshap(bart_unify,  obs1 , verbose = F, interactions = F)
  unif_list2 = rbind( unif_list2, treeshap1$shaps)
  print(i)
}
toc()


par(mfrow=c(1,1))
barplot(colMeans(unif_list2), main = paste("Shapely Value Average Obs: ", obs_ids_names[j]))
summary(colMeans(unif_list2))

#Per function
par(mfrow=c(2,4))
shape_spahely_total=c()
for (i in 1:num_predictors) {
  barplot(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)],
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ),
          main = paste("Obs-", obs_id, var_names2[i]),
          ylim=c(min(colMeans(unif_list2)), max(colMeans(unif_list2))))

  shape_spahely_total[i] = sum(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)])
}



df_shapely_total = tibble(name = var_names2, shapely = shape_spahely_total)

#Figure 9.a)
shap_predictor_plot(df_shapely_total)

#Figure 9.c)
shap_loc_mean_plot(12, test_obs, shap_name="Eccentricity" )


### Get all the individual function positions

shap_loc_mean_plot(1, test_obs, shap_name="Area" )
shap_loc_mean_plot(2, test_obs, shap_name="Perimeter" )
shap_loc_mean_plot(3, test_obs, shap_name="Mean Radius")
shap_loc_mean_plot(4, test_obs, shap_name="SD Radius")
shap_loc_mean_plot(5, test_obs, shap_name="Min Radius")
shap_loc_mean_plot(6, test_obs, shap_name="Max Radius")
shap_loc_mean_plot(7, test_obs, shap_name="Centroid Size" )
shap_loc_mean_plot(8, test_obs, shap_name="Sphericity" )
shap_loc_mean_plot(9, test_obs, shap_name="Convexity" )
shap_loc_mean_plot(10, test_obs, shap_name="Roundness" )
shap_loc_mean_plot(11, test_obs, shap_name="Circularity" )
shap_loc_mean_plot(12, test_obs, shap_name="Eccentricity" )
shap_loc_mean_plot(13, test_obs, shap_name="Elongation" )
shap_loc_mean_plot(14, test_obs, shap_name="Rectangularity" )
shap_loc_mean_plot(15, test_obs, shap_name="Solidity" )
#######################################################################################################################################












