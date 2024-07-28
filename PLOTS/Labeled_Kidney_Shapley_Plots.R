# Shapely Value plots for Labeled Kidney feature importance

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

#Data cleaning
#colSums(is.na(goodk))
#goodk[is.na(goodk)] <- mean(goodk$sphericity1, na.rm=TRUE)
#colSums(is.na(badk))


set_good <- split(goodk,goodk$name)
set_bad <- split(badk,badk$name)
factor_good = factor(goodk$name)
factor_bad = factor(badk$name)
#Gives the plan name

set_names <- c(levels(factor_good), levels(factor_bad))

#Extracts the coefficients at specif intervals across all funcitons. This case gives 11 values
#per shape feature

grid = seq(.02, from=0, to=1)
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


#####################################################################################################
#####################################################################################################
# View Shape functions

# See Perimeter!
par(mfrow=c(1,1))
plot(set_good[[3]][[3]], xlim=c(0,80), ylim=c(0,4000),
     type = "l", xlab="Slice Index", ylab = "Area", main = "Area per Image Slice")
for (j in 1:length(set_good)) {

  lines(set_good[[j]][[3]], col = "green")

}
for (j in 1:length(set_bad)) {
  lines(set_bad[[j]][[3]], col = "red")
}




num_functions = length(var_names2)
tempbasis = create.bspline.basis(rangeval = c(0,1), nbasis=16)
grid <- seq(0,1,1/(length(data_set[1,1:numPoinst])-1))

templist      = vector("list",12)
templist[[1]] = rep(1,312)

# For all 11 shape features in the data set
for (i in 0:(num_functions-1)) {
  tempSmooth= smooth.basis(grid, t(data_set[,(numPoinst*i+1):(numPoinst*i+numPoinst)]), tempbasis)
  templist[[i+2]] = tempSmooth$fd

}

#Figure 1.c)
outcomes = c(rep(0,260), rep(1, 52))


dev.off()
par(cex.lab = 1.5,  # Axis labels size
    cex.axis = 1.75, # Axis tick labels size
    cex.main = 2, # Main title size
    cex.sub = 2)  # Subtitle size

plot(templist[[2]],
     ylab = "Area",
     xlab = "Slice Percentage",
     cex.lab = 10,   # Increase the size of the axis labels
     cex.axis = 2,  # Increase the size of the axis tick labels
     cex.main = 2,  # Increase the size of the main title
     cex.sub = 2    # Increase the size of the subtitle (if any)
)
for(i in 1:length(outcomes)){
  col1 = "green"
  if(outcomes[i] == 1){
    col1 = "red"
  }
  lines(templist[[2]][i,], col =col1 )
}








#####################################################################################################
#####################################################################################################





#create Train and Test sets
df = as_tibble(data_set)

outcomes = c(rep(0,260), rep(1, 52))
outcomes_class = c(rep("Acceptable",260), rep("Unacceptable", 52))
true_y = outcomes

# Testing out FBART



num_intvls = dim(data_set)[2]/num_predictors

points = data_set

n_point =  dim(data_set)[2]/num_predictors

#c(ntree=70, ndpost=500, binaryOffset=20, kk=2, power=2, base=0.5)

#Logistic Regression

set.seed(1)
sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.8,0.2))
train  <- as.matrix(df[sample, ])
test   <- as.matrix(df[!sample, ])
y_train = true_y[sample]
y_test = true_y[!sample]

y_test_class = outcomes_class[!sample]

which(!sample)

### Best Case!
binaryOffset = NULL
rho_int = 1
sigma_int = 1

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
                       mean_prob = m00$prob.test.mean,
                       upper_ci = upper_ci,
                       lower_ci =lower_ci,
                       stdev = stdev,
                       class = y_test_class)


df_test_probs



par(mfrow=c(1,1))
barplot(colMeans(m00$varprob_fp),  main="fBART functional predictors prob.",
        space=0, names.arg  = paste0(var_names2), las=2)

var_names2[12]= "eccentricity"

df_fp = tibble( names = var_names2, probability= colMeans(m00$varprob_fp) )
# Create a bar plot with car names tilted at a 45-degree angle

ggplot(df_fp, aes(x=names, y=probability)) +
  geom_bar(stat="identity") +
  ggtitle("Functional Predictor Importance")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=16))



var_names2

df_int_14 = tibble( names = as.character( round((1:num_intvls)/num_intvls, 2 ) )
                    , probability= colMeans(m00$varprob_times[[14]]))
# Create a bar plot with car names tilted at a 45-degree angle

selected_int = as.character( round((1:num_intvls)/num_intvls, 2 ))[c(1, 10*(1:5), num_intvls ) ]

ggplot(df_int_14, aes(x=names, y=probability)) +
  geom_bar(stat="identity") +
  ggtitle("Rectangularity")+
  xlab("Slice Percentage")+
  ylim(0,0.3)+
  theme_bw()+
  scale_x_discrete(breaks = selected_int ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=14))




df_int_7 = tibble( names = as.character( round((1:num_intvls)/num_intvls, 2 ) )
                   , probability= colMeans(m00$varprob_times[[7]]))
# Create a bar plot with car names tilted at a 45-degree angle

selected_int = as.character( round((1:num_intvls)/num_intvls, 2 ))[c(1, 10*(1:5), num_intvls ) ]

ggplot(df_int_7, aes(x=names, y=probability)) +
  geom_bar(stat="identity") +
  ggtitle("Centroid Size")+
  xlab("Slice Percentage")+
  ylim(0,0.3)+
  theme_bw()+
  scale_x_discrete(breaks = selected_int ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=14))


df_int_1 = tibble( names = as.character( round((1:num_intvls)/num_intvls, 2 ) )
                   , probability= colMeans(m00$varprob_times[[1]]))
# Create a bar plot with car names tilted at a 45-degree angle

selected_int = as.character( round((1:num_intvls)/num_intvls, 2 ))[c(1, 10*(1:5), num_intvls ) ]

ggplot(df_int_1, aes(x=names, y=probability)) +
  geom_bar(stat="identity") +
  ggtitle("Area")+
  xlab("Slice Percentage")+
  ylim(0,0.3)+
  theme_bw()+
  scale_x_discrete(breaks = selected_int ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text=element_text(size=14))



library(pROC)
rocobj <- roc(y_test, m00$prob.test.mean)
coords(rocobj)
best_threshold = coords(rocobj, "best")
best_threshold


table(y_test, m00$prob.test.mean>0.5)

#spec
46/46
#acc
49/51
#sens
3/(5)

table(y_test, m00$prob.test.mean>best_threshold$threshold)
#spec
39/46
#acc
44/51
#sens
5/(5)


ggplot(df_test_probs, aes(x=index, y=mean_prob,
                          color= class))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1,
                position=position_dodge(0.05))+
  geom_hline(yintercept=c(0.5, best_threshold$threshold),
             linetype='dashed', color=c( 'black', 'grey'))+
  ylim(-0.05, 1)+
  ggtitle("Prediction Probabilities with CI 95%")+
  xlab("Plan Index")+
  ylab("Probability")+
  theme_classic()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        axis.text=element_text(size=16))+
  scale_color_manual(values=c("green", "red"))


ggplot(df_test_probs, aes(x=index, y=mean_prob,
                          color= class))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_prob - stdev, ymax=mean_prob+stdev), width=1,
                position=position_dodge(0.05))+
  geom_hline(yintercept=c(0.5, best_threshold$threshold),
             linetype='dashed', color=c( 'black', 'grey'))+
  ylim(-0.05, 1)+
  ggtitle("Prediction Probabilities with St.Dev.")+
  xlab("Plan Index")+
  ylab("Probability")+
  theme_classic()+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        axis.text=element_text(size=16))+
  scale_color_manual(values=c("green", "red"))


#See Global

par(mfrow=c(2,4))
for (i in 1:num_predictors) {
  barplot(colMeans(m00$varprob_times[[i]]), main=paste("", var_names2[i]),
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ) , space=0)
}

######################################################################################################################################################################
######################################################################################################################################################################
#Plot Functions

shap_predictor_plot = function(df_shapely_total){
  ggplot(df_shapely_total, aes(x=name, y=shapely, fill = shapely>0)) +
    geom_bar(stat="identity") +
    theme_bw()+
    coord_flip()+
    ggtitle("Predictor Importance")+
    ylab("")+
    xlab("")+
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

#####################################################################################################################################
#####################################################################################################################################
# FIGURE 4

# Get the names
obs_names =   set_names[which(!sample)[which( m00$prob.test.mean>0.5)]]
obs_names
obs_ids  =  which(!sample)[ which( m00$prob.test.mean>0.5) ]
obs_ids
obs_id = obs_ids[1]

obs1= data.frame(df)[obs_id,]
unif_list2=NULL
tic()
for (i in 1:ndpost) {
  bart_unify = BART.unify(m00, data.frame(df), ntree, ndpost, i )
  treeshap1 <- treeshap(bart_unify,  obs1 , verbose = F, interactions = F)
  unif_list2 = rbind( unif_list2, treeshap1$shaps)
  print(i)
}
toc()

barplot(colMeans(unif_list2), main = paste("Shapely Value Average Obs: ", obs_id))
summary(colMeans(unif_list2))

#See each location shapely value per shape function
shape_spahely_total=c()
par(mfrow=c(2,4))
for (i in 1:num_predictors) {
  barplot(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)],
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ),
          main = paste("Obs-",obs_id, var_names2[i]),
          ylim=c(min(colMeans(unif_list2)), max(colMeans(unif_list2))))

  shape_spahely_total[i] = sum(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)])
}
par(mfrow=c(1,1))

#Figure 4.a
df_shapely_total = tibble(name = var_names2, shapely = shape_spahely_total)
shap_predictor_plot(df_shapely_total)


#Figure 4.c
fp_num=7
test_obs = which( m00$prob.test.mean>0.5)[1]
shap_loc_mean_plot(fp_num, test_obs, shap_name="Centroid Size" )


#####################################################################################################################################
#####################################################################################################################################
# FIGURE 5

# Get the names
obs_names =   set_names[which(!sample)[which( m00$prob.test.mean>0.5)]]
obs_ids  =  which(!sample)[ which( m00$prob.test.mean>0.5) ]
obs_id = obs_ids[3]

obs1= data.frame(df)[obs_id,]
unif_list2=NULL
#See how long shapley values take to calculate
tic()
for (i in 1:ndpost) {
  bart_unify = BART.unify(m00, data.frame(df), ntree, ndpost, i )
  treeshap1 <- treeshap(bart_unify,  obs1 , verbose = F, interactions = F)
  unif_list2 = rbind( unif_list2, treeshap1$shaps)
  print(i)
}
toc()


#See each location shapely value per shape function
shape_spahely_total=c()
par(mfrow=c(2,4))
for (i in 1:num_predictors) {
  barplot(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)],
          names.arg  = as.character( round((1:num_intvls)/num_intvls, 2 ) ),
          main = paste("Obs-",obs_id, var_names2[i]),
          ylim=c(min(colMeans(unif_list2)), max(colMeans(unif_list2))))

  shape_spahely_total[i] = sum(colMeans(unif_list2)[(num_intvls*(i-1)+1):(num_intvls*i)])
}
par(mfrow=c(1,1))

#Figure 5.a)
df_shapely_total = tibble(name = var_names2, shapely = shape_spahely_total)
shap_predictor_plot(df_shapely_total)

#Figure 5.c)
fp_num=14
test_obs = which( m00$prob.test.mean>0.5)[3]

shap_loc_mean_plot(fp_num, test_obs, shap_name="Rectangularity" )




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

df_shape_loc11$line = df_shape_loc11$line - 0.7
df_shape_loc11$mean_line = df_shape_loc11$mean_line - 0.7



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
                     labels = function(x) x + 0.7,
                     sec.axis = sec_axis(~./multiplier, name=""))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        axis.text=element_text(size=18))







########################################

df_shape_loc11

ggplot(df_shape_loc11, aes(x = percentage)) +
  geom_bar(aes(y = shapely * multiplier, fill = shapely > 0), stat = "identity", alpha = 0.7) +
  geom_line(aes(y = mean_line), color = "green", linetype = 2, linewidth = 2) +
  geom_line(aes(y = line), color = "red", linewidth = 1.5) +
  scale_fill_manual(values = c("darkblue", "skyblue"),
                    name = "Sign",
                    breaks = c(TRUE, FALSE),
                    labels = c("Positive", "Negative")) +
  ggtitle("") +
  xlab("Slice Percentage") +
  scale_y_continuous(name = shap_name,
                     limits = c(0, 0.1),  # This keeps the primary axis unchanged
                     breaks = seq(0, 0.1, by = 0.02),
                     labels = function(x) x + 0.7,  # Add 0.7 to labels
                     sec.axis = sec_axis(~ . / multiplier)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        axis.text = element_text(size = 18))












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

