
# Set-up ---------------------------------------------------------

install.packages("tidyverse")
install.packages("xgboost")
install.packages("Ckmeans.1d.dp")
install.packages("recipes")
install.packages("caret")
install.packages("vip")
install.packages("ROCR")
install.packages("rsample")
install.packages("Hmisc")
install.packages("data.table")
install.packages("doParallel")
install.packages("pROC")
library(tidyverse)
library(xgboost)
library(Ckmeans.1d.dp)
library(recipes)
library(caret)
library(vip)
library(ROCR)
library(rsample)
library(Hmisc)
library(data.table)
library(doParallel)
library(pROC)
getwd()
original_data <- read.csv("Data/original_data.csv")
## Shuffle the original data
set.seed(19891116)
shuffled_data <- original_data[sample(nrow(original_data)),]
## Split sample
table(shuffled_data$death) %>% prop.table()
set.seed(19891116)
split_strat <- initial_split(shuffled_data, prop = 0.7,
                             strata = "death") # initial_split() split data to training/test data
train_data <- training(split_strat) # rsample::training makes training data
test_data <- testing(split_strat) # rsample::testing makes test data
table(train_data$death) %>% prop.table()
table(test_data$death) %>% prop.table()
## Converting data to suitable format for xgboost package
xgb_prep<- recipe(death ~ ., data = train_data) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = train_data, retain= TRUE) %>%
  juice () 
X <- as.matrix(xgb_prep[setdiff(names(xgb_prep), "death")]) # by setdiff(), extract columns except for "death"
Y <- xgb_prep$death

# XGBoost -----------------------------------------------------------------

set.seed(19891116)
## Setting-up hyperparameters
### eta: Learning rate, or step size shrinkage in update (0-1). Default is 0.3. Typical final values were 0.01-0.2.
### gamma: Minimum amount of splitting (0-∞). Default is 0.
### max_depth: Maximum levels deep that each tree can grow (0-∞). Default is 6. Typical values were 3-10.
### min_child_weight: minimum degree of impurity needed in a node before attempting to split it (0-∞). Default is 1. 
### subsample: The proportion of cases to be randomly sampled (without replacement) for each tree (0-1). Typical values were 0.5-1.
### colsample_bytree: The proportion of predictor variables sampled for each tree (0-1). Default is 1.
### lambda: L2 regularization term on weights. Default is 1.
### alpha: L1 regularization term on weights. Default is 0.
### gamma: minimum loss reduction required to make another leaf partition in a tree.
### tree_method string: Tree construction algorithm. Default is auto.
### scale_pos_weight: Balance of positive and negative weights. Default is 1.
### nrounds: maximum number of sequentially built trees in the model.
### eval_metric: Type of residual error/loss function. 
## Grid search for hyperparameters
hyperparameters <- xgb.cv(data = X,
                          label= Y,
                          nrounds = 5000,
                          objective= "binary:logistic",
                          early_stopping_rounds= 50,
                          nfold = 4,
                          params = list(eta= 0.01,
                                        max_depth= 3,
                                        min_child_weight= 3,
                                        subsample= 0.5,
                                        colsample_bytree= 0.5),
                          metrics = "auc",
                          verbose= 0)
max(hyperparameters[["evaluation_log"]][["test_auc_mean"]]) # Here, check the best AUC in cv
hyper_grid <- expand.grid(eta= 0.1,
                          max_depth= c(2,4,6,8,10),
                          min_child_weight= c(1,1.5,2,2.5,3),
                          subsample= 1,
                          colsample_bytree= 1,
                          gamma = 0,
                          lambda = 1,
                          alpha = 0,
                          auc = 0, 
                          trees= 0 )
for(i in seq_len(nrow(hyper_grid))){
  set.seed (19891116)
  m <- xgb.cv(data = X,
              label= Y,
              nrounds = 4000,
              objective= "binary:logistic",
              early_stopping_rounds= 50,
              nfold = 4,
              verbose= 0,
              params = list(eta = hyper_grid$eta[i],
                            max_depth = hyper_grid$max_depth[i],
                            min_child_weight = hyper_grid$min_child_weight[i],
                            subsample = hyper_grid$subsample[i],
                            colsample_bytree = hyper_grid$colsample_bytree[i],
                            gamma = hyper_grid$gamma[i],
                            lambda = hyper_grid$lambda[i],
                            alpha = hyper_grid$alpha[i]),
              metrics = "auc"
  )
  hyper_grid$auc[i] <- max(m[["evaluation_log"]][["test_auc_mean"]])
  hyper_grid$trees[i]<-m$best_iteration
} 
## Results
hyper_grid %>% 
  arrange(desc(auc)) %>%
  glimpse ()
## Visualizing the tuning
plot(log(hyper_grid$auc),
     type = "l",
     cex.lab = 0.7,
     col = "darkorange",
     xaxt = "n",
     yaxt = "n",
     xlab = "Number of trees",
     ylab = "AUC")
axis(1, col = "black",
     col.axis = "black",
     cex.axis = 0.7)
axis(2, col = "black",
     col.axis = "black",
     cex.axis = 0.7)
## Setup the optimal parameter list 
params <- list(
  eta= 0.1,
  max_depth= 2,
  min_child_weight= 2,
  subsample= 1,
  colsample_bytree= 1,
  gamma= 0,
  lambda= 1,
  alpha= 0) 
## Test in the train dataset 
xgb.train.final <- xgboost(
  params = params,
  data = X,
  label= Y,
  nrounds = 35,
  objective= "binary:logistic",
  verbose= 0)
train_predict <- predict(xgb.train.final, X)
## Performance 
AUC <- function(xb.hat,y){
  n<-length(xb.hat)
  n1<-sum(y)
  mean.rank <- mean(rank(xb.hat)[y == 1])
  AUC<-(mean.rank - (n1 + 1)/2)/(n - n1)
  return(AUC) }
AUC(train_predict, Y)
cstatNo <- rcorr.cens(train_predict, Y) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
cstatYes <- rcorr.cens(train_predict, Y) 
cat(cstatYes[1], "[", cstatYes[1]-1.96/2*cstatYes[3], " - ", cstatYes[1]+1.96/2*cstatYes[3],"]")
train_predict <- ifelse(train_predict <= 0.5, -1, 1)
confusionMatrix(as.factor(train_predict), as.factor(ifelse(Y == 0, -1, 1)))
## Test in the test dataset
xgb_test<- recipe(death ~ ., data = test_data) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = test_data, retain= TRUE) %>%
  juice () 
XX <- as.matrix(xgb_test[setdiff(names(xgb_test), "death")]) 
YY <- xgb_test$death
xgb_test <- xgboost(
  params = params,
  data = XX,
  label= YY,
  nrounds = 33,
  objective= "binary:logistic",
  verbose= 0)
test_predict <- predict(xgb_test, XX)
## Performance
AUC(test_predict, YY)
cstatNo <- rcorr.cens(test_predict, YY) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
## ROC
roc3 <- roc(YY, test_predict)
plot(roc3)
## Plot importance of variables 
vip (xgb_test)
feature_name <- c("Sex",
                  "Alterned mental Status",
                  "Blood urea nitrogen",
                  "Blood eosinophil count",
                  "Activity of daily living",
                  "Respiratory rate",
                  "Systolic blood pressure",
                  "Diastolic blood pressure",
                  "Heart rate",
                  "Age")
importance_matrix <- xgb.importance(feature_name,
                                    model = xgb_test)
xgb.plot.importance(importance_matrix[1:10, ])

