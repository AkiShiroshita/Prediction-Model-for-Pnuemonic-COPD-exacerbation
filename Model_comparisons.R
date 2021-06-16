# ROC curve
test_data <- test_data %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_curb = ifelse(bun >= 25, 1, 0)) %>%
  mutate(bun_bap = ifelse(bun >= 19, 1, 0)) %>% 
  mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0)) %>%
  mutate(hr = ifelse(hr >= 109, 1, 0))
test_data <-test_data %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(curb_score = con + bun_curb + rr + bp + age) %>% 
  mutate(bap_score = bun_bap + con + hr + age) %>%
  mutate(curb_score = as.factor(curb_score)) %>% 
  mutate(bap_score = as.factor(bap_score))
test_bap <- glm(death ~ offset(1*(bun_bap + con + hr + age)),
                data=test_data, x=TRUE, y=TRUE) 
test_curb <- glm(death ~ offset(1*(con + bun_curb + rr + bp + age)),
                 data=test_data, x=TRUE, y=TRUE) 
roc1 <- roc(test_bap$y, test_bap$linear.predictors)
roc2 <- roc(test_curb$y, test_curb$linear.predictors)
roc3 <- roc(YY, test_predict)
plot(roc1, lty=1, legacy.axes = TRUE)
plot(roc2, lty = 2, add = TRUE)
plot(roc3, lty = 3, add = TRUE)
legend(x = 0.2, y = 0.3, bty = "n", lwd = 3, lty = 1:10,
       legend = c("BAP-65","CURB-65","GXBoost"))
roc.test(roc1, roc3, method = "bootstrap", boot.n = 2000)
roc.test(roc2, roc3, method = "bootstrap", boot.n = 2000)
# Comparison between BAP-65 and XGBoost model
library(fbroc)
test_data_comparison<- test_data[complete.cases(test_data[c("age", "bun", "hr", "con", "death")]),]
xgb_test_comparison<- recipe(death ~ ., data = test_data_comparison) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = test_data_comparison, retain= TRUE) %>%
  juice () 
XXX <- as.matrix(xgb_test_comparison[setdiff(names(xgb_test_comparison), "death")]) 
YYY <- xgb_test_comparison$death
xgb_test_comparison_fit <- xgboost(
  params = params,
  data = XXX,
  label= YYY,
  nrounds = 33,
  objective= "binary:logistic",
  verbose= 0)
test_predict_comparison <- predict(xgb_test_comparison_fit, XXX)
test_data_comparison <- test_data_comparison %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_bap = ifelse(bun >= 19, 1, 0)) %>%
  mutate(hr = ifelse(hr >= 109, 1, 0))
glimpse(test_data_comparison)
test_data_comparison<- test_data_comparison %>%
  mutate(age = as.logical(age),
         hr = as.logical(hr),
         con = as.logical(con),
         bun_bap = as.logical(bun_bap),
         death = as.logical(death))
test_bap_comparison <- glm(death ~ offset(1*(bun_bap + con + hr + age)),
                data=test_data_comparison, x=TRUE, y=TRUE) 
YYY<- as.logical(YYY)
result.boot <- boot.paired.roc(test_bap_comparison$linear.predictors, test_predict_comparison, 
                               YYY, n.boot = 100)
perf(result.boot, "auc")
# Comparison between CURB-65 and XGBoost model
test_data_comparison<- test_data[complete.cases(test_data[c("age", "bun", "rr", "con", "sbp", "dbp", "death")]),]
xgb_test_comparison<- recipe(death ~ ., data = test_data_comparison) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = test_data_comparison, retain= TRUE) %>%
  juice () 
XXX <- as.matrix(xgb_test_comparison[setdiff(names(xgb_test_comparison), "death")]) 
YYY <- xgb_test_comparison$death
xgb_test_comparison_fit <- xgboost(
  params = params,
  data = XXX,
  label= YYY,
  nrounds = 33,
  objective= "binary:logistic",
  verbose= 0)
test_predict <- predict(xgb_test_comparison_fit, XXX)
test_data_comparison <- test_data_comparison %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_curb = ifelse(bun >= 25, 1, 0)) %>%
  mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0)) 
test_data_comparison<- test_data_comparison %>%
  mutate(age = as.logical(age),
         rr = as.logical(hr),
         bp = as.logical(bp),
         con = as.logical(con),
         bun_curb = as.logical(bun_curb),
         death = as.logical(death))
test_curb_comparison <- glm(death ~ offset(1*(con + bun_curb + rr + bp + age)),
                            data=test_data_comparison, x=TRUE, y=TRUE) 
YYY<- as.logical(YYY)
result.boot <- boot.paired.roc(test_curb_comparison$linear.predictors, test_predict, 
                               YYY, n.boot = 100)
perf(result.boot, "auc")
