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
roc2 <- roc(test_bap$y, test_bap$linear.predictors)
roc3 <- roc(YY, test_predict)
plot(roc1, lty=1, legacy.axes = TRUE)
plot(roc2, lty = 2, add = TRUE)
plot(roc3, lty = 3, add = TRUE)
legend(x = 0.2, y = 0.3, bty = "n", lwd = 3, lty = 1:10,
       legend = c("BAP-65","CURB-65","GXBoost"))
roc.test(roc1, roc3, method = "bootstrap", boot.n = 2000)
roc.test(roc2, roc3, method = "bootstrap", boot.n = 2000)
# COmparison of ROCS in another methods
library(fbroc)
result.boot <- boot.paired.roc(test_bap$linear.predictors, test_predict, 
                               test_bap$y, n.boot = 100)
