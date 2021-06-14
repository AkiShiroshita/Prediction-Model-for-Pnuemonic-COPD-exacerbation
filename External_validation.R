
# Setting up datasets -----------------------------------------------------

#### Reference: https://statinfer.com/203-4-2-calculating-sensitivity-and-specificity-in-r/
library(tidyverse)
library(DescTools) #Winsorize(x)
library(gt)
library(mice)
library(rms)
library(measures)
library(mfp)
library(mgcv)
library(VIM)
library(naniar)
library(Hmisc)
library(ROCR)
library(tableone)
library(caret)
library(pROC)
library(svglite)
url <- "https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz"
pkgFile <- "norm2_2.0.3.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(norm2)
### ams, bun, rr, sbp, dbp, hr, age
### BAP-65: BUN ≥ 25 mg/dL, Pulse ≥ 109 beats/min, Age ≥ 65 years
### CURB-65: BUN > 19 mg/dL, Respiratory Rate ≥ 30, Systolic BP < 90 mmHg or Diastolic BP ≤ 60 mmHg, Age ≥ 65
kobe <- read.csv("Data/kobe_cleaned.csv")
agmc <- read.csv("Data/agmc_cleaned.csv")
ichinishi <- read.csv("Data/ichinishi_cleaned.csv")
saiseikai <- read.csv("Data/saiseikai_cleaned.csv")
kmc <- read.csv("Data/kmc_cleaned.csv")
dat1 <- kobe %>% 
  select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death, los, intubation) %>% 
  mutate_all(.funs = ~ as.character(.))
dat2 <- agmc %>% 
  select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death, los, intubation) %>% 
  mutate_all(.funs = ~ as.character(.))
dat3 <- ichinishi %>% 
  select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death, los, intubation) %>% 
  mutate_all(.funs = ~ as.character(.))
dat4 <- saiseikai %>% 
  select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death, los, intubation) %>% 
  mutate_all(.funs = ~ as.character(.))
dat5 <- kmc %>% 
  select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death, los, intubation) %>% 
  mutate_all(.funs = ~ as.character(.))
dat <- bind_rows(dat1, dat2)
dat <- bind_rows(dat, dat3)
dat <- bind_rows(dat, dat4)
dat <- bind_rows(dat, dat5)
original_data <- dat %>% 
  select(-los, -intubation)
original_data %>% write_csv("original_data.csv")

# Descriptive analysis --------------------------------------------------

### Table one
dat %>% glimpse()
dat <- dat %>% 
  mutate_all(.funs = ~ as.numeric(.))
vars <- c("con", "adl", "sex", "bun", "rr", "sbp", "eo", "dbp", "hr", "age", "death", "los", "intubation")
factorVars <- c("con", "death", "sex", "intubation", "adl")
table1 <- CreateTableOne(vars = vars,
                         data = dat,
                         includeNA = TRUE,
                         factorVars = factorVars)
table1 %>% 
  print(nonnormal = "los")
### Missing pattern
miss <- miss_var_summary(dat)
miss
na.patterns <- naclus(dat)
plot(na.patterns, ylab="Fraction of NAs in common")
naplot(na.patterns)
na.pattern(dat)
aggr(dat)
### Table 1 by prognosis
table2 <- CreateTableOne(vars = vars,
                         data = dat,
                         includeNA = TRUE,
                         strata = "death",
                         factorVars = factorVars)
table2 %>% 
  print(nonnormal = "los") %>% 

# Complete case analysis --------------------------------------------------

dat %>% glimpse()
dat_cc <- dat %>% 
  na.omit()
dat_cc %>% glimpse()
dat_cc_cat <- dat_cc %>% 
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_curb = ifelse(bun >= 25, 1, 0)) %>%
  mutate(bun_bap = ifelse(bun >= 19, 1, 0)) %>% 
  mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0)) %>%
  mutate(hr = ifelse(hr >= 109, 1, 0)) %>% 
  select(death, intubation, los, age, bun_curb, bun_bap, rr, bp, con, hr) %>% 
  mutate_all(.funs = ~ as.logical(.))
dat_cc_cat %>% glimpse()
cat_vars <- c("con", "bun_curb", "bun_bap", "rr", "bp", "hr", "age", "death", "los", "intubation")
cat_factorVars <- c("con", "bun_curb", "bun_bap", "rr", "bp", "hr", "age", "death", "intubation")
table3 <- CreateTableOne(vars = cat_vars,
                         data = dat_cc_cat,
                         includeNA = TRUE,
                         strata = "death",
                         factorVars = cat_factorVars)
table3 %>% 
  print(nonnormal = "los")
### Discrimination 
### Calibration
#### Model fit
#### Brier score
#### BAP-65
ext_cc_bap <- glm(death ~ offset(1*(bun_bap + con + hr + age)),
                  data=dat_cc_cat, x=TRUE, y=TRUE) 
summary(ext_cc_bap)
probabilities = ext_cc_bap$linear.predictor
response = as.factor(as.numeric(probabilities > 0.5))
truth = ext_cc_bap$y
positive = 1
negative = 0
BrierScaled(probabilities, truth, negative, positive)
#### CURB-65
ext_cc_curb <- glm(death ~ offset(1*(con + bun_curb + rr + bp + age)),
                data=dat_cc_cat, x=TRUE, y=TRUE) 
summary(ext_cc_curb)
### Risk tables
table_curb <- CreateTableOne(vars = "death",
                             strata = "curb_score",
                             factorVars = "death",
                             data = dat_cc_cat)
table_curb
table_bap <- CreateTableOne(vars = "death",
                            strata = "bap_score",
                            factorVars = "death",
                            data = dat_cc_cat)
table_bap
#### Calibration plot
res1 <- val.prob(logit=ext_cc_bap$linear.predictor, y=ext_cc_bap$y) 
res2 <- val.prob(logit=ext_cc_curb$linear.predictor, y=ext_cc_curb$y) 
res3 <- val.prob(logit=train_predict, y=Y) 
res1[c("Intercept","Slope")]
res2[c("Intercept","Slope")]
res3[c("Intercept","Slope")]
library(generalhoslem)
logitgof(train_predict, Y, g = 10)
#### AUC
AUC <- function(xb.hat,y){
  n<-length(xb.hat)
  n1<-sum(y)
  mean.rank <- mean(rank(xb.hat)[y == 1])
  AUC<-(mean.rank - (n1 + 1)/2)/(n - n1)
  return(AUC) }
AUC(ext_cc_bap$linear.predictor, ext_cc_bap$y)
cstatNo <- rcorr.cens(ext_cc_bap$linear.predictors, ext_cc_bap$y) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
cstatYes <- rcorr.cens(ext_cc_bap$linear.predictors, ext_cc_bap$y) 
cat(cstatYes[1], "[", cstatYes[1]-1.96/2*cstatYes[3], " - ", cstatYes[1]+1.96/2*cstatYes[3],"]")
AUC(ext_cc_curb$linear.predictor, ext_cc_curb$y)
cstatNo <- rcorr.cens(ext_cc_curb$linear.predictors, ext_cc_curb$y) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
cstatYes <- rcorr.cens(ext_cc_curb$linear.predictors, ext_cc_curb$y) 
cat(cstatYes[1], "[", cstatYes[1]-1.96/2*cstatYes[3], " - ", cstatYes[1]+1.96/2*cstatYes[3],"]")
#### ROC curves
dat_cc_cat <-dat_cc_cat %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(curb_score = con + bun_curb + rr + bp + age) %>% 
  mutate(bap_score = bun_bap + con + hr + age) %>%
  mutate(curb_score = as.factor(curb_score)) %>% 
  mutate(bap_score = as.factor(bap_score))
dat_cc_cat %>% glimpse()
dat_cc_cat <-dat_cc_cat %>% 
  mutate_all(.funs = ~ as.numeric(.))
roc1 <- roc(ext_cc_bap$y, ext_cc_bap$linear.predictors)
roc2 <- roc(ext_cc_curb$y, ext_cc_curb$linear.predictors)
roc3 <- roc(Y, train_predict)
plot(roc1, lty=1, legacy.axes = TRUE)
plot(roc2, lty = 2, add = TRUE)
plot(roc3, lty = 3, add = TRUE)
legend(x = 0.2, y = 0.3, bty = "n", lwd = 3, lty = 1:10,
       legend = c("BAP-65","CURB-65","GXBoost"))
roc.test(roc1, roc3, method = "bootstrap", boot.n = 2000)
roc.test(roc2, roc3, method = "bootstrap", boot.n = 2000)
#### Sensitivity and Specificity
dat_cc_cat %>% glimpse()
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$curb_score >= 0, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$curb_score > 0, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$curb_score > 2, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$curb_score > 3, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$curb_score > 4, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
predicted_values <- ifelse(dat_cc_cat$curb_score > 5, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat %>% glimpse()
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$bap_score >= 0, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$bap_score > 0, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$bap_score > 1, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$bap_score > 2, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
dat_cc_cat <-dat_cc_cat %>%
  mutate_all(.funs = ~ as.numeric(.))
predicted_values <- ifelse(dat_cc_cat$bap_score > 3, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)
predicted_values <- ifelse(dat_cc_cat$bap_score > 4, 1, 0)
actual_values <- dat_cc_cat$death
conf_matrix <- table(predicted_values, actual_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)

# Multiple impuation ------------------------------------------------------

dat %>% glimpse()
dat <- dat %>% 
  mutate(sex = as.logical(sex),
         adl = as.logical(adl),
         con = as.logical(con),
         intubation = as.logical(intubation),
         death = as.logical(death))
dat1 <- mice(dat, maxit = 0)
dat1$method
dat1$predictorMatrix
dat2 <- mice(dat, m = 10, maxit = 20, printFlag = FALSE, seed = 1234)
plot(dat2)
dat_i <- mice(dat, m = 100, maxit = 20, print = FALSE, seed = 1234)
plot(dat_i)
se_curb <- auc_curb <- se_bap <- auc_bap <- vector(length = dat_i$m, mode = "list")
for (i in 1:dat_i$m) {
  com <- complete(dat_i, i)
  com <- com %>% 
    mutate(age = ifelse(age >= 65, 1, 0)) %>% 
    mutate(bun_curb = ifelse(bun >= 25, 1, 0)) %>%
    mutate(bun_bap = ifelse(bun >= 19, 1, 0)) %>% 
    mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
    mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0)) %>%
    mutate(hr = ifelse(hr >= 109, 1, 0)) %>% 
    mutate_all(.funs = ~ as.logical(.))
  fit_curb <- glm(death ~ offset(1*(con + bun_curb + rr + bp + age)),
                     data=com, x=TRUE, y=TRUE) 
  fit_bap <- glm(death ~ offset(1*(bun_bap + con + hr + age)),
                    data=com, x=TRUE, y=TRUE) 
  auc_curb[[i]] <- AUC(fit_curb$linear.predictor, fit_curb$y)
  auc_bap[[i]] <- AUC(fit_bap$linear.predictor, fit_bap$y)
  cstatYes_curb <- rcorr.cens(fit_curb$linear.predictors, fit_curb$y) 
  cstatYes_bap <- rcorr.cens(fit_bap$linear.predictors, fit_bap$y) 
  se_curb[[i]] <- 1.96/2*cstatYes_curb[3]
  se_bap[[i]] <- 1.96/2*cstatYes_bap[3]
}
miinf_curb <- miInference(auc_curb, se_curb)
miinf_curb
miinf_bap <- miInference(auc_bap, se_bap)
miinf_bap
