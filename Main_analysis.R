
# Set-up ------------------------------------------------------------------

library(tidyverse)
library(xgboost)
library(Ckmeans.1d.dp)
library(recipes)
library(gt)
library(pracma)
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
library(rsample)
library(data.table)
library(doParallel)
library(pROC)
library(norm2)
library(boot)
library(rms)
## Preparing data
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

# Complete case analysis preparation --------------------------------------------------

## Table one
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
## Missing pattern
miss <- miss_var_summary(dat)
miss
na.patterns <- naclus(dat)
plot(na.patterns, ylab="Fraction of NAs in common")
naplot(na.patterns)
na.pattern(dat)
aggr(dat)
## Table 1 by prognosis
table2 <- CreateTableOne(vars = vars,
                         data = dat,
                         includeNA = TRUE,
                         strata = "death",
                         factorVars = factorVars)
table2 %>% 
  print(nonnormal = "los") 


# Complete case - BAP-65 --------------------------------------------------

## BAP-65: BUN ≥ 25 mg/dL, Pulse ≥ 109 beats/min, Age ≥ 65 years
dat %>% glimpse()
dat_bap_cc<- dat[complete.cases(dat[c("age", "bun", "hr", "con", "death")]),]
dat_bap_cc %>% glimpse()
dat_bap_cc <- dat_bap_cc %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_bap = ifelse(bun >= 25, 1, 0)) %>% 
  mutate(hr = ifelse(hr >= 109, 1, 0))
dat_bap_cc %>% glimpse()
dat_bap_cc <- dat_bap_cc %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(bap_score = bun_bap + con + hr + age) %>%
  mutate(bap_score = as.factor(bap_score))
dat_bap_cc %>% glimpse()
cat_vars <- c("con", "bun_bap", "hr", "age", "bap_score", "death", "los", "intubation")
cat_factorVars <- c("con", "bun_bap", "hr", "age", "bap_score", "death", "intubation")
table3 <- CreateTableOne(vars = cat_vars,
                         data = dat_bap_cc,
                         includeNA = TRUE,
                         strata = "death",
                         factorVars = cat_factorVars)
table3 %>% 
  print(nonnormal = "los")
## Discrimination of BAP-65 
dat_bap_cc <-dat_bap_cc %>%
  mutate_all(.funs = ~ as.numeric(.))
### Sensitivity and Specificity at different cut-off points
# Cut-off = 0
predicted_values0 <- ifelse(dat_bap_cc$bap_score >= 0, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix0 <- table(factor(predicted_values0, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr0 = conf_matrix0[2,2]/(conf_matrix0[2,1] + conf_matrix0[2,2])
fpr0 = conf_matrix0[2,1]/(conf_matrix0[2,1] + conf_matrix0[2,2])
x0 <- sensitivity(conf_matrix0)
y0 <- 1-specificity(conf_matrix0)
# Cut-off = 1
predicted_values1 <- ifelse(dat_bap_cc$bap_score >= 1, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix1 <- table(factor(predicted_values1, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr1 = conf_matrix1[2,2]/(conf_matrix1[2,1] + conf_matrix1[2,2])
fpr1 = conf_matrix1[2,1]/(conf_matrix1[2,1] + conf_matrix1[2,2])
x1 <- sensitivity(conf_matrix1)
y1 <- 1-specificity(conf_matrix1)
# Cut-off = 2
predicted_values2 <- ifelse(dat_bap_cc$bap_score >= 2, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix2 <- table(factor(predicted_values2, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr2 = conf_matrix2[2,2]/(conf_matrix2[2,1] + conf_matrix2[2,2])
fpr2 = conf_matrix2[2,1]/(conf_matrix2[2,1] + conf_matrix2[2,2])
x2 <- sensitivity(conf_matrix2)
y2 <- 1-specificity(conf_matrix2)
# Cut-off = 3
predicted_values3 <- ifelse(dat_bap_cc$bap_score >= 3, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix3 <- table(factor(predicted_values3, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr3 = conf_matrix3[2,2]/(conf_matrix3[2,1] + conf_matrix3[2,2])
fpr3 = conf_matrix3[2,1]/(conf_matrix3[2,1] + conf_matrix3[2,2])
x3 <- sensitivity(conf_matrix3)
y3 <- 1-specificity(conf_matrix3)
# Cut-off = 4
predicted_values4 <- ifelse(dat_bap_cc$bap_score >= 4, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix4 <- table(factor(predicted_values4, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr4 = conf_matrix4[2,2]/(conf_matrix4[2,1] + conf_matrix4[2,2])
fpr4 = conf_matrix4[2,1]/(conf_matrix4[2,1] + conf_matrix4[2,2])
x4 <- sensitivity(conf_matrix4)
y4 <- 1- specificity(conf_matrix4)
# Cut-off = 5
predicted_values5 <- ifelse(dat_bap_cc$bap_score >= 5, 1, 0)
actual_values <- dat_bap_cc$death
conf_matrix5 <- table(factor(predicted_values5, levels = 0:1),
                      factor(actual_values, levels = 0:1))
tpr5 = conf_matrix5[2,2]/(conf_matrix5[2,1] + conf_matrix5[2,2])
fpr5 = conf_matrix5[2,1]/(conf_matrix5[2,1] + conf_matrix5[2,2])
x5 <- sensitivity(conf_matrix5)
y5 <- 1-specificity(conf_matrix5)
# Describing ROC curve
roc_bap <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5),
                      Sensitivity = c(x0, x1, x2, x3, x4, x5),
                      Specificity = c(y0, y1, y2, y3, y4, y5))
trapz(roc_bap$Specificity, roc_bap$Sensitivity) # Trapezoid integration
plot(Sensitivity ~ Specificity,
     data = roc_bap,
     xlim = c(0,1),
     ylim = c(0,1),
     type="b",
     bg = "black") #Just for confirmation!! Look at the x-axis.
## 95% confidence interval
boot_num <- 2000
roc_bap_list <- c()
for(i in seq_len(boot_num)){
  set.seed(i)
  bap_data_bootstrap <-
    dat_bap_cc[sample(dim(dat_bap_cc)[1],dim(dat_bap_cc)[1],replace=TRUE),]  
  # Cut-off = 0
  predicted_values0 <- ifelse(bap_data_bootstrap$bap_score >= 0, 1, 0)
  actual_values <- bap_data_bootstrap$death
  conf_matrix0 <- table(factor(predicted_values0, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr0 = conf_matrix0[2,2]/(conf_matrix0[2,1] + conf_matrix0[2,2])
  fpr0 = conf_matrix0[2,1]/(conf_matrix0[2,1] + conf_matrix0[2,2])
  x0 <- sensitivity(conf_matrix0)
  y0 <- 1-specificity(conf_matrix0)
  # Cut-off = 1
  predicted_values1 <- ifelse(bap_data_bootstrap$bap_score >= 1, 1, 0)
  conf_matrix1 <- table(factor(predicted_values1, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr1 = conf_matrix1[2,2]/(conf_matrix1[2,1] + conf_matrix1[2,2])
  fpr1 = conf_matrix1[2,1]/(conf_matrix1[2,1] + conf_matrix1[2,2])
  x1 <- sensitivity(conf_matrix1)
  y1 <- 1-specificity(conf_matrix1)
  # Cut-off = 2
  predicted_values2 <- ifelse(bap_data_bootstrap$bap_score >= 2, 1, 0)
  conf_matrix2 <- table(factor(predicted_values2, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr2 = conf_matrix2[2,2]/(conf_matrix2[2,1] + conf_matrix2[2,2])
  fpr2 = conf_matrix2[2,1]/(conf_matrix2[2,1] + conf_matrix2[2,2])
  x2 <- sensitivity(conf_matrix2)
  y2 <- 1-specificity(conf_matrix2)
  # Cut-off = 3
  predicted_values3 <- ifelse(bap_data_bootstrap$bap_score >= 3, 1, 0)
  conf_matrix3 <- table(factor(predicted_values3, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr3 = conf_matrix3[2,2]/(conf_matrix3[2,1] + conf_matrix3[2,2])
  fpr3 = conf_matrix3[2,1]/(conf_matrix3[2,1] + conf_matrix3[2,2])
  x3 <- sensitivity(conf_matrix3)
  y3 <- 1-specificity(conf_matrix3)
  # Cut-off = 4
  predicted_values4 <- ifelse(bap_data_bootstrap$bap_score >= 4, 1, 0)
  conf_matrix4 <- table(factor(predicted_values4, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr4 = conf_matrix4[2,2]/(conf_matrix4[2,1] + conf_matrix4[2,2])
  fpr4 = conf_matrix4[2,1]/(conf_matrix4[2,1] + conf_matrix4[2,2])
  x4 <- sensitivity(conf_matrix4)
  y4 <- 1- specificity(conf_matrix4)
  # Cut-off = 5
  predicted_values5 <- ifelse(bap_data_bootstrap$bap_score >= 5, 1, 0)
  conf_matrix5 <- table(factor(predicted_values5, levels = 0:1),
                        factor(actual_values, levels = 0:1))
  tpr5 = conf_matrix5[2,2]/(conf_matrix5[2,1] + conf_matrix5[2,2])
  fpr5 = conf_matrix5[2,1]/(conf_matrix5[2,1] + conf_matrix5[2,2])
  x5 <- sensitivity(conf_matrix5)
  y5 <- 1-specificity(conf_matrix5)
  # AUC
  roc_bap <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5),
                        Sensitivity = c(x0, x1, x2, x3, x4, x5),
                        Specificity = c(y0, y1, y2, y3, y4, y5))
  roc_bap_list <- c(roc_bap_list, trapz(roc_bap$Specificity, roc_bap$Sensitivity)) 
}
quantile(roc_bap_list,c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))


# Complete case - CURB-65 -------------------------------------------------

## CURB-65: BUN > 19 mg/dL, Respiratory Rate ≥ 30, Systolic BP < 90 mmHg or Diastolic BP ≤ 60 mmHg, Age ≥ 65
dat %>% glimpse()
dat_curb_cc<- dat[complete.cases(dat[c("age", "bun", "rr", "sbp", "dbp", "death")]),]
dat_curb_cc %>% glimpse()
dat_curb_cc <- dat_curb_cc %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_curb = ifelse(bun >= 19, 1, 0)) %>% 
  mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0))
dat_curb_cc %>% glimpse()
dat_curb_cc <- dat_curb_cc %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(curb_score = con + bun_curb + rr + bp + age) %>%
  mutate(curb_score = as.factor(curb_score))
dat_curb_cc %>% glimpse()
cat_vars <- c("bun_curb", "rr", "bp", "age", "curb_score", "death", "los", "intubation")
cat_factorVars <- c("con", "bun_curb", "rr", "bp", "hr", "age", "bap_score", "death", "intubation")
table4 <- CreateTableOne(vars = cat_vars,
                         data = dat_curb_cc,
                         includeNA = TRUE,
                         strata = "death",
                         factorVars = cat_factorVars)
table4 %>% 
  print(nonnormal = "los")
## Discrimination of CURB-65 
dat_curb_cc <-dat_curb_cc %>%
  mutate_all(.funs = ~ as.numeric(.))
### true positive and false positive rate at different cut-off points
# Cut-off = 0
predicted_values00 <- ifelse(dat_curb_cc$curb_score >= 0, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix00 <- table(factor(predicted_values00, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr00 = conf_matrix00[2,2]/(conf_matrix00[2,1] + conf_matrix00[2,2])
fpr00 = conf_matrix00[2,1]/(conf_matrix00[2,1] + conf_matrix00[2,2])
x00 <- sensitivity(conf_matrix00)
y00 <- 1-specificity(conf_matrix00)
# Cut-off = 1
predicted_values01 <- ifelse(dat_curb_cc$curb_score >= 1, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix01 <- table(factor(predicted_values01, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr01 = conf_matrix01[2,2]/(conf_matrix01[2,1] + conf_matrix01[2,2])
fpr01 = conf_matrix01[2,1]/(conf_matrix01[2,1] + conf_matrix01[2,2])
x01 <- sensitivity(conf_matrix01)
y01 <- 1-specificity(conf_matrix01)
# Cut-off = 2
predicted_values02 <- ifelse(dat_curb_cc$curb_score >= 2, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix02 <- table(factor(predicted_values02, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr02 = conf_matrix02[2,2]/(conf_matrix02[2,1] + conf_matrix02[2,2])
fpr02 = conf_matrix02[2,1]/(conf_matrix02[2,1] + conf_matrix02[2,2])
x02 <- sensitivity(conf_matrix02)
y02 <- 1-specificity(conf_matrix02)
# Cut-off = 3
predicted_values03 <- ifelse(dat_curb_cc$curb_score >= 3, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix03 <- table(factor(predicted_values03, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr03 = conf_matrix03[2,2]/(conf_matrix03[2,1] + conf_matrix03[2,2])
fpr03 = conf_matrix03[2,1]/(conf_matrix03[2,1] + conf_matrix03[2,2])
x03 <- sensitivity(conf_matrix03)
y03 <- 1-specificity(conf_matrix03)
# Cut-off = 4
predicted_values04 <- ifelse(dat_curb_cc$curb_score >= 4, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix04 <- table(factor(predicted_values04, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr04 = conf_matrix04[2,2]/(conf_matrix04[2,1] + conf_matrix04[2,2])
fpr04 = conf_matrix04[2,1]/(conf_matrix04[2,1] + conf_matrix04[2,2])
x04 <- sensitivity(conf_matrix04)
y04 <- 1- specificity(conf_matrix04)
# Cut-off = 5
predicted_values05 <- ifelse(dat_curb_cc$curb_score >= 5, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix05 <- table(factor(predicted_values05, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr05 = conf_matrix05[2,2]/(conf_matrix05[2,1] + conf_matrix05[2,2])
fpr05 = conf_matrix05[2,1]/(conf_matrix05[2,1] + conf_matrix05[2,2])
x05 <- sensitivity(conf_matrix05)
y05 <- 1-specificity(conf_matrix05)
# Cut-off = 6
predicted_values06 <- ifelse(dat_curb_cc$curb_score >= 6, 1, 0)
actual_values <- dat_curb_cc$death
conf_matrix06 <- table(factor(predicted_values06, levels = 0:1),
                       factor(actual_values, levels = 0:1))
tpr06 = conf_matrix06[2,2]/(conf_matrix06[2,1] + conf_matrix06[2,2])
fpr06 = conf_matrix06[2,1]/(conf_matrix06[2,1] + conf_matrix06[2,2])
y06 <- sensitivity(conf_matrix06)
x06 <- 1- specificity(conf_matrix06)
# Describing ROC curve
roc_curb <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5, 6),
                       Sensitivity = c(x00, x01, x02, x03, x04, x05, x06),
                       Specificity = c(y00, y01, y02, y03, y04, y05, y06))
trapz(roc_curb$Specificity, roc_curb$Sensitivity) # Trapezoid integration
plot(Sensitivity ~ Specificity,
     data = roc_curb,
     xlim = c(0,1),
     ylim = c(0,1),
     type="b",
     bg = "black") #Just for confirmation!! Look at the x-axis.
## 95% CI
## 95% confidence interval
boot_num <- 2000
roc_curb_list <- c()
for(i in seq_len(boot_num)){
  set.seed(i)
  curb_data_bootstrap <-
    dat_curb_cc[sample(dim(dat_curb_cc)[1],dim(dat_curb_cc)[1],replace=TRUE),]  
  # Cut-off = 0
  predicted_values00 <- ifelse(curb_data_bootstrap$curb_score >= 0, 1, 0)
  actual_values <- curb_data_bootstrap$death
  conf_matrix00 <- table(factor(predicted_values00, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr00 = conf_matrix00[2,2]/(conf_matrix00[2,1] + conf_matrix00[2,2])
  fpr00 = conf_matrix00[2,1]/(conf_matrix00[2,1] + conf_matrix00[2,2])
  x00 <- sensitivity(conf_matrix00)
  y00 <- 1-specificity(conf_matrix00)
  # Cut-off = 1
  predicted_values01 <- ifelse(curb_data_bootstrap$curb_score >= 1, 1, 0)
  conf_matrix01 <- table(factor(predicted_values01, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr01 = conf_matrix01[2,2]/(conf_matrix01[2,1] + conf_matrix01[2,2])
  fpr01 = conf_matrix01[2,1]/(conf_matrix01[2,1] + conf_matrix01[2,2])
  x01 <- sensitivity(conf_matrix01)
  y01 <- 1-specificity(conf_matrix01)
  # Cut-off = 2
  predicted_values02 <- ifelse(curb_data_bootstrap$curb_score >= 2, 1, 0)
  conf_matrix02 <- table(factor(predicted_values02, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr02 = conf_matrix02[2,2]/(conf_matrix02[2,1] + conf_matrix02[2,2])
  fpr02 = conf_matrix02[2,1]/(conf_matrix02[2,1] + conf_matrix02[2,2])
  x02 <- sensitivity(conf_matrix02)
  y02 <- 1-specificity(conf_matrix02)
  # Cut-off = 3
  predicted_values03 <- ifelse(curb_data_bootstrap$curb_score >= 3, 1, 0)
  conf_matrix03 <- table(factor(predicted_values03, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr03 = conf_matrix03[2,2]/(conf_matrix03[2,1] + conf_matrix03[2,2])
  fpr03 = conf_matrix03[2,1]/(conf_matrix03[2,1] + conf_matrix03[2,2])
  x03 <- sensitivity(conf_matrix03)
  y03 <- 1-specificity(conf_matrix03)
  # Cut-off = 4
  predicted_values04 <- ifelse(curb_data_bootstrap$curb_score >= 4, 1, 0)
  conf_matrix04 <- table(factor(predicted_values04, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr04 = conf_matrix04[2,2]/(conf_matrix04[2,1] + conf_matrix04[2,2])
  fpr04 = conf_matrix04[2,1]/(conf_matrix04[2,1] + conf_matrix04[2,2])
  x04 <- sensitivity(conf_matrix04)
  y04 <- 1- specificity(conf_matrix04)
  # Cut-off = 5
  predicted_values05 <- ifelse(curb_data_bootstrap$curb_score >= 5, 1, 0)
  conf_matrix05 <- table(factor(predicted_values05, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr05 = conf_matrix05[2,2]/(conf_matrix05[2,1] + conf_matrix05[2,2])
  fpr05 = conf_matrix05[2,1]/(conf_matrix05[2,1] + conf_matrix05[2,2])
  x05 <- sensitivity(conf_matrix05)
  y05 <- 1-specificity(conf_matrix05)
  # Cut-off = 6
  predicted_values06 <- ifelse(curb_data_bootstrap$curb_score >= 6, 1, 0)
  conf_matrix06 <- table(factor(predicted_values06, levels = 0:1),
                         factor(actual_values, levels = 0:1))
  tpr06 = conf_matrix06[2,2]/(conf_matrix06[2,1] + conf_matrix06[2,2])
  fpr06 = conf_matrix06[2,1]/(conf_matrix06[2,1] + conf_matrix06[2,2])
  y06 <- sensitivity(conf_matrix06)
  x06 <- 1- specificity(conf_matrix06)
  # AUC
  roc_curb <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5, 6),
                         Sensitivity = c(x00, x01, x02, x03, x04, x05, x06),
                         Specificity = c(y00, y01, y02, y03, y04, y05, y06))
  trapz(roc_curb$Specificity, roc_curb$Sensitivity) # Trapezoid integration
  roc_curb_list <- c(roc_curb_list, trapz(roc_curb$Specificity, roc_curb$Sensitivity)) 
}
quantile(roc_curb_list,c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))


# Multiple imputation -----------------------------------------------------

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
auc_bap <- se_bap <- auc_curb <- se_curb <- vector(length = dat_i$m, mode = "list")
## BAP-65
# Start
for (i in 1:dat_i$m) {
  # Data preparation
  com <- complete(dat_i, i)
  com <- com %>% 
    mutate_all(.funs = ~ as.numeric(.)) %>%  
    mutate(age = ifelse(age >= 65, 1, 0)) %>% 
    mutate(bun_bap = ifelse(bun >= 25, 1, 0)) %>% 
    mutate(hr = ifelse(hr >= 109, 1, 0)) %>% 
    mutate_all(.funs = ~ as.numeric(.)) %>% 
    mutate(bap_score = bun_bap + con + hr + age) 
  # Discrimination of BAP-65 
  # AUC
  # Bootstrap
  boot_num <- 2000
  roc_bap_list_b <- c()
  for(s in seq_len(boot_num)){
    set.seed(s)
    bap_data_bootstrap <-
      com[sample(dim(com)[1],dim(com)[1],replace=TRUE),]  
    # Cut-off = 0
    predicted_values_b0 <- ifelse(bap_data_bootstrap$bap_score >= 0, 1, 0)
    actual_values_b <- bap_data_bootstrap$death
    conf_matrix_b0 <- table(factor(predicted_values_b0, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr_b0 = conf_matrix_b0[2,2]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    fpr_b0 = conf_matrix_b0[2,1]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    x_b0 <- sensitivity(conf_matrix_b0)
    y_b0 <- 1-specificity(conf_matrix_b0)
    # Cut-off = 1
    predicted_values_b1 <- ifelse(bap_data_bootstrap$bap_score >= 1, 1, 0)
    conf_matrix_b1 <- table(factor(predicted_values_b1, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr_b1 = conf_matrix_b1[2,2]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    fpr_b1 = conf_matrix_b1[2,1]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    x_b1 <- sensitivity(conf_matrix_b1)
    y_b1 <- 1-specificity(conf_matrix_b1)
    # Cut-off = 2
    predicted_values_b2 <- ifelse(bap_data_bootstrap$bap_score >= 2, 1, 0)
    conf_matrix_b2 <- table(factor(predicted_values_b2, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr_b2 = conf_matrix_b2[2,2]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    fpr_b2 = conf_matrix_b2[2,1]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    x_b2 <- sensitivity(conf_matrix_b2)
    y_b2 <- 1-specificity(conf_matrix_b2)
    # Cut-off = 3
    predicted_values_b3 <- ifelse(bap_data_bootstrap$bap_score >= 3, 1, 0)
    conf_matrix_b3 <- table(factor(predicted_values_b3, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr_b3 = conf_matrix_b3[2,2]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    fpr_b3 = conf_matrix_b3[2,1]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    x_b3 <- sensitivity(conf_matrix_b3)
    y_b3 <- 1-specificity(conf_matrix_b3)
    # Cut-off = 4
    predicted_values_b4 <- ifelse(bap_data_bootstrap$bap_score >= 4, 1, 0)
    conf_matrix_b4 <- table(factor(predicted_values_b4, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr_b4 = conf_matrix_b4[2,2]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    fpr_b4 = conf_matrix_b4[2,1]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    x_b4 <- sensitivity(conf_matrix_b4)
    y_b4 <- 1- specificity(conf_matrix_b4)
    # Cut-off = 5
    predicted_values_b5 <- ifelse(bap_data_bootstrap$bap_score >= 5, 1, 0)
    conf_matrix_b5 <- table(factor(predicted_values_b5, levels = 0:1),
                          factor(actual_values_b, levels = 0:1))
    tpr5 = conf_matrix_b5[2,2]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    fpr5 = conf_matrix_b5[2,1]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    x_b5 <- sensitivity(conf_matrix_b5)
    y_b5 <- 1-specificity(conf_matrix_b5)
    # AUC
    roc_bap_b <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5),
                          Sensitivity = c(x_b0, x_b1, x_b2, x_b3, x_b4, x_b5),
                          Specificity = c(y_b0, y_b1, y_b2, y_b3, y_b4, y_b5))
    roc_bap_list_b <- c(roc_bap_list_b, trapz(roc_bap_b$Specificity, roc_bap_b$Sensitivity)) 
  }
  res <- quantile(roc_bap_list_b,c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))
  auc_bap[[i]] <- res[2]
  se_bap[[i]] <- (res[2]-res[1])/1.96
}
# End
miinf_bap <- miInference(auc_bap, se_bap)
miinf_bap

#Est       SE Est/SE       df p Pct.mis
#50% 0.69361 0.030342  22.86 321349.9 0     1.8

## CURB-65
# Start
for (i in 1:dat_i$m) {
  # Data preparation
  com <- complete(dat_i, i)
  com <- com %>% 
    mutate_all(.funs = ~ as.numeric(.)) %>%  
    mutate(age = ifelse(age >= 65, 1, 0)) %>% 
    mutate(bun_curb = ifelse(bun >= 19, 1, 0)) %>% 
    mutate(rr = ifelse(rr >= 30, 1, 0)) %>% 
    mutate(bp = ifelse(sbp < 90 | dbp <= 60, 1, 0)) %>% 
    mutate(curb_score = con + bun_curb + rr + bp + age) 
  # Discrimination of CURB-65 
  # AUC
  # Bootstrap
  boot_num <- 2000
  roc_curb_list_b <- c()
  for(s in seq_len(boot_num)){
    set.seed(s)
    curb_data_bootstrap <-
      com[sample(dim(com)[1],dim(com)[1],replace=TRUE),]  
    # Cut-off = 0
    predicted_values_b0 <- ifelse(curb_data_bootstrap$curb_score >= 0, 1, 0)
    actual_values_b <- curb_data_bootstrap$death
    conf_matrix_b0 <- table(factor(predicted_values_b0, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b0 = conf_matrix_b0[2,2]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    fpr_b0 = conf_matrix_b0[2,1]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    x_b0 <- sensitivity(conf_matrix_b0)
    y_b0 <- 1-specificity(conf_matrix_b0)
    # Cut-off = 1
    predicted_values_b1 <- ifelse(curb_data_bootstrap$curb_score >= 1, 1, 0)
    conf_matrix_b1 <- table(factor(predicted_values_b1, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b1 = conf_matrix_b1[2,2]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    fpr_b1 = conf_matrix_b1[2,1]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    x_b1 <- sensitivity(conf_matrix_b1)
    y_b1 <- 1-specificity(conf_matrix_b1)
    # Cut-off = 2
    predicted_values_b2 <- ifelse(curb_data_bootstrap$curb_score >= 2, 1, 0)
    conf_matrix_b2 <- table(factor(predicted_values_b2, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b2 = conf_matrix_b2[2,2]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    fpr_b2 = conf_matrix_b2[2,1]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    x_b2 <- sensitivity(conf_matrix_b2)
    y_b2 <- 1-specificity(conf_matrix_b2)
    # Cut-off = 3
    predicted_values_b3 <- ifelse(curb_data_bootstrap$curb_score >= 3, 1, 0)
    conf_matrix_b3 <- table(factor(predicted_values_b3, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b3 = conf_matrix_b3[2,2]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    fpr_b3 = conf_matrix_b3[2,1]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    x_b3 <- sensitivity(conf_matrix_b3)
    y_b3 <- 1-specificity(conf_matrix_b3)
    # Cut-off = 4
    predicted_values_b4 <- ifelse(curb_data_bootstrap$curb_score >= 4, 1, 0)
    conf_matrix_b4 <- table(factor(predicted_values_b4, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b4 = conf_matrix_b4[2,2]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    fpr_b4 = conf_matrix_b4[2,1]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    x_b4 <- sensitivity(conf_matrix_b4)
    y_b4 <- 1- specificity(conf_matrix_b4)
    # Cut-off = 5
    predicted_values_b5 <- ifelse(curb_data_bootstrap$curb_score >= 5, 1, 0)
    conf_matrix_b5 <- table(factor(predicted_values_b5, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b5 = conf_matrix_b5[2,2]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    fpr_b5 = conf_matrix_b5[2,1]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    x_b5 <- sensitivity(conf_matrix_b5)
    y_b5 <- 1-specificity(conf_matrix_b5)
    # Cut-off = 6
    predicted_values_b6 <- ifelse(curb_data_bootstrap$curb_score >= 6, 1, 0)
    conf_matrix_b6 <- table(factor(predicted_values_b6, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b6 = conf_matrix_b6[2,2]/(conf_matrix_b6[2,1] + conf_matrix_b6[2,2])
    fpr_b6 = conf_matrix_b6[2,1]/(conf_matrix_b6[2,1] + conf_matrix_b6[2,2])
    x_b6 <- sensitivity(conf_matrix_b6)
    y_b6 <- 1-specificity(conf_matrix_b6)
    # AUC
    roc_curb_b <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5, 6),
                            Sensitivity = c(x_b0, x_b1, x_b2, x_b3, x_b4, x_b5, x_b6),
                            Specificity = c(y_b0, y_b1, y_b2, y_b3, y_b4, y_b5, y_b6))
    roc_curb_list_b <- c(roc_curb_list_b, trapz(roc_curb_b$Specificity, roc_curb_b$Sensitivity)) 
  }
  res <- quantile(roc_curb_list_b,c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))
  auc_curb[[i]] <- res[2]
  se_curb[[i]] <- (res[2]-res[1])/1.96
}
# End
miinf_curb <- miInference(auc_curb, se_curb)
miinf_curb
#Est       SE Est/SE     df p Pct.mis
#50% 0.68609 0.033115 20.718 130242 0     2.8

# XGBoost model -----------------------------------------------------------

original_data <- read.csv("Data/original_data.csv")
## Shuffle the original data
set.seed(19891116)
shuffled_data <- original_data[sample(nrow(original_data)),]
## Split sample
table(shuffled_data$death) %>% prop.table()
#0          1 
#0.92605042 0.07394958 
set.seed(19891116)
split_strat <- initial_split(shuffled_data, prop = 0.7,
                             strata = "death") # initial_split() split data to training/test data
train_data <- training(split_strat) # rsample::training makes training data
test_data <- testing(split_strat) # rsample::testing makes test data
table(train_data$death) %>% prop.table()
#0          1 
#0.92677071 0.07322929 
table(test_data$death) %>% prop.table()
#0          1 
#0.92436975 0.07563025 
## Converting data to suitable format for xgboost package
xgb_prep<- recipe(death ~ ., data = train_data) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = train_data, retain= TRUE) %>%
  juice () 
X <- as.matrix(xgb_prep[setdiff(names(xgb_prep), "death")]) # by setdiff(), extract columns except for "death"
Y <- xgb_prep$death
## Development
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
## Tuning max_depth and min_child weight
hyper_grid <- expand.grid(eta= 0.1,
                          max_depth= c(2,4,6,8,10),
                          min_child_weight= c(1,2,3,4,5),
                          subsample= 0.8,
                          colsample_bytree= 0.8,
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

#Rows: 25
#Columns: 10
#$ eta              <dbl> 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0~
#$ max_depth        <dbl> 4, 6, 10, 8, 4, 6, 4, 2, 2, 8, 10, 4, 8, 10, 6, 4, 2, 8, 10, 8, 10, 6, 6, 2, 2
#$ min_child_weight <dbl> 2, 2, 2, 2, 1, 1, 3, 2, 1, 1, 1, 4, 4, 4, 5, 5, 3, 3, 3, 5, 5, 3, 4, 5, 4
#$ subsample        <dbl> 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0~
#$ colsample_bytree <dbl> 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0~
#$ gamma            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#$ lambda           <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
#$ alpha            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#$ auc              <dbl> 0.8111747, 0.8110622, 0.8095992, 0.8071860, 0.8007230, 0.7997142, 0.7987108, 0.7984435, 0.798137~
#$ trees            <dbl> 35, 41, 40, 40, 38, 23, 32, 39, 40, 36, 23, 27, 22, 22, 29, 30, 36, 17, 17, 29, 29, 27, 27, 32, ~

## Cross-validation 
## Hyperparameters other than "nround" were fixed 
hyper_grid <- expand.grid(eta= 0.1,
                          max_depth= 4,
                          min_child_weight= 2,
                          subsample= 0.8,
                          colsample_bytree= 0.8,
                          gamma = 0,
                          lambda = 1,
                          alpha = 1,
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
hyper_grid %>% 
  arrange(desc(auc)) %>%
  glimpse ()
max(hyperparameters[["evaluation_log"]][["test_auc_mean"]])
#Rows: 1
#Columns: 10
#$ eta              <dbl> 0.1
#$ max_depth        <dbl> 4
#$ min_child_weight <dbl> 2
#$ subsample        <dbl> 0.8
#$ colsample_bytree <dbl> 0.8
#$ gamma            <dbl> 0
#$ lambda           <dbl> 1
#$ alpha            <dbl> 1
#$ auc              <dbl> 0.7903223
#$ trees            <dbl> 37

## Setup the optimal parameter list 
params <- list(
  eta= 0.1,
  max_depth= 4,
  min_child_weight= 2,
  subsample= 0.8,
  colsample_bytree= 0.8,
  gamma = 0,
  lambda = 1,
  alpha = 1) 
## Test in the train dataset 
xgb.train.final <- xgboost(
  params = params,
  data = X,
  label= Y,
  nrounds = 44,
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
# 0.974157 [ 0.95958  -  0.9887339 ]
train_predict <- ifelse(train_predict <= 0.5, -1, 1)
confusionMatrix(as.factor(train_predict), as.factor(ifelse(Y == 0, -1, 1)))
## Test in the test dataset
xgb_test<- recipe(death ~ ., data = test_data) %>% 
  step_integer(all_nominal ()) %>%
  prep(training = test_data, retain= TRUE) %>%
  juice () 
XX <- as.matrix(xgb_test[setdiff(names(xgb_test), "death")]) 
YY <- xgb_test$death
test_predict <- predict(xgb.train.final, XX)
## Performance
AUC(test_predict, YY)
cstatNo <- rcorr.cens(test_predict, YY) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
# 0.6816498 [ 0.5762702  -  0.7870295 ]
## Visualizing the tuning
plot(hyper_grid$auc, 
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
  subsample= 0.8,
  colsample_bytree= 0.8,
  gamma = 0,
  lambda = 1,
  alpha = 0) 
## Test in the train dataset 
xgb.train.final <- xgboost(
  params = params,
  data = X,
  label= Y,
  nrounds = 39,
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
test_predict <- predict(xgb.train.final, XX)
## Performance
AUC(test_predict, YY)
cstatNo <- rcorr.cens(test_predict, YY) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
## Plot importance of variables 
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
                                    model = xgb.train.final)
xgb.plot.importance(importance_matrix[1:10, ])

# Model comparison --------------------------------------------------------

# ROC curve
test_data_roc <- test_data %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>%  
  mutate(age_cat = ifelse(age >= 65, 1, 0)) %>% 
  mutate(bun_curb = ifelse(bun >= 25, 1, 0)) %>%
  mutate(bun_bap = ifelse(bun >= 19, 1, 0)) %>% 
  mutate(rr_cat = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp_cat = ifelse(sbp < 90 | dbp <= 60, 1, 0)) %>%
  mutate(hr_cat = ifelse(hr >= 109, 1, 0))
test_data_roc <-test_data_roc %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(curb_score = con + bun_curb + rr_cat + bp_cat + age_cat) %>% 
  mutate(bap_score = bun_bap + con + hr_cat + age_cat) 
## BAP-65
# Cut-off = 0
predicted_values_c0 <- ifelse(test_data_roc$bap_score >= 0, 1, 0)
actual_values_c <- test_data_roc$death
conf_matrix_c0 <- table(factor(predicted_values_c0, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr0 = conf_matrix_c0[2,2]/(conf_matrix_c0[2,1] + conf_matrix_c0[2,2])
fpr0 = conf_matrix_c0[2,1]/(conf_matrix_c0[2,1] + conf_matrix_c0[2,2])
x_c0 <- sensitivity(conf_matrix_c0)
y_c0 <- 1-specificity(conf_matrix_c0)
# Cut-off = 1
predicted_values_c1 <- ifelse(test_data_roc$bap_score >= 1, 1, 0)
conf_matrix_c1 <- table(factor(predicted_values_c1, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr_c1 = conf_matrix_c1[2,2]/(conf_matrix_c1[2,1] + conf_matrix_c1[2,2])
fpr_c1 = conf_matrix_c1[2,1]/(conf_matrix_c1[2,1] + conf_matrix_c1[2,2])
x_c1 <- sensitivity(conf_matrix_c1)
y_c1 <- 1-specificity(conf_matrix_c1)
# Cut-off = 2
predicted_values_c2 <- ifelse(test_data_roc$bap_score >= 2, 1, 0)
conf_matrix_c2 <- table(factor(predicted_values_c2, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr_c2 = conf_matrix_c2[2,2]/(conf_matrix_c2[2,1] + conf_matrix_c2[2,2])
fpr_c2 = conf_matrix_c2[2,1]/(conf_matrix_c2[2,1] + conf_matrix_c2[2,2])
x_c2 <- sensitivity(conf_matrix_c2)
y_c2 <- 1-specificity(conf_matrix_c2)
# Cut-off = 3
predicted_values_c3 <- ifelse(test_data_roc$bap_score >= 3, 1, 0)
conf_matrix_c3 <- table(factor(predicted_values_c3, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr_c3 = conf_matrix_c3[2,2]/(conf_matrix_c3[2,1] + conf_matrix_c3[2,2])
fpr_c3 = conf_matrix_c3[2,1]/(conf_matrix_c3[2,1] + conf_matrix_c3[2,2])
x_c3 <- sensitivity(conf_matrix_c3)
y_c3 <- 1-specificity(conf_matrix_c3)
# Cut-off = 4
predicted_values_c4 <- ifelse(test_data_roc$bap_score >= 4, 1, 0)
conf_matrix_c4 <- table(factor(predicted_values_c4, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr_c4 = conf_matrix_c4[2,2]/(conf_matrix_c4[2,1] + conf_matrix_c4[2,2])
fpr_c4 = conf_matrix_c4[2,1]/(conf_matrix_c4[2,1] + conf_matrix_c4[2,2])
x_c4 <- sensitivity(conf_matrix_c4)
y_c4 <- 1- specificity(conf_matrix_c4)
# Cut-off = 5
predicted_values_c5 <- ifelse(test_data_roc$bap_score >= 5, 1, 0)
conf_matrix_c5 <- table(factor(predicted_values_c5, levels = 0:1),
                      factor(actual_values_c, levels = 0:1))
tpr_c5 = conf_matrix_c5[2,2]/(conf_matrix_c5[2,1] + conf_matrix_c5[2,2])
fpr_c5 = conf_matrix_c5[2,1]/(conf_matrix_c5[2,1] + conf_matrix_c5[2,2])
x_c5 <- sensitivity(conf_matrix_c5)
y_c5 <- 1-specificity(conf_matrix_c5)
# Describing ROC curve
roc_bap_c <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5),
                      Sensitivity = c(x_c0, x_c1, x_c2, x_c3, x_c4, x_c5),
                      Specificity = c(y_c0, y_c1, y_c2, y_c3, y_c4, y_c5))
trapz(roc_bap_c$Specificity, roc_bap_c$Sensitivity) # Trapezoid integration
plot(Sensitivity ~ Specificity,
     data = roc_bap_c,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "1-Specificity",
     ylab = "Sensitivity",
     type="l",
     lwd = 2)
## CURB-65
# Cut-off = 0
predicted_values_c00 <- ifelse(test_data_roc$curb_score >= 0, 1, 0)
conf_matrix_c00 <- table(factor(predicted_values_c00, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c00 = conf_matrix_c00[2,2]/(conf_matrix_c00[2,1] + conf_matrix_c00[2,2])
fpr_c00 = conf_matrix_c00[2,1]/(conf_matrix_c00[2,1] + conf_matrix_c00[2,2])
x_c00 <- sensitivity(conf_matrix_c00)
y_c00 <- 1-specificity(conf_matrix_c00)
# Cut-off = 1
predicted_values_c01 <- ifelse(test_data_roc$curb_score >= 1, 1, 0)
conf_matrix_c01 <- table(factor(predicted_values_c01, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c01 = conf_matrix_c01[2,2]/(conf_matrix_c01[2,1] + conf_matrix_c01[2,2])
fpr_c01 = conf_matrix_c01[2,1]/(conf_matrix_c01[2,1] + conf_matrix_c01[2,2])
x_c01 <- sensitivity(conf_matrix_c01)
y_c01 <- 1-specificity(conf_matrix_c01)
# Cut-off = 2
predicted_values_c02 <- ifelse(test_data_roc$curb_score >= 2, 1, 0)
conf_matrix_c02 <- table(factor(predicted_values_c02, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c02 = conf_matrix_c02[2,2]/(conf_matrix_c02[2,1] + conf_matrix_c02[2,2])
fpr_c02 = conf_matrix_c02[2,1]/(conf_matrix_c02[2,1] + conf_matrix_c02[2,2])
x_c02 <- sensitivity(conf_matrix_c02)
y_c02 <- 1-specificity(conf_matrix_c02)
# Cut-off = 3
predicted_values_c03 <- ifelse(test_data_roc$curb_score >= 3, 1, 0)
conf_matrix_c03 <- table(factor(predicted_values_c03, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c03 = conf_matrix_c03[2,2]/(conf_matrix_c03[2,1] + conf_matrix_c03[2,2])
fpr_c03 = conf_matrix_c03[2,1]/(conf_matrix_c03[2,1] + conf_matrix_c03[2,2])
x_c03 <- sensitivity(conf_matrix_c03)
y_c03 <- 1-specificity(conf_matrix_c03)
# Cut-off = 4
predicted_values_c04 <- ifelse(test_data_roc$curb_score >= 4, 1, 0)
conf_matrix_c04 <- table(factor(predicted_values_c04, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c04 = conf_matrix_c04[2,2]/(conf_matrix_c04[2,1] + conf_matrix_c04[2,2])
fpr_c04 = conf_matrix_c04[2,1]/(conf_matrix_c04[2,1] + conf_matrix_c04[2,2])
x_c04 <- sensitivity(conf_matrix_c04)
y_c04 <- 1- specificity(conf_matrix_c04)
# Cut-off = 5
predicted_values_c05 <- ifelse(test_data_roc$curb_score >= 5, 1, 0)
conf_matrix_c05 <- table(factor(predicted_values_c05, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c05 = conf_matrix_c05[2,2]/(conf_matrix_c05[2,1] + conf_matrix_c05[2,2])
fpr_c05 = conf_matrix_c05[2,1]/(conf_matrix_c05[2,1] + conf_matrix_c05[2,2])
x_c05 <- sensitivity(conf_matrix_c05)
y_c05 <- 1-specificity(conf_matrix_c05)
# Cut-off = 6
predicted_values_c06 <- ifelse(test_data_roc$curb_score >= 6, 1, 0)
conf_matrix_c06 <- table(factor(predicted_values_c06, levels = 0:1),
                       factor(actual_values_c, levels = 0:1))
tpr_c06 = conf_matrix_c06[2,2]/(conf_matrix_c06[2,1] + conf_matrix_c06[2,2])
fpr_c06 = conf_matrix_c06[2,1]/(conf_matrix_c06[2,1] + conf_matrix_c06[2,2])
y_c06 <- sensitivity(conf_matrix_c06)
x_c06 <- 1- specificity(conf_matrix_c06)
# Describing ROC curve
roc_curb_c <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5, 6),
                       Sensitivity = c(x_c00, x_c01, x_c02, x_c03, x_c04, x_c05, x_c06),
                       Specificity = c(y_c00, y_c01, y_c02, y_c03, y_c04, y_c05, y_c06))
trapz(roc_curb_c$Specificity, roc_curb_c$Sensitivity) # Trapezoid integration
par(new=T)
plot(Sensitivity ~ Specificity,
     data = roc_curb_c,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "1-Specificity",
     ylab = "Sensitivity",
     lty = "dashed",
     lwd = 2,
     type="l")
## XGBoost model
roc3 <- roc(YY, test_predict)
## ROC total
specificity_inv <- 1- roc3[["specificities"]]
par(new=T)
plot(roc3[["sensitivities"]] ~ specificity_inv,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "1-Specificity",
     ylab = "Sensitivity",
     lwd = 2,
     lty = "dotted",
     type="l")
abline(0, 1, col = "grey")
legend(x = 0.8, y = 0.3, bty = "n", lwd = 2, lty = c("solid","dashed","dotted"),
       legend = c("BAP-65","CURB-65","GXBoost"))
## AUC comparison
boot_num <- 2000 
auc_bap_comp <- auc_curb_comp <- auc_xgb_comp <- vector(length = dat_i$m, mode = "list")
set.seed(198911116)
for(s in seq_len(boot_num)){
    data_bootstrap <-
      test_data_roc[sample(dim(test_data_roc)[1],dim(test_data_roc)[1],replace=TRUE),]  
    # AUC 0f BAP-65
    ## Cut-off = 0
    predicted_values_b0 <- ifelse(data_bootstrap$bap_score >= 0, 1, 0)
    actual_values_b <- data_bootstrap$death
    conf_matrix_b0 <- table(factor(predicted_values_b0, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b0 = conf_matrix_b0[2,2]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    fpr_b0 = conf_matrix_b0[2,1]/(conf_matrix_b0[2,1] + conf_matrix_b0[2,2])
    x_b0 <- sensitivity(conf_matrix_b0)
    y_b0 <- 1-specificity(conf_matrix_b0)
    ## Cut-off = 1
    predicted_values_b1 <- ifelse(data_bootstrap$bap_score >= 1, 1, 0)
    conf_matrix_b1 <- table(factor(predicted_values_b1, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b1 = conf_matrix_b1[2,2]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    fpr_b1 = conf_matrix_b1[2,1]/(conf_matrix_b1[2,1] + conf_matrix_b1[2,2])
    x_b1 <- sensitivity(conf_matrix_b1)
    y_b1 <- 1-specificity(conf_matrix_b1)
    ## Cut-off = 2
    predicted_values_b2 <- ifelse(data_bootstrap$bap_score >= 2, 1, 0)
    conf_matrix_b2 <- table(factor(predicted_values_b2, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b2 = conf_matrix_b2[2,2]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    fpr_b2 = conf_matrix_b2[2,1]/(conf_matrix_b2[2,1] + conf_matrix_b2[2,2])
    x_b2 <- sensitivity(conf_matrix_b2)
    y_b2 <- 1-specificity(conf_matrix_b2)
    ## Cut-off = 3
    predicted_values_b3 <- ifelse(data_bootstrap$bap_score >= 3, 1, 0)
    conf_matrix_b3 <- table(factor(predicted_values_b3, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b3 = conf_matrix_b3[2,2]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    fpr_b3 = conf_matrix_b3[2,1]/(conf_matrix_b3[2,1] + conf_matrix_b3[2,2])
    x_b3 <- sensitivity(conf_matrix_b3)
    y_b3 <- 1-specificity(conf_matrix_b3)
    ## Cut-off = 4
    predicted_values_b4 <- ifelse(data_bootstrap$bap_score >= 4, 1, 0)
    conf_matrix_b4 <- table(factor(predicted_values_b4, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b4 = conf_matrix_b4[2,2]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    fpr_b4 = conf_matrix_b4[2,1]/(conf_matrix_b4[2,1] + conf_matrix_b4[2,2])
    x_b4 <- sensitivity(conf_matrix_b4)
    y_b4 <- 1- specificity(conf_matrix_b4)
    ## Cut-off = 5
    predicted_values_b5 <- ifelse(data_bootstrap$bap_score >= 5, 1, 0)
    conf_matrix_b5 <- table(factor(predicted_values_b5, levels = 0:1),
                            factor(actual_values_b, levels = 0:1))
    tpr_b5 = conf_matrix_b5[2,2]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    fpr_b5 = conf_matrix_b5[2,1]/(conf_matrix_b5[2,1] + conf_matrix_b5[2,2])
    x_b5 <- sensitivity(conf_matrix_b5)
    y_b5 <- 1-specificity(conf_matrix_b5)
    # AUC 
    roc_bap_boot <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5),
                               Sensitivity = c(x_b0, x_b1, x_b2, x_b3, x_b4, x_b5),
                               Specificity = c(y_b0, y_b1, y_b2, y_b3, y_b4, y_b5))
    auc_bap_comp <- c(auc_bap_comp, trapz(roc_bap_boot$Specificity, roc_bap_boot$Sensitivity)) 
    # AUC of CURB-65
    # Cut-off = 0
    predicted_values_bb0 <- ifelse(data_bootstrap$curb_score >= 0, 1, 0)
    actual_values_bb <- data_bootstrap$death
    conf_matrix_bb0 <- table(factor(predicted_values_bb0, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb0 = conf_matrix_bb0[2,2]/(conf_matrix_bb0[2,1] + conf_matrix_bb0[2,2])
    fpr_bb0 = conf_matrix_bb0[2,1]/(conf_matrix_bb0[2,1] + conf_matrix_bb0[2,2])
    x_bb0 <- sensitivity(conf_matrix_bb0)
    y_bb0 <- 1-specificity(conf_matrix_bb0)
    # Cut-off = 1
    predicted_values_bb1 <- ifelse(data_bootstrap$curb_score >= 1, 1, 0)
    conf_matrix_bb1 <- table(factor(predicted_values_bb1, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb1 = conf_matrix_bb1[2,2]/(conf_matrix_bb1[2,1] + conf_matrix_bb1[2,2])
    fpr_bb1 = conf_matrix_bb1[2,1]/(conf_matrix_bb1[2,1] + conf_matrix_bb1[2,2])
    x_bb1 <- sensitivity(conf_matrix_bb1)
    y_bb1 <- 1-specificity(conf_matrix_bb1)
    # Cut-off = 2
    predicted_values_bb2 <- ifelse(data_bootstrap$curb_score >= 2, 1, 0)
    conf_matrix_bb2 <- table(factor(predicted_values_bb2, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb2 = conf_matrix_bb2[2,2]/(conf_matrix_bb2[2,1] + conf_matrix_bb2[2,2])
    fpr_bb2 = conf_matrix_bb2[2,1]/(conf_matrix_bb2[2,1] + conf_matrix_bb2[2,2])
    x_bb2 <- sensitivity(conf_matrix_bb2)
    y_bb2 <- 1-specificity(conf_matrix_bb2)
    # Cut-off = 3
    predicted_values_bb3 <- ifelse(data_bootstrap$curb_score >= 3, 1, 0)
    conf_matrix_bb3 <- table(factor(predicted_values_bb3, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb3 = conf_matrix_bb3[2,2]/(conf_matrix_bb3[2,1] + conf_matrix_bb3[2,2])
    fpr_bb3 = conf_matrix_bb3[2,1]/(conf_matrix_bb3[2,1] + conf_matrix_bb3[2,2])
    x_bb3 <- sensitivity(conf_matrix_bb3)
    y_bb3 <- 1-specificity(conf_matrix_bb3)
    # Cut-off = 4
    predicted_values_bb4 <- ifelse(data_bootstrap$curb_score >= 4, 1, 0)
    conf_matrix_bb4 <- table(factor(predicted_values_bb4, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb4 = conf_matrix_bb4[2,2]/(conf_matrix_bb4[2,1] + conf_matrix_bb4[2,2])
    fpr_bb4 = conf_matrix_bb4[2,1]/(conf_matrix_bb4[2,1] + conf_matrix_bb4[2,2])
    x_bb4 <- sensitivity(conf_matrix_bb4)
    y_bb4 <- 1- specificity(conf_matrix_bb4)
    # Cut-off = 5
    predicted_values_bb5 <- ifelse(data_bootstrap$curb_score >= 5, 1, 0)
    conf_matrix_bb5 <- table(factor(predicted_values_bb5, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb5 = conf_matrix_bb5[2,2]/(conf_matrix_bb5[2,1] + conf_matrix_bb5[2,2])
    fpr_bb5 = conf_matrix_bb5[2,1]/(conf_matrix_bb5[2,1] + conf_matrix_bb5[2,2])
    x_bb5 <- sensitivity(conf_matrix_bb5)
    y_bb5 <- 1-specificity(conf_matrix_bb5)
    # Cut-off = 6
    predicted_values_bb6 <- ifelse(data_bootstrap$curb_score >= 6, 1, 0)
    conf_matrix_bb6 <- table(factor(predicted_values_bb6, levels = 0:1),
                             factor(actual_values_bb, levels = 0:1))
    tpr_bb6 = conf_matrix_bb6[2,2]/(conf_matrix_bb6[2,1] + conf_matrix_bb6[2,2])
    fpr_bb6 = conf_matrix_bb6[2,1]/(conf_matrix_bb6[2,1] + conf_matrix_bb6[2,2])
    x_bb6 <- sensitivity(conf_matrix_bb6)
    y_bb6 <- 1-specificity(conf_matrix_bb6)
    # AUC
    roc_curb_boot <- data.frame(cutpoint = c(0, 1, 2, 3, 4, 5, 6),
                                Sensitivity = c(x_bb0, x_bb1, x_bb2, x_bb3, x_bb4, x_bb5, x_bb6),
                                Specificity = c(y_bb0, y_bb1, y_bb2, y_bb3, y_bb4, y_bb5, y_bb6))
    auc_curb_comp <- c(auc_curb_comp, trapz(roc_curb_boot$Specificity, roc_curb_boot$Sensitivity))  
    # AUC of XGBoost
    data_bootstrap_xgb <- data_bootstrap %>% 
      select(sex, con, bun, eo, adl, rr, sbp, dbp, hr, age, death)
    boot_prep<- 
      recipe(death ~ ., data = data_bootstrap_xgb) %>% 
      step_integer(all_nominal ()) %>%
      prep(training = data_bootstrap_xgb, retain= TRUE) %>%
      juice () 
    boot_xgb_X <- as.matrix(boot_prep[setdiff(names(boot_prep), "death")]) 
    boot_xgb_Y <- data_bootstrap$death
    pred_boot_xgb <- predict(xgb.train.final,boot_xgb_X)
    roc_xgb <- roc(boot_xgb_Y, pred_boot_xgb)
    auc_xgb_comp <- c(auc_xgb_comp, roc_xgb$auc)
    } # End
bap_vec <- unlist(auc_bap_comp, use.names=FALSE)
curb_vec <- unlist(auc_curb_comp, use.names=FALSE)
xgb_vec <- unlist(auc_xgb_comp, use.names=FALSE)
## BAP-65 vs XGBoost
diff_bap_xgb <- bap_vec - xgb_vec
quantile(diff_bap_xgb, c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))
## CURB-65 vs XGBoost
diff_curb_xgb <- curb_vec - xgb_vec
quantile(diff_curb_xgb, c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))
diff_bap_curb <- bap_vec - curb_vec
quantile(diff_bap_curb, c(0+(1-0.95)/2, .5, 1-(1-0.95)/2))