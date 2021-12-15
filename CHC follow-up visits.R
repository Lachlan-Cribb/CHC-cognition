# load packages

library(haven)
library(tidyverse)
library(lavaan)
library(psych)
library(semTools)
library(mice)

#### Baseline CFA data ####

## Load data 

arcli <- read_sav("/home/lachlan/Dropbox/ARCLI data entry/BLINDED DATA SETS/BLINDED _DO NOT TOUCH/ARCLI_V1_Baseline_ALL_7-SEPT-2021.sav")

# remove SPSS missing data indicators

arcli <- arcli %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., 999)))

# remove blank rows

arcli <- arcli %>% 
  filter(!IDno == "")

# Rename problematic variable names

arcli <- arcli %>% 
  rename(CDR_ImmWRCorrect_v1 = `CDR_ImmWRCorrect#_v1`) %>% 
  rename(Serial3Correct_v1 = `Serial3sCorrectresponses#_v1`) %>% 
  rename(Serial3Total_v1 = `Serial3sTotalResponses#_v1`) %>% 
  rename(Serial7Correct_v1 = `Serial7sCorrectresponses#_v1`) %>% 
  rename(Serial7Total_v1 = `Serial7sTotalResponses#_v1`)

arcli <- arcli %>% 
  rename(CDR_DelayedWRCorrect_v1 = `CDR_DelayedWRCorrect#_v1`)

# create new overall score

arcli$CDR_NumericWMOverallAcc_v1 <- (arcli$CDR_NumericWMOrigStimAcc_v1 + arcli$CDR_NumericWMNewStimAcc_v1) / 2

# set sex variable to factor

arcli$Sex <- as.factor(arcli$Sex)

levels(arcli$Sex) <- c("male", "female")

# select relevant variables 

cog_vars <- arcli %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v1, Initialage,CDR_RVIPacc_v1,CDR_ImmWRCorrect_v1, CDR_DelayedWRCorrect_v1, 
                CDR_SimpleRTms_v1, CDR_ChoiceRTms_v1, CDR_NumericWMOverallAcc_v1, 
                CDR_DigitVigRTms_v1, SCB_ConStroopRTms_v1, SCB_InconStroopRTms_v1, 
                Jenson8ChoiceDecisionTimems_v1,Jenson4ChoiceDecisionTimems_v1,
                Jenson1ChoiceDecisionTimems_v1,Serial3sAcc_v1, Serial7sAcc_v1)


## Remove rows without cog results

cog_vars <- cog_vars[which(rowMeans(is.na(cog_vars)) < 0.80),]

# log transform RT variables

RT_vars <- cog_vars %>% select(
  CDR_RVIPRTms_v1,CDR_SimpleRTms_v1, CDR_ChoiceRTms_v1, 
  CDR_DigitVigRTms_v1, SCB_ConStroopRTms_v1, SCB_InconStroopRTms_v1, 
  Jenson8ChoiceDecisionTimems_v1,Jenson4ChoiceDecisionTimems_v1,
  Jenson1ChoiceDecisionTimems_v1) %>% 
  names() # variable names

cog_vars <- cog_vars %>% 
  mutate(across(all_of(RT_vars), ~ log2(.), .names = "log{.col}")) # log 2 transform

# check histograms #

cog_vars[,c(19:27)] %>% select(where(is.numeric)) %>% 
  multi.hist(global = F)

cog_vars %>% select(CDR_ImmWRCorrect_v1, CDR_DelayedWRCorrect_v1,
                               CDR_NumericWMOverallAcc_v1,
                               Serial3sAcc_v1, Serial7sAcc_v1) %>% 
  multi.hist(global = F)


# winsorise outliers for RT vars

winsorise <- function(x, sd){
  maxin <- max(x[scale(x) < sd],na.rm=T)
  minin <- min(x[scale(x) > -sd],na.rm=T)
  x <- ifelse(x > maxin & !is.na(x), maxin, x)
  x <- ifelse(x < minin & !is.na(x), minin, x)
  x
}

# Winsorise outliers 

cog_vars <- cog_vars %>% # winsorise everything except age
  mutate(across(c(where(is.numeric), -Initialage), winsorise, sd = 3, .names = "w{.col}"))

cog_vars2 <- cog_vars %>% select(IDno, Initialage, Sex, starts_with("wlog"), 
                                 wCDR_ImmWRCorrect_v1,wCDR_DelayedWRCorrect_v1, 
                                 wCDR_NumericWMOverallAcc_v1, wSerial3sAcc_v1,
                                 wSerial7sAcc_v1) # select relevant variables

# histograms of winsorised variables

cog_vars2 %>% select(starts_with("wlog")) %>% 
  multi.hist(global = F)

cog_vars2 %>% select(wCDR_ImmWRCorrect_v1, wCDR_DelayedWRCorrect_v1,
                    wCDR_NumericWMOverallAcc_v1,wSerial3sAcc_v1, wSerial7sAcc_v1) %>% 
  multi.hist(global = F)

## remove cases with more than half cog vars missing

nrow(cog_vars2)

cog_vars2 <- cog_vars2[which(rowMeans(is.na(cog_vars2)) < 0.50),]

nrow(cog_vars2)

# check missing data per variable

n_miss <- colSums(is.na(cog_vars2))

median(n_miss / nrow(cog_vars2))

sort(n_miss / nrow(cog_vars2))

# Reverse score reaction time measures

reverse_RT <- function(x){
  x <- x * -1
  x <- x + max(abs(x),na.rm=T) + 1
  x
}

cog_vars2 <- cog_vars2 %>% 
  mutate(across(starts_with("wlog"), reverse_RT, .names = "{.col}_rev"))

## rescale accuracy variables (10% units)

cog_vars2 <- cog_vars2 %>% mutate(across(starts_with("wSerial"), ~  . / 10))

cog_vars2 <- cog_vars2 %>% mutate(wCDR_NumericWMOverallAcc_v1 = wCDR_NumericWMOverallAcc_v1 / 10)

# select relevant variables #

cog_vars3 <- cog_vars2 %>% 
  select(IDno, Sex, Initialage, ends_with("_rev"), wCDR_ImmWRCorrect_v1,
         wCDR_DelayedWRCorrect_v1, wCDR_NumericWMOverallAcc_v1,
         wSerial3sAcc_v1, wSerial7sAcc_v1)

# drop v1 from variable names

names(cog_vars3) <- sub("_v1", "", names(cog_vars3))

## Multiple imputation using predictive mean matching

predictormatrix <- quickpred(cog_vars3, 
                           include=c("Sex", "Initialage"),
                           exclude="IDno",
                           mincor = 0.1)

cog_vars3_imp <- mice(cog_vars3, m = 20, method = "pmm",
                      predictorMatrix = predictormatrix, diagnostics = T,
                      seed = 10)

imp_dat <- complete(cog_vars3_imp, include = F, 
                    action = "all")

# check imputation for each variable 

stripplot(cog_vars3_imp, wlogCDR_SimpleRTms_v1_rev +
                  wlogCDR_ChoiceRTms_v1_rev +
                  wlogCDR_DigitVigRTms_v1_rev +
                  wlogJenson8ChoiceDecisionTimems_v1_rev +
                  wlogJenson4ChoiceDecisionTimems_v1_rev +
                  wlogJenson1ChoiceDecisionTimems_v1_rev ~ .imp)

stripplot(cog_vars3_imp, wlogSCB_ConStroopRTms_v1_rev +
            wlogSCB_InconStroopRTms_v1_rev  ~ .imp)


stripplot(cog_vars3_imp, wCDR_ImmWRCorrect_v1 +
            wCDR_DelayedWRCorrect_v1 +
            wSerial3sAcc_v1 +
            wSerial7sAcc_v1 +
            wCDR_NumericWMOverallAcc_v1  ~ .imp)


#### Fit CFA model using baseline data ####

# define CFA model 

model_bl <- 
  
' # model 
  
## Factor loadings 
  
Gt =~ wlogCDR_SimpleRTms_rev +
       wlogCDR_SimpleRTms_rev +
       wlogCDR_ChoiceRTms_rev +
       wlogCDR_DigitVigRTms_rev +
       wlogJenson8ChoiceDecisionTimems_rev +
       wlogJenson4ChoiceDecisionTimems_rev +
       wlogJenson1ChoiceDecisionTimems_rev
       
Gs =~  wlogSCB_ConStroopRTms_rev +
       wlogSCB_ConStroopRTms_rev +
       wlogSCB_InconStroopRTms_rev +
       # cross loading
       wlogCDR_ChoiceRTms_rev 

       
Glr_m6 =~ wCDR_ImmWRCorrect +
       wCDR_ImmWRCorrect +
       wCDR_DelayedWRCorrect 
       

Gwm_ac =~ wSerial3sAcc +
        wSerial3sAcc +
        wSerial7sAcc +
        wCDR_NumericWMOverallAcc
        

## correlated residuals

wlogJenson4ChoiceDecisionTimems_rev ~~ wlogJenson1ChoiceDecisionTimems_rev
wlogJenson8ChoiceDecisionTimems_rev ~~ wlogJenson1ChoiceDecisionTimems_rev
wlogJenson4ChoiceDecisionTimems_rev ~~ wlogJenson8ChoiceDecisionTimems_rev
'

# fit model

test_bl <- cfa.mi(model_bl, data = imp_dat, estimator = "MLR", std.lv = F)

# results

summary(test_bl, fit.measures=T, standardized = T, asymptotic = T)

# Save CHC factors to data frame
### PROBLEM ###
## NO predict method for extracting factor scores from lavaan.mi objects ##




#### 3 month CFA data ####

v2 <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\BLINDED DATA SETS\\BLINDED _DO NOT TOUCH\\ARCLI_V2_3 Month_ALL_6-SEPT-2021.sav")

# remove SPSS missing data indicators

v2 <- v2 %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., 999)))

# remove blank rows

v2 <- v2 %>% 
  filter(!IDno == "")

# Rename problematic variable names

v2 <- v2 %>% 
  rename(CDR_ImmWRCorrect_v2 = `CDR_ImmWRCorrect#_v2`)

v2 <- v2 %>% 
  rename(CDR_DelayedWRCorrect_v2 = `CDR_DelayedWRCorrect#_v2`)

# create new overall score

v2$CDR_NumericWMOverallAcc_v2 <- (v2$CDR_NumericWMOrigStimAcc_v2 + v2$CDR_NumericWMNewStimAcc_v2) / 2

# set sex variable to factor

v2$Sex <- as.factor(v2$Sex)

levels(v2$Sex) <- c("male", "female")

# Add age variable to v2 dataset 

age_bl <- arcli %>% select(IDno, Initialage)

v2 <- v2 %>% 
  left_join(age_bl, by = "IDno")

# select relevant variables 

cog_vars_v2 <- v2 %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v2, Initialage,CDR_RVIPacc_v2,CDR_ImmWRCorrect_v2, CDR_DelayedWRCorrect_v2, 
                CDR_SimpleRTms_v2, CDR_ChoiceRTms_v2, CDR_NumericWMOverallAcc_v2, 
                CDR_DigitVigRTms_v2, SCB_ConStroopRTms_v2, SCB_InconStroopRTms_v2, 
                Jenson8ChoiceDecisionTimems_v2,Jenson4ChoiceDecisionTimems_v2,
                Jenson1ChoiceDecisionTimems_v2,Serial3sAcc_v2, Serial7sAcc_v2)


## Remove rows without cog results

cog_vars_v2 <- cog_vars_v2[which(rowMeans(is.na(cog_vars_v2)) < 0.80),]

# log transform RT variables

RT_vars <- cog_vars_v2 %>% select(
  CDR_RVIPRTms_v2,CDR_SimpleRTms_v2, CDR_ChoiceRTms_v2, 
  CDR_DigitVigRTms_v2, SCB_ConStroopRTms_v2, SCB_InconStroopRTms_v2, 
  Jenson8ChoiceDecisionTimems_v2,Jenson4ChoiceDecisionTimems_v2,
  Jenson1ChoiceDecisionTimems_v2) %>% 
  names() # variable names

cog_vars_v2 <- cog_vars_v2 %>% 
  mutate(across(all_of(RT_vars), ~ log2(.), .names = "log{.col}")) # log 2 transform

# check histograms #

cog_vars_v2[,c(19:27)] %>% select(where(is.numeric)) %>% 
  multi.hist(global = F)

cog_vars_v2 %>% select(CDR_ImmWRCorrect_v2, CDR_DelayedWRCorrect_v2,
                    CDR_NumericWMOverallAcc_v2,
                    Serial3sAcc_v2, Serial7sAcc_v2) %>% 
  multi.hist(global = F)

# Winsorise outliers 

cog_vars_v2 <- cog_vars_v2 %>% # winsorise everything except age
  mutate(across(c(where(is.numeric), -Initialage), winsorise, sd = 3, .names = "w{.col}"))

cog_vars_v22 <- cog_vars_v2 %>% select(IDno, Initialage, Sex, starts_with("wlog"), 
                                 wCDR_ImmWRCorrect_v2,wCDR_DelayedWRCorrect_v2, 
                                 wCDR_NumericWMOverallAcc_v2, wSerial3sAcc_v2,
                                 wSerial7sAcc_v2) # select relevant variables

# histograms of winsorised variables

cog_vars_v22 %>% select(starts_with("wlog")) %>% 
  multi.hist(global = F)

cog_vars_v22 %>% select(wCDR_ImmWRCorrect_v2, wCDR_DelayedWRCorrect_v2,
                     wCDR_NumericWMOverallAcc_v2,wSerial3sAcc_v2, wSerial7sAcc_v2) %>% 
  multi.hist(global = F)

## remove cases with more than half cog vars missing

nrow(cog_vars_v22)

cog_vars_v22 <- cog_vars_v22[which(rowMeans(is.na(cog_vars_v22)) < 0.50),]

nrow(cog_vars_v22)

# check missing data per variable

n_miss <- colSums(is.na(cog_vars_v22))

median(n_miss / nrow(cog_vars_v22))

sort(n_miss / nrow(cog_vars_v22))

# Reverse score reaction time measures

cog_vars_v22 <- cog_vars_v22 %>% 
  mutate(across(starts_with("wlog"), reverse_RT, .names = "{.col}_rev"))

## rescale accuracy variables (10% units)

cog_vars_v22 <- cog_vars_v22 %>% mutate(across(starts_with("wSerial"), ~  . / 10))

cog_vars_v22 <- cog_vars_v22 %>% mutate(wCDR_NumericWMOverallAcc_v2 = wCDR_NumericWMOverallAcc_v2 / 10)

# select relevant variables #

cog_vars_v23 <- cog_vars_v22 %>% 
  select(IDno, Sex, Initialage, ends_with("_rev"), wCDR_ImmWRCorrect_v2,
         wCDR_DelayedWRCorrect_v2, wCDR_NumericWMOverallAcc_v2,
         wSerial3sAcc_v2, wSerial7sAcc_v2)

## single imputation

# predictive mean matching

predictormatrix <- quickpred(cog_vars_v23, 
                             include=c("Sex", "Initialage"),
                             exclude="IDno",
                             mincor = 0.1)

cog_vars_v23_imp <- mice(cog_vars_v23, m = 1, method = "pmm",
                      predictorMatrix = predictormatrix, diagnostics = T,
                      seed = 10)

imp_dat_v2 <- complete(cog_vars_v23_imp) # imputed data 

# check imputation for each variable 

stripplot(cog_vars_v23_imp, wlogCDR_SimpleRTms_v2_rev +
            wlogCDR_ChoiceRTms_v2_rev +
            wlogCDR_DigitVigRTms_v2_rev +
            wlogJenson8ChoiceDecisionTimems_v2_rev +
            wlogJenson4ChoiceDecisionTimems_v2_rev +
            wlogJenson1ChoiceDecisionTimems_v2_rev ~ .imp)

stripplot(cog_vars_v23_imp, wlogSCB_ConStroopRTms_v2_rev +
            wlogSCB_InconStroopRTms_v2_rev  ~ .imp)


stripplot(cog_vars_v23_imp, wCDR_ImmWRCorrect_v2 +
            wCDR_DelayedWRCorrect_v2 +
            wSerial3sAcc_v2 +
            wSerial7sAcc_v2 +
            wCDR_NumericWMOverallAcc_v2  ~ .imp)

# data for CFA 

varTable(imp_dat_v2)

#### six months CFA data ####

v3 <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\BLINDED DATA SETS\\BLINDED _DO NOT TOUCH\\ARCLI_V3_6 Month_ALL_6-SEPT-2021.sav")

# remove SPSS missing data indicators

v3 <- v3 %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., 999)))

# remove blank rows

v3 <- v3 %>% 
  filter(!IDno == "")

# Rename problematic variable names

v3 <- v3 %>% 
  rename(CDR_ImmWRCorrect_v3 = `CDR_ImmWRCorrect#_v3`)

v3 <- v3 %>% 
  rename(CDR_DelayedWRCorrect_v3 = `CDR_DelayedWRCorrect#_v3`)

# create new overall score

v3$CDR_NumericWMOverallAcc_v3 <- (v3$CDR_NumericWMOrigStimAcc_v3 + v3$CDR_NumericWMNewStimAcc_v3) / 2

# set sex variable to factor

v3$Sex <- as.factor(v3$Sex)

levels(v3$Sex) <- c("male", "female")

# Add age variable to v3 dataset 

age_bl <- arcli %>% select(IDno, Initialage)

v3 <- v3 %>% 
  left_join(age_bl, by = "IDno")

# select relevant variables 

cog_vars_v3 <- v3 %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v3, Initialage,CDR_RVIPacc_v3,CDR_ImmWRCorrect_v3, CDR_DelayedWRCorrect_v3, 
                CDR_SimpleRTms_v3, CDR_ChoiceRTms_v3, CDR_NumericWMOverallAcc_v3, 
                CDR_DigitVigRTms_v3, SCB_ConStroopRTms_v3, SCB_InconStroopRTms_v3, 
                Jenson8ChoiceDecisionTimems_v3,Jenson4ChoiceDecisionTimems_v3,
                Jenson1ChoiceDecisionTimems_v3,Serial3sAcc_v3, Serial7sAcc_v3)

## Remove rows without cog results

cog_vars_v3 <- cog_vars_v3[which(rowMeans(is.na(cog_vars_v3)) < 0.80),]

# Remove error from congruent stroop RT variable

cog_vars_v3 <- cog_vars_v3 %>% 
  mutate(SCB_ConStroopRTms_v3 = na_if(SCB_ConStroopRTms_v3, 0))

# log transform RT variables

RT_vars <- cog_vars_v3 %>% select(
  CDR_RVIPRTms_v3,CDR_SimpleRTms_v3, CDR_ChoiceRTms_v3, 
  CDR_DigitVigRTms_v3, SCB_ConStroopRTms_v3, SCB_InconStroopRTms_v3, 
  Jenson8ChoiceDecisionTimems_v3,Jenson4ChoiceDecisionTimems_v3,
  Jenson1ChoiceDecisionTimems_v3) %>% 
  names() # variable names

cog_vars_v3 <- cog_vars_v3 %>% 
  mutate(across(all_of(RT_vars), ~ log2(.), .names = "log{.col}")) # log 2 transform

# check histograms #

cog_vars_v3 %>% 
  select(starts_with("log")) %>% 
  multi.hist(global = F)

cog_vars_v3 %>% select(CDR_ImmWRCorrect_v3, CDR_DelayedWRCorrect_v3,
                       CDR_NumericWMOverallAcc_v3,
                       Serial3sAcc_v3, Serial7sAcc_v3) %>% 
  multi.hist(global = F)

# Winsorise outliers 

cog_vars_v3 <- cog_vars_v3 %>% # winsorise everything except age
  mutate(across(c(where(is.numeric), -Initialage), winsorise, sd = 3, .names = "w{.col}"))

cog_vars_v32 <- cog_vars_v3 %>% select(IDno, Initialage, Sex, starts_with("wlog"), 
                                       wCDR_ImmWRCorrect_v3,wCDR_DelayedWRCorrect_v3, 
                                       wCDR_NumericWMOverallAcc_v3, wSerial3sAcc_v3,
                                       wSerial7sAcc_v3) # select relevant variables

# histograms of winsorised variables

cog_vars_v32 %>% select(starts_with("wlog")) %>% 
  multi.hist(global = F)

cog_vars_v32 %>% select(wCDR_ImmWRCorrect_v3, wCDR_DelayedWRCorrect_v3,
                        wCDR_NumericWMOverallAcc_v3,wSerial3sAcc_v3, wSerial7sAcc_v3) %>% 
  multi.hist(global = F)

## remove cases with more than half cog vars missing

nrow(cog_vars_v32)

cog_vars_v32 <- cog_vars_v32[which(rowMeans(is.na(cog_vars_v32)) < 0.50),]

nrow(cog_vars_v32)

# check missing data per variable

n_miss <- colSums(is.na(cog_vars_v32))

median(n_miss / nrow(cog_vars_v32))

sort(n_miss / nrow(cog_vars_v32))

# Reverse score reaction time measures

cog_vars_v32 <- cog_vars_v32 %>% 
  mutate(across(starts_with("wlog"), reverse_RT, .names = "{.col}_rev"))

## rescale accuracy variables (10% units)

cog_vars_v32 <- cog_vars_v32 %>% mutate(across(starts_with("wSerial"), ~  . / 10))

cog_vars_v32 <- cog_vars_v32 %>% mutate(wCDR_NumericWMOverallAcc_v3 = wCDR_NumericWMOverallAcc_v3 / 10)

# select relevant variables #

cog_vars_v33 <- cog_vars_v32 %>% 
  select(IDno, Sex, Initialage, ends_with("_rev"), wCDR_ImmWRCorrect_v3,
         wCDR_DelayedWRCorrect_v3, wCDR_NumericWMOverallAcc_v3,
         wSerial3sAcc_v3, wSerial7sAcc_v3)

## single imputation

# predictive mean matching

predictormatrix <- quickpred(cog_vars_v33, 
                             include=c("Sex", "Initialage"),
                             exclude="IDno",
                             mincor = 0.1)

cog_vars_v33_imp <- mice(cog_vars_v33, m = 1, method = "pmm",
                         predictorMatrix = predictormatrix, diagnostics = T,
                         seed = 10)

imp_dat_v3 <- complete(cog_vars_v33_imp) # imputed data 

# check imputation for each variable 

stripplot(cog_vars_v33_imp, wlogCDR_SimpleRTms_v3_rev +
            wlogCDR_ChoiceRTms_v3_rev +
            wlogCDR_DigitVigRTms_v3_rev +
            wlogJenson8ChoiceDecisionTimems_v3_rev +
            wlogJenson4ChoiceDecisionTimems_v3_rev +
            wlogJenson1ChoiceDecisionTimems_v3_rev ~ .imp)

stripplot(cog_vars_v33_imp, wlogSCB_ConStroopRTms_v3_rev +
            wlogSCB_InconStroopRTms_v3_rev  ~ .imp)


stripplot(cog_vars_v33_imp, wCDR_ImmWRCorrect_v3 +
            wCDR_DelayedWRCorrect_v3 +
            wSerial3sAcc_v3 +
            wSerial7sAcc_v3 +
            wCDR_NumericWMOverallAcc_v3  ~ .imp)

# data for CFA 

varTable(imp_dat_v3)

#### 12 months CFA data ####

v4 <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\BLINDED DATA SETS\\BLINDED _DO NOT TOUCH\\ARCLI_V4_12 Month_ALL_6-SEPT-2021.sav")

# remove SPSS missing data indicators

v4 <- v4 %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -99))) %>% 
  mutate(across(where(is.numeric), ~ na_if(., 999)))

# remove blank rows

v4 <- v4 %>% 
  filter(!IDno == "")

# Rename problematic variable names

v4 <- v4 %>% 
  rename(CDR_ImmWRCorrect_v4 = `CDR_ImmWRCorrect#_v4`)

v4 <- v4 %>% 
  rename(CDR_DelayedWRCorrect_v4 = `CDR_DelayedWRCorrect#_v4`)

v4 <- v4 %>% 
  rename(Serial3sAcc_v4 = Serial3sAcc_v41)

# create new overall score

v4$CDR_NumericWMOverallAcc_v4 <- (v4$CDR_NumericWMOrigStimAcc_v4 + v4$CDR_NumericWMNewStimAcc_v4) / 2

# set sex variable to factor

v4$Sex <- as.factor(v4$Sex)

levels(v4$Sex) <- c("male", "female")

# Add age variable to v4 dataset 

age_bl <- arcli %>% select(IDno, Initialage)

v4 <- v4 %>% 
  left_join(age_bl, by = "IDno")

# select relevant variables 

cog_vars_v4 <- v4 %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v4, Initialage,CDR_RVIPacc_v4,CDR_ImmWRCorrect_v4, CDR_DelayedWRCorrect_v4, 
                CDR_SimpleRTms_v4, CDR_ChoiceRTms_v4, CDR_NumericWMOverallAcc_v4, 
                CDR_DigitVigRTms_v4, SCB_ConStroopRTms_v4, SCB_InconStroopRTms_v4, 
                Jenson8ChoiceDecisionTimems_v4,Jenson4ChoiceDecisionTimems_v4,
                Jenson1ChoiceDecisionTimems_v4,Serial3sAcc_v4, Serial7sAcc_v4)

## Remove rows without cog results

cog_vars_v4 <- cog_vars_v4[which(rowMeans(is.na(cog_vars_v4)) < 0.80),]

# log transform RT variables

RT_vars <- cog_vars_v4 %>% select(
  CDR_RVIPRTms_v4,CDR_SimpleRTms_v4, CDR_ChoiceRTms_v4, 
  CDR_DigitVigRTms_v4, SCB_ConStroopRTms_v4, SCB_InconStroopRTms_v4, 
  Jenson8ChoiceDecisionTimems_v4,Jenson4ChoiceDecisionTimems_v4,
  Jenson1ChoiceDecisionTimems_v4) %>% 
  names() # variable names

cog_vars_v4 <- cog_vars_v4 %>% 
  mutate(across(all_of(RT_vars), ~ log2(.), .names = "log{.col}")) # log 2 transform

# check histograms #

cog_vars_v4 %>% 
  select(starts_with("log")) %>% 
  multi.hist(global = F)

cog_vars_v4 %>% select(CDR_ImmWRCorrect_v4, CDR_DelayedWRCorrect_v4,
                       CDR_NumericWMOverallAcc_v4,
                       Serial3sAcc_v4, Serial7sAcc_v4) %>% 
  multi.hist(global = F)

# Winsorise outliers 

cog_vars_v4 <- cog_vars_v4 %>% # winsorise everything except age
  mutate(across(c(where(is.numeric), -Initialage), winsorise, sd = 3, .names = "w{.col}"))

cog_vars_v42 <- cog_vars_v4 %>% select(IDno, Initialage, Sex, starts_with("wlog"), 
                                       wCDR_ImmWRCorrect_v4,wCDR_DelayedWRCorrect_v4, 
                                       wCDR_NumericWMOverallAcc_v4, wSerial3sAcc_v4,
                                       wSerial7sAcc_v4) # select relevant variables

# histograms of winsorised variables

cog_vars_v42 %>% select(starts_with("wlog")) %>% 
  multi.hist(global = F)

cog_vars_v42 %>% select(wCDR_ImmWRCorrect_v4, wCDR_DelayedWRCorrect_v4,
                        wCDR_NumericWMOverallAcc_v4,wSerial3sAcc_v4, wSerial7sAcc_v4) %>% 
  multi.hist(global = F)

## remove cases with more than half cog vars missing

nrow(cog_vars_v42)

cog_vars_v42 <- cog_vars_v42[which(rowMeans(is.na(cog_vars_v42)) < 0.50),]

nrow(cog_vars_v42)

# check missing data per variable

n_miss <- colSums(is.na(cog_vars_v42))

median(n_miss / nrow(cog_vars_v42))

sort(n_miss / nrow(cog_vars_v42))

# Reverse score reaction time measures

cog_vars_v42 <- cog_vars_v42 %>% 
  mutate(across(starts_with("wlog"), reverse_RT, .names = "{.col}_rev"))

## rescale accuracy variables (10% units)

cog_vars_v42 <- cog_vars_v42 %>% mutate(across(starts_with("wSerial"), ~  . / 10))

cog_vars_v42 <- cog_vars_v42 %>% mutate(wCDR_NumericWMOverallAcc_v4 = wCDR_NumericWMOverallAcc_v4 / 10)

# select relevant variables #

cog_vars_v43 <- cog_vars_v42 %>% 
  select(IDno, Sex, Initialage, ends_with("_rev"), wCDR_ImmWRCorrect_v4,
         wCDR_DelayedWRCorrect_v4, wCDR_NumericWMOverallAcc_v4,
         wSerial3sAcc_v4, wSerial7sAcc_v4)

## single imputation

# predictive mean matching

predictormatrix <- quickpred(cog_vars_v43, 
                             include=c("Sex", "Initialage"),
                             exclude="IDno",
                             mincor = 0.1)

cog_vars_v43_imp <- mice(cog_vars_v43, m = 1, method = "pmm",
                         predictorMatrix = predictormatrix, diagnostics = T,
                         seed = 10)

imp_dat_v4 <- complete(cog_vars_v43_imp) # imputed data 

# check imputation for each variable 

stripplot(cog_vars_v43_imp, wlogCDR_SimpleRTms_v4_rev +
            wlogCDR_ChoiceRTms_v4_rev +
            wlogCDR_DigitVigRTms_v4_rev +
            wlogJenson8ChoiceDecisionTimems_v4_rev +
            wlogJenson4ChoiceDecisionTimems_v4_rev +
            wlogJenson1ChoiceDecisionTimems_v4_rev ~ .imp)

stripplot(cog_vars_v43_imp, wlogSCB_ConStroopRTms_v4_rev +
            wlogSCB_InconStroopRTms_v4_rev  ~ .imp)


stripplot(cog_vars_v43_imp, wCDR_ImmWRCorrect_v4 +
            wCDR_DelayedWRCorrect_v4 +
            wSerial3sAcc_v4 +
            wSerial7sAcc_v4 +
            wCDR_NumericWMOverallAcc_v4  ~ .imp)

# data for CFA 

varTable(imp_dat_v4)


#### Combine data from each visit ####

cfa_data <- imp_dat %>% 
  full_join(imp_dat_v2, by = c("IDno", "Initialage" , "Sex")) %>% 
  full_join(imp_dat_v3, by = c("IDno", "Initialage" , "Sex")) %>% 
  full_join(imp_dat_v4, by = c("IDno", "Initialage" , "Sex"))


# subtract mean from each variable

cfa_data <- cfa_data %>% 
  mutate(across(where(is.numeric), ~ . - mean(., na.rm=T)))

#### CFA ####

## longitudinal covariances

d <- tibble(
  a = c("wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v2_rev",
  "wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v3_rev",
  "wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev",
  "wlogCDR_SimpleRTms_v2_rev ~~ g1*wlogCDR_SimpleRTms_v3_rev",
  "wlogCDR_SimpleRTms_v2_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev",
  "wlogCDR_SimpleRTms_v3_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev"))


long_cov <- d %>% 
  mutate(b = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogCDR_ChoiceRTms_")) %>% 
  mutate(b = str_replace_all(b, "g1", "g2")) %>% 
  mutate(c = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogCDR_DigitVigRTms_")) %>% 
  mutate(c = str_replace_all(c, "g1", "g3")) %>% 
  mutate(d = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogJenson8ChoiceDecisionTimems_")) %>% 
  mutate(d = str_replace_all(d, "g1", "g4")) %>% 
  mutate(e = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogJenson4ChoiceDecisionTimems_")) %>% 
  mutate(e = str_replace_all(e, "g1", "g5")) %>% 
  mutate(f = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogJenson1ChoiceDecisionTimems_")) %>% 
  mutate(f = str_replace_all(f, "g1", "g6")) %>% 
  mutate(g = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogSCB_ConStroopRTms_")) %>% 
  mutate(g = str_replace_all(g, "g1", "g7")) %>%
  mutate(h = str_replace_all(a, "wlogCDR_SimpleRTms_", "wlogSCB_InconStroopRTms_")) %>% 
  mutate(h = str_replace_all(h, "g1", "g8")) %>%
  mutate(i = str_replace_all(a, "wlogCDR_SimpleRTms_", "wCDR_ImmWRCorrect_")) %>% 
  mutate(i = str_replace_all(i, "g1", "g9")) %>%
  mutate(j = str_replace_all(a, "wlogCDR_SimpleRTms_", "wCDR_DelayedWRCorrect_")) %>% 
  mutate(j = str_replace_all(j, "g1", "g10")) %>%
  mutate(k = str_replace_all(a, "wlogCDR_SimpleRTms_", "wSerial3sAcc_")) %>% 
  mutate(k = str_replace_all(k, "g1", "g11")) %>%
  mutate(l = str_replace_all(a, "wlogCDR_SimpleRTms_", "wSerial7sAcc_")) %>% 
  mutate(l = str_replace_all(l, "g1", "g12")) %>%
  mutate(m = str_replace_all(a, "wlogCDR_SimpleRTms_", "wCDR_NumericWMOverallAcc_")) %>% 
  mutate(m = str_replace_all(m, "g1", "g13")) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = case_when(
    str_detect(value, "wCDR_ImmWRCorrect_") ~ str_replace_all(value, "_rev", ""),
    str_detect(value, "wCDR_DelayedWRCorrect_") ~ str_replace_all(value, "_rev", ""),
    str_detect(value, "wSerial3sAcc_") ~ str_replace_all(value, "_rev", ""),
    str_detect(value, "wSerial7sAcc_") ~ str_replace_all(value, "_rev", ""),
    str_detect(value, "wCDR_NumericWMOverallAcc_") ~ str_replace_all(value, "_rev", ""),
    TRUE ~ value
  ))

writeClipboard(long_cov$value)

## Latent variables covariances

row.names(cov(cfa_data[,-c(1:3)]))

colnames(cov(cfa_data[,-c(1:3)]))

covar <- row.names(lavaan::vcov(test)) %>% as_tibble()

covar <- covar[161:296,]

covar <- covar %>% 
  mutate(coef = case_when(
    str_detect(value, "Gt") & str_detect(value, "Gs") ~ "v1",
    str_detect(value, "Gt") & str_detect(value, "Glr") ~ "v2",
    str_detect(value, "Gt") & str_detect(value, "Gwm") ~ "v3",
    str_detect(value, "Gs") & str_detect(value, "Glr") ~ "v4",
    str_detect(value, "Gs") & str_detect(value, "Gwm") ~ "v5",
    str_detect(value, "Glr") & str_detect(value, "Gwm") ~ "v6",
    TRUE ~ ""
  )) %>% 
  separate(value, into = c("a", "b"), sep = "~~") %>% 
  mutate(covar1 = if_else(!coef == "", paste(coef, b, sep = "*"), b)) %>% 
  mutate(covar2 = paste(a, covar1, sep = " ~~ ")) %>% 
  filter(!str_detect(covar2, "~~ G"))

writeClipboard(covar$covar2)

## Shared method variances

d <- expand.grid(
  c("wlogJenson4ChoiceDecisionTimems_v1_rev", "wlogJenson4ChoiceDecisionTimems_v2_rev", "wlogJenson4ChoiceDecisionTimems_v3_rev", "wlogJenson4ChoiceDecisionTimems_v4_rev"),
  c("wlogJenson1ChoiceDecisionTimems_v1_rev", "wlogJenson1ChoiceDecisionTimems_v2_rev", "wlogJenson1ChoiceDecisionTimems_v3_rev", "wlogJenson1ChoiceDecisionTimems_v4_rev")
) %>% 
  as_tibble() %>% 
  mutate(a = paste(Var1, Var2, sep = " ~~ m1*")) %>% 
  select(-Var1, -Var2)

d

jens <- d %>% 
  mutate(b = str_replace(a, "wlogJenson4ChoiceDecisionTimems_", "wlogJenson8ChoiceDecisionTimems_")) %>% 
  mutate(b = str_replace(b, "m1", "m2")) %>% 
  mutate(c = str_replace(a, "wlogJenson1ChoiceDecisionTimems_", "wlogJenson8ChoiceDecisionTimems_")) %>% 
  mutate(c = str_replace(c, "m1", "m3")) %>% 
  pivot_longer(everything())


writeClipboard(jens$value)

# define CFA model 

model <- 
  
' # model 
  
## Factor loadings 
  
Gt_v1 =~ NA*wlogCDR_SimpleRTms_v1_rev +
       a*wlogCDR_SimpleRTms_v1_rev +
       a1*wlogCDR_ChoiceRTms_v1_rev +
       a2*wlogCDR_DigitVigRTms_v1_rev +
       a3*wlogJenson8ChoiceDecisionTimems_v1_rev +
       a4*wlogJenson4ChoiceDecisionTimems_v1_rev +
       a5*wlogJenson1ChoiceDecisionTimems_v1_rev
       
Gs_v1 =~  NA*wlogSCB_ConStroopRTms_v1_rev +
       b*wlogSCB_ConStroopRTms_v1_rev +
       b1*wlogSCB_InconStroopRTms_v1_rev +
       # cross loading
       b2*wlogCDR_ChoiceRTms_v1_rev 

       
Glr_m6_v1 =~ NA*wCDR_ImmWRCorrect_v1 +
       c*wCDR_ImmWRCorrect_v1 +
       c1*wCDR_DelayedWRCorrect_v1 
       

Gwm_ac_v1 =~ NA*wSerial3sAcc_v1 +
        d*wSerial3sAcc_v1 +
        d1*wSerial7sAcc_v1 +
        d2*wCDR_NumericWMOverallAcc_v1
        
Gt_v2 =~ NA*wlogCDR_SimpleRTms_v2_rev +
       a*wlogCDR_SimpleRTms_v2_rev +
       a1*wlogCDR_ChoiceRTms_v2_rev +
       a2*wlogCDR_DigitVigRTms_v2_rev +
       a3*wlogJenson8ChoiceDecisionTimems_v2_rev +
       a4*wlogJenson4ChoiceDecisionTimems_v2_rev +
       a5*wlogJenson1ChoiceDecisionTimems_v2_rev
       
Gs_v2 =~  NA*wlogSCB_ConStroopRTms_v2_rev +
       b*wlogSCB_ConStroopRTms_v2_rev +
       b1*wlogSCB_InconStroopRTms_v2_rev +
       # cross loading
       b2*wlogCDR_ChoiceRTms_v2_rev 

       
Glr_m6_v2 =~ NA*wCDR_ImmWRCorrect_v2 +
       c*wCDR_ImmWRCorrect_v2 +
       c1*wCDR_DelayedWRCorrect_v2 
       

Gwm_ac_v2 =~ NA*wSerial3sAcc_v2 +
        d*wSerial3sAcc_v2 +
        d1*wSerial7sAcc_v2 +
        d2*wCDR_NumericWMOverallAcc_v2
        
Gt_v3 =~ NA*wlogCDR_SimpleRTms_v3_rev +
       a*wlogCDR_SimpleRTms_v3_rev +
       a1*wlogCDR_ChoiceRTms_v3_rev +
       a2*wlogCDR_DigitVigRTms_v3_rev +
       a3*wlogJenson8ChoiceDecisionTimems_v3_rev +
       a4*wlogJenson4ChoiceDecisionTimems_v3_rev +
       a5*wlogJenson1ChoiceDecisionTimems_v3_rev
       
Gs_v3 =~ NA*wlogSCB_ConStroopRTms_v3_rev +
       b*wlogSCB_ConStroopRTms_v3_rev +
       b1*wlogSCB_InconStroopRTms_v3_rev +
       # cross loading
       b2*wlogCDR_ChoiceRTms_v3_rev 

       
Glr_m6_v3 =~ NA*wCDR_ImmWRCorrect_v3 +
       c*wCDR_ImmWRCorrect_v3 +
       c1*wCDR_DelayedWRCorrect_v3 
       

Gwm_ac_v3 =~ NA*wSerial3sAcc_v3 +
        d*wSerial3sAcc_v3 +
        d1*wSerial7sAcc_v3 +
        d2*wCDR_NumericWMOverallAcc_v3
        

Gt_v4 =~ NA*wlogCDR_SimpleRTms_v4_rev +
       a*wlogCDR_SimpleRTms_v4_rev +
       a1*wlogCDR_ChoiceRTms_v4_rev +
       a2*wlogCDR_DigitVigRTms_v4_rev +
       a3*wlogJenson8ChoiceDecisionTimems_v4_rev +
       a4*wlogJenson4ChoiceDecisionTimems_v4_rev +
       a5*wlogJenson1ChoiceDecisionTimems_v4_rev
       
Gs_v4 =~ NA*wlogSCB_ConStroopRTms_v4_rev +
       b*wlogSCB_ConStroopRTms_v4_rev +
       b1*wlogSCB_InconStroopRTms_v4_rev +
       # cross loading
       b2*wlogCDR_ChoiceRTms_v4_rev 

       
Glr_m6_v4 =~ NA*wCDR_ImmWRCorrect_v4 +
       c*wCDR_ImmWRCorrect_v4 +
       c1*wCDR_DelayedWRCorrect_v4 
       

Gwm_ac_v4 =~ NA*wSerial3sAcc_v4 +
        d*wSerial3sAcc_v4 +
        d1*wSerial7sAcc_v4 +
        d2*wCDR_NumericWMOverallAcc_v4       
        
## Manifest intercepts

# v1

wlogCDR_SimpleRTms_v1_rev ~ i1*1
wlogCDR_ChoiceRTms_v1_rev ~ i2*1
wlogCDR_DigitVigRTms_v1_rev ~ i3*1
wlogJenson8ChoiceDecisionTimems_v1_rev ~ i4*1
wlogJenson4ChoiceDecisionTimems_v1_rev ~ i5*1
wlogJenson1ChoiceDecisionTimems_v1_rev ~ i6*1
wlogSCB_ConStroopRTms_v1_rev ~ i7*1
wlogSCB_InconStroopRTms_v1_rev ~ i8*1
wCDR_ImmWRCorrect_v1 ~ i9*1
wCDR_DelayedWRCorrect_v1 ~ i10*1
wSerial3sAcc_v1 ~ i11*1
wSerial7sAcc_v1 ~ i12*1
wCDR_NumericWMOverallAcc_v1 ~ i13*1

#v2

wlogCDR_SimpleRTms_v2_rev ~ i1*1
wlogCDR_ChoiceRTms_v2_rev ~ i2*1
wlogCDR_DigitVigRTms_v2_rev ~ i3*1
wlogJenson8ChoiceDecisionTimems_v2_rev ~ i4*1
wlogJenson4ChoiceDecisionTimems_v2_rev ~ i5*1
wlogJenson1ChoiceDecisionTimems_v2_rev ~ i6*1
wlogSCB_ConStroopRTms_v2_rev ~ i7*1
wlogSCB_InconStroopRTms_v2_rev ~ i8*1
wCDR_ImmWRCorrect_v2 ~ i9*1
wCDR_DelayedWRCorrect_v2 ~ i10*1
wSerial3sAcc_v2 ~ i11*1
wSerial7sAcc_v2 ~ i12*1
wCDR_NumericWMOverallAcc_v2 ~ i13*1

#v3

wlogCDR_SimpleRTms_v3_rev ~ i1*1
wlogCDR_ChoiceRTms_v3_rev ~ i2*1
wlogCDR_DigitVigRTms_v3_rev ~ i3*1
wlogJenson8ChoiceDecisionTimems_v3_rev ~ i4*1
wlogJenson4ChoiceDecisionTimems_v3_rev ~ i5*1
wlogJenson1ChoiceDecisionTimems_v3_rev ~ i6*1
wlogSCB_ConStroopRTms_v3_rev ~ i7*1
wlogSCB_InconStroopRTms_v3_rev ~ i8*1
wCDR_ImmWRCorrect_v3 ~ i9*1
wCDR_DelayedWRCorrect_v3 ~ i10*1
wSerial3sAcc_v3 ~ i11*1
wSerial7sAcc_v3 ~ i12*1
wCDR_NumericWMOverallAcc_v3 ~ i13*1

#v4

wlogCDR_SimpleRTms_v4_rev ~ i1*1
wlogCDR_ChoiceRTms_v4_rev ~ i2*1
wlogCDR_DigitVigRTms_v4_rev ~ i3*1
wlogJenson8ChoiceDecisionTimems_v4_rev ~ i4*1
wlogJenson4ChoiceDecisionTimems_v4_rev ~ i5*1
wlogJenson1ChoiceDecisionTimems_v4_rev ~ i6*1
wlogSCB_ConStroopRTms_v4_rev ~ i7*1
wlogSCB_InconStroopRTms_v4_rev ~ i8*1
wCDR_ImmWRCorrect_v4 ~ i9*1
wCDR_DelayedWRCorrect_v4 ~ i10*1
wSerial3sAcc_v4 ~ i11*1
wSerial7sAcc_v4 ~ i12*1
wCDR_NumericWMOverallAcc_v4 ~ i13*1


## latent means

Gt_v1 ~ 0*1
Gt_v2 ~ 1
Gt_v3 ~ 1
Gt_v4 ~ 1

Gs_v1 ~ 0*1
Gs_v2 ~ 1
Gs_v3 ~ 1
Gs_v4 ~ 1

Gwm_ac_v1 ~ 0*1
Gwm_ac_v2 ~ 1
Gwm_ac_v3 ~ 1
Gwm_ac_v4 ~ 1

Glr_m6_v1 ~ 0*1
Glr_m6_v2 ~ 1
Glr_m6_v3 ~ 1
Glr_m6_v4 ~ 1


## Observed variances

# v1

wlogCDR_SimpleRTms_v1_rev ~~ u1*wlogCDR_SimpleRTms_v1_rev
wlogCDR_ChoiceRTms_v1_rev ~~ u2*wlogCDR_ChoiceRTms_v1_rev
wlogCDR_DigitVigRTms_v1_rev ~~ u3*wlogCDR_DigitVigRTms_v1_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ u4*wlogJenson8ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ u5*wlogJenson4ChoiceDecisionTimems_v1_rev
wlogJenson1ChoiceDecisionTimems_v1_rev ~~ u6*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogSCB_ConStroopRTms_v1_rev ~~ u7*wlogSCB_ConStroopRTms_v1_rev
wlogSCB_InconStroopRTms_v1_rev ~~ u8*wlogSCB_InconStroopRTms_v1_rev
wCDR_ImmWRCorrect_v1 ~~ u9*wCDR_ImmWRCorrect_v1
wCDR_DelayedWRCorrect_v1 ~~ u10*wCDR_DelayedWRCorrect_v1
wSerial3sAcc_v1 ~~ u11*wSerial3sAcc_v1
wSerial7sAcc_v1 ~~ u12*wSerial7sAcc_v1
wCDR_NumericWMOverallAcc_v1 ~~ u13*wCDR_NumericWMOverallAcc_v1

# v2

wlogCDR_SimpleRTms_v2_rev ~~ u1*wlogCDR_SimpleRTms_v2_rev
wlogCDR_ChoiceRTms_v2_rev ~~ u2*wlogCDR_ChoiceRTms_v2_rev
wlogCDR_DigitVigRTms_v2_rev ~~ u3*wlogCDR_DigitVigRTms_v2_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ u4*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ u5*wlogJenson4ChoiceDecisionTimems_v2_rev
wlogJenson1ChoiceDecisionTimems_v2_rev ~~ u6*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogSCB_ConStroopRTms_v2_rev ~~ u7*wlogSCB_ConStroopRTms_v2_rev
wlogSCB_InconStroopRTms_v2_rev ~~ u8*wlogSCB_InconStroopRTms_v2_rev
wCDR_ImmWRCorrect_v2 ~~ u9*wCDR_ImmWRCorrect_v2
wCDR_DelayedWRCorrect_v2 ~~ u10*wCDR_DelayedWRCorrect_v2
wSerial3sAcc_v2 ~~ u11*wSerial3sAcc_v2
wSerial7sAcc_v2 ~~ u12*wSerial7sAcc_v2
wCDR_NumericWMOverallAcc_v2 ~~ u13*wCDR_NumericWMOverallAcc_v2

# v3

wlogCDR_SimpleRTms_v3_rev ~~ u1*wlogCDR_SimpleRTms_v3_rev
wlogCDR_ChoiceRTms_v3_rev ~~ u2*wlogCDR_ChoiceRTms_v3_rev
wlogCDR_DigitVigRTms_v3_rev ~~ u3*wlogCDR_DigitVigRTms_v3_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ u4*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ u5*wlogJenson4ChoiceDecisionTimems_v3_rev
wlogJenson1ChoiceDecisionTimems_v3_rev ~~ u6*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogSCB_ConStroopRTms_v3_rev ~~ u7*wlogSCB_ConStroopRTms_v3_rev
wlogSCB_InconStroopRTms_v3_rev ~~ u8*wlogSCB_InconStroopRTms_v3_rev
wCDR_ImmWRCorrect_v3 ~~ u9*wCDR_ImmWRCorrect_v3
wCDR_DelayedWRCorrect_v3 ~~ u10*wCDR_DelayedWRCorrect_v3
wSerial3sAcc_v3 ~~ u11*wSerial3sAcc_v3
wSerial7sAcc_v3 ~~ u12*wSerial7sAcc_v3
wCDR_NumericWMOverallAcc_v3 ~~ u13*wCDR_NumericWMOverallAcc_v3

# v4

wlogCDR_SimpleRTms_v4_rev ~~ u1*wlogCDR_SimpleRTms_v4_rev
wlogCDR_ChoiceRTms_v4_rev ~~ u2*wlogCDR_ChoiceRTms_v4_rev
wlogCDR_DigitVigRTms_v4_rev ~~ u3*wlogCDR_DigitVigRTms_v4_rev
wlogJenson8ChoiceDecisionTimems_v4_rev ~~ u4*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ u5*wlogJenson4ChoiceDecisionTimems_v4_rev
wlogJenson1ChoiceDecisionTimems_v4_rev ~~ u6*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogSCB_ConStroopRTms_v4_rev ~~ u7*wlogSCB_ConStroopRTms_v4_rev
wlogSCB_InconStroopRTms_v4_rev ~~ u8*wlogSCB_InconStroopRTms_v4_rev
wCDR_ImmWRCorrect_v4 ~~ u9*wCDR_ImmWRCorrect_v4
wCDR_DelayedWRCorrect_v4 ~~ u10*wCDR_DelayedWRCorrect_v4
wSerial3sAcc_v4 ~~ u11*wSerial3sAcc_v4
wSerial7sAcc_v4 ~~ u12*wSerial7sAcc_v4
wCDR_NumericWMOverallAcc_v4 ~~ u13*wCDR_NumericWMOverallAcc_v4

## longitudinal covariances

wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v2_rev
wlogCDR_ChoiceRTms_v1_rev ~~ g2*wlogCDR_ChoiceRTms_v2_rev
wlogCDR_DigitVigRTms_v1_rev ~~ g3*wlogCDR_DigitVigRTms_v2_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v2_rev
wlogJenson1ChoiceDecisionTimems_v1_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogSCB_ConStroopRTms_v1_rev ~~ g7*wlogSCB_ConStroopRTms_v2_rev
wlogSCB_InconStroopRTms_v1_rev ~~ g8*wlogSCB_InconStroopRTms_v2_rev
wCDR_ImmWRCorrect_v1 ~~ g9*wCDR_ImmWRCorrect_v2
wCDR_DelayedWRCorrect_v1 ~~ g10*wCDR_DelayedWRCorrect_v2
wSerial3sAcc_v1 ~~ g11*wSerial3sAcc_v2
wSerial7sAcc_v1 ~~ g12*wSerial7sAcc_v2
wCDR_NumericWMOverallAcc_v1 ~~ g13*wCDR_NumericWMOverallAcc_v2
wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v3_rev
wlogCDR_ChoiceRTms_v1_rev ~~ g2*wlogCDR_ChoiceRTms_v3_rev
wlogCDR_DigitVigRTms_v1_rev ~~ g3*wlogCDR_DigitVigRTms_v3_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v3_rev
wlogJenson1ChoiceDecisionTimems_v1_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogSCB_ConStroopRTms_v1_rev ~~ g7*wlogSCB_ConStroopRTms_v3_rev
wlogSCB_InconStroopRTms_v1_rev ~~ g8*wlogSCB_InconStroopRTms_v3_rev
wCDR_ImmWRCorrect_v1 ~~ g9*wCDR_ImmWRCorrect_v3
wCDR_DelayedWRCorrect_v1 ~~ g10*wCDR_DelayedWRCorrect_v3
wSerial3sAcc_v1 ~~ g11*wSerial3sAcc_v3
wSerial7sAcc_v1 ~~ g12*wSerial7sAcc_v3
wCDR_NumericWMOverallAcc_v1 ~~ g13*wCDR_NumericWMOverallAcc_v3
wlogCDR_SimpleRTms_v1_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev
wlogCDR_ChoiceRTms_v1_rev ~~ g2*wlogCDR_ChoiceRTms_v4_rev
wlogCDR_DigitVigRTms_v1_rev ~~ g3*wlogCDR_DigitVigRTms_v4_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v4_rev
wlogJenson1ChoiceDecisionTimems_v1_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogSCB_ConStroopRTms_v1_rev ~~ g7*wlogSCB_ConStroopRTms_v4_rev
wlogSCB_InconStroopRTms_v1_rev ~~ g8*wlogSCB_InconStroopRTms_v4_rev
wCDR_ImmWRCorrect_v1 ~~ g9*wCDR_ImmWRCorrect_v4
wCDR_DelayedWRCorrect_v1 ~~ g10*wCDR_DelayedWRCorrect_v4
wSerial3sAcc_v1 ~~ g11*wSerial3sAcc_v4
wSerial7sAcc_v1 ~~ g12*wSerial7sAcc_v4
wCDR_NumericWMOverallAcc_v1 ~~ g13*wCDR_NumericWMOverallAcc_v4
wlogCDR_SimpleRTms_v2_rev ~~ g1*wlogCDR_SimpleRTms_v3_rev
wlogCDR_ChoiceRTms_v2_rev ~~ g2*wlogCDR_ChoiceRTms_v3_rev
wlogCDR_DigitVigRTms_v2_rev ~~ g3*wlogCDR_DigitVigRTms_v3_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v3_rev
wlogJenson1ChoiceDecisionTimems_v2_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogSCB_ConStroopRTms_v2_rev ~~ g7*wlogSCB_ConStroopRTms_v3_rev
wlogSCB_InconStroopRTms_v2_rev ~~ g8*wlogSCB_InconStroopRTms_v3_rev
wCDR_ImmWRCorrect_v2 ~~ g9*wCDR_ImmWRCorrect_v3
wCDR_DelayedWRCorrect_v2 ~~ g10*wCDR_DelayedWRCorrect_v3
wSerial3sAcc_v2 ~~ g11*wSerial3sAcc_v3
wSerial7sAcc_v2 ~~ g12*wSerial7sAcc_v3
wCDR_NumericWMOverallAcc_v2 ~~ g13*wCDR_NumericWMOverallAcc_v3
wlogCDR_SimpleRTms_v2_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev
wlogCDR_ChoiceRTms_v2_rev ~~ g2*wlogCDR_ChoiceRTms_v4_rev
wlogCDR_DigitVigRTms_v2_rev ~~ g3*wlogCDR_DigitVigRTms_v4_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v4_rev
wlogJenson1ChoiceDecisionTimems_v2_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogSCB_ConStroopRTms_v2_rev ~~ g7*wlogSCB_ConStroopRTms_v4_rev
wlogSCB_InconStroopRTms_v2_rev ~~ g8*wlogSCB_InconStroopRTms_v4_rev
wCDR_ImmWRCorrect_v2 ~~ g9*wCDR_ImmWRCorrect_v4
wCDR_DelayedWRCorrect_v2 ~~ g10*wCDR_DelayedWRCorrect_v4
wSerial3sAcc_v2 ~~ g11*wSerial3sAcc_v4
wSerial7sAcc_v2 ~~ g12*wSerial7sAcc_v4
wCDR_NumericWMOverallAcc_v2 ~~ g13*wCDR_NumericWMOverallAcc_v4
wlogCDR_SimpleRTms_v3_rev ~~ g1*wlogCDR_SimpleRTms_v4_rev
wlogCDR_ChoiceRTms_v3_rev ~~ g2*wlogCDR_ChoiceRTms_v4_rev
wlogCDR_DigitVigRTms_v3_rev ~~ g3*wlogCDR_DigitVigRTms_v4_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ g4*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ g5*wlogJenson4ChoiceDecisionTimems_v4_rev
wlogJenson1ChoiceDecisionTimems_v3_rev ~~ g6*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogSCB_ConStroopRTms_v3_rev ~~ g7*wlogSCB_ConStroopRTms_v4_rev
wlogSCB_InconStroopRTms_v3_rev ~~ g8*wlogSCB_InconStroopRTms_v4_rev
wCDR_ImmWRCorrect_v3 ~~ g9*wCDR_ImmWRCorrect_v4
wCDR_DelayedWRCorrect_v3 ~~ g10*wCDR_DelayedWRCorrect_v4
wSerial3sAcc_v3 ~~ g11*wSerial3sAcc_v4
wSerial7sAcc_v3 ~~ g12*wSerial7sAcc_v4
wCDR_NumericWMOverallAcc_v3 ~~ g13*wCDR_NumericWMOverallAcc_v4

## latent covariances 

Gt_v1 ~~ Gs_v1
Gt_v1 ~~ Glr_m6_v1
Gt_v1 ~~ Gwm_ac_v1
Gt_v1 ~~ Gs_v2
Gt_v1 ~~ Glr_m6_v2
Gt_v1 ~~ Gwm_ac_v2
Gt_v1 ~~ Gs_v3
Gt_v1 ~~ Glr_m6_v3
Gt_v1 ~~ Gwm_ac_v3
Gt_v1 ~~ Gs_v4
Gt_v1 ~~ Glr_m6_v4
Gt_v1 ~~ Gwm_ac_v4
Gs_v1 ~~ Glr_m6_v1
Gs_v1 ~~ Gwm_ac_v1
Gs_v1 ~~ Gt_v2
Gs_v1 ~~ Glr_m6_v2
Gs_v1 ~~ Gwm_ac_v2
Gs_v1 ~~ Gt_v3
Gs_v1 ~~ Glr_m6_v3
Gs_v1 ~~ Gwm_ac_v3
Gs_v1 ~~ Gt_v4
Gs_v1 ~~ Glr_m6_v4
Gs_v1 ~~ Gwm_ac_v4
Glr_m6_v1 ~~ Gwm_ac_v1
Glr_m6_v1 ~~ Gt_v2
Glr_m6_v1 ~~ Gs_v2
Glr_m6_v1 ~~ Gwm_ac_v2
Glr_m6_v1 ~~ Gt_v3
Glr_m6_v1 ~~ Gs_v3
Glr_m6_v1 ~~ Gwm_ac_v3
Glr_m6_v1 ~~ Gt_v4
Glr_m6_v1 ~~ Gs_v4
Glr_m6_v1 ~~ Gwm_ac_v4
Gwm_ac_v1 ~~ Gt_v2
Gwm_ac_v1 ~~ Gs_v2
Gwm_ac_v1 ~~ Glr_m6_v2
Gwm_ac_v1 ~~ Gt_v3
Gwm_ac_v1 ~~ Gs_v3
Gwm_ac_v1 ~~ Glr_m6_v3
Gwm_ac_v1 ~~ Gt_v4
Gwm_ac_v1 ~~ Gs_v4
Gwm_ac_v1 ~~ Glr_m6_v4
Gt_v2 ~~ Gs_v2
Gt_v2 ~~ Glr_m6_v2
Gt_v2 ~~ Gwm_ac_v2
Gt_v2 ~~ Gs_v3
Gt_v2 ~~ Glr_m6_v3
Gt_v2 ~~ Gwm_ac_v3
Gt_v2 ~~ Gs_v4
Gt_v2 ~~ Glr_m6_v4
Gt_v2 ~~ Gwm_ac_v4
Gs_v2 ~~ Glr_m6_v2
Gs_v2 ~~ Gwm_ac_v2
Gs_v2 ~~ Gt_v3
Gs_v2 ~~ Glr_m6_v3
Gs_v2 ~~ Gwm_ac_v3
Gs_v2 ~~ Gt_v4
Gs_v2 ~~ Glr_m6_v4
Gs_v2 ~~ Gwm_ac_v4
Glr_m6_v2 ~~ Gwm_ac_v2
Glr_m6_v2 ~~ Gt_v3
Glr_m6_v2 ~~ Gs_v3
Glr_m6_v2 ~~ Gwm_ac_v3
Glr_m6_v2 ~~ Gt_v4
Glr_m6_v2 ~~ Gs_v4
Glr_m6_v2 ~~ Gwm_ac_v4
Gwm_ac_v2 ~~ Gt_v3
Gwm_ac_v2 ~~ Gs_v3
Gwm_ac_v2 ~~ Glr_m6_v3
Gwm_ac_v2 ~~ Gt_v4
Gwm_ac_v2 ~~ Gs_v4
Gwm_ac_v2 ~~ Glr_m6_v4
Gt_v3 ~~ Gs_v3
Gt_v3 ~~ Glr_m6_v3
Gt_v3 ~~ Gwm_ac_v3
Gt_v3 ~~ Gs_v4
Gt_v3 ~~ Glr_m6_v4
Gt_v3 ~~ Gwm_ac_v4
Gs_v3 ~~ Glr_m6_v3
Gs_v3 ~~ Gwm_ac_v3
Gs_v3 ~~ Gt_v4
Gs_v3 ~~ Glr_m6_v4
Gs_v3 ~~ Gwm_ac_v4
Glr_m6_v3 ~~ Gwm_ac_v3
Glr_m6_v3 ~~ Gt_v4
Glr_m6_v3 ~~ Gs_v4
Glr_m6_v3 ~~ Gwm_ac_v4
Gwm_ac_v3 ~~ Gt_v4
Gwm_ac_v3 ~~ Gs_v4
Gwm_ac_v3 ~~ Glr_m6_v4
Gt_v4 ~~ Gs_v4
Gt_v4 ~~ Glr_m6_v4
Gt_v4 ~~ Gwm_ac_v4
Gs_v4 ~~ Glr_m6_v4
Gs_v4 ~~ Gwm_ac_v4
Glr_m6_v4 ~~ Gwm_ac_v4

## latent variances

Gt_v1 ~~ lv1*Gt_v1
Gt_v2 ~~ lv1*Gt_v2
Gt_v3 ~~ lv1*Gt_v3
Gt_v4 ~~ lv1*Gt_v4

Gs_v1 ~~ lv2*Gs_v1
Gs_v2 ~~ lv2*Gs_v2
Gs_v3 ~~ lv2*Gs_v3
Gs_v4 ~~ lv2*Gs_v4

Gwm_ac_v1 ~~ lv3*Gwm_ac_v1
Gwm_ac_v2 ~~ lv3*Gwm_ac_v2
Gwm_ac_v3 ~~ lv3*Gwm_ac_v3
Gwm_ac_v4 ~~ lv3*Gwm_ac_v4

Glr_m6_v1 ~~ lv4*Glr_m6_v1
Glr_m6_v2 ~~ lv4*Glr_m6_v2
Glr_m6_v3 ~~ lv4*Glr_m6_v3
Glr_m6_v4 ~~ lv4*Glr_m6_v4

## correlated residuals

wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v4_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v1_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson8ChoiceDecisionTimems_v4_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v2_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson8ChoiceDecisionTimems_v4_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v3_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson8ChoiceDecisionTimems_v2_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v2_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson8ChoiceDecisionTimems_v3_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v3_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m1*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson8ChoiceDecisionTimems_v4_rev ~~ m2*wlogJenson1ChoiceDecisionTimems_v4_rev
wlogJenson4ChoiceDecisionTimems_v4_rev ~~ m3*wlogJenson8ChoiceDecisionTimems_v4_rev
'

# fit model

test <- cfa(model, data = cfa_data, estimator = "MLR", std.lv = F)

# results

summary(test, fit.measures=T, standardized = T)

## Save factor scores ##

cfa_data <- cfa_data %>% 
  drop_na()

preds <- as_tibble(lavPredict(test, newdata = cfa_data, type = "lv",
                              method = "Bartlett"))


chc_factors <- cbind(cfa_data, preds)


#### Add into unblinded masterfile ####

master <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\ARCLI GROUPS\\UNBLINDED DATA SETS\\ARCLI_MASTERFILE_23-SEPT_2021.sav")

chc_factors <- chc_factors %>% 
  select(IDno, Gt_v1:Gwm_ac_v4)

master <- master %>% 
  full_join(chc_factors, by = "IDno")

write_sav(master, "master with chc.sav")
