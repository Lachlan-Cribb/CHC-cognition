library(haven)
library(tidyverse)
library(lavaan)
library(psych)
library(semTools)
library(mice)


## Load data 

arcli <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\BLINDED DATA SETS\\BLINDED _DO NOT TOUCH\\ARCLI_V1_Baseline_ALL_13-JUL-2021.sav")

# remove -999

arcli <- na_if(arcli, -999)

arcli <- na_if(arcli, -99)


# Rename problematic variable names


arcli <- arcli %>% 
  rename(CDR_ImmWRCorrect_v1 = `CDR_ImmWRCorrect#_v1`)

arcli <- arcli %>% 
  rename(CDR_DelayedWRCorrect_v1 = `CDR_DelayedWRCorrect#_v1`)


# create new overall score

arcli$CDR_NumericWMOverallAcc_v1 <- (arcli$CDR_NumericWMOrigStimAcc_v1 + arcli$CDR_NumericWMNewStimAcc_v1) / 2

# set sex variable to factor

arcli$Sex <- as.factor(arcli$Sex)

levels(arcli$Sex) <- c("male", "female")


# select relevant variables 

cog_vars <- arcli %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v1, Initialage,CDR_RVIPacc_v1,TMT_Bms_v1, CDR_ImmWRCorrect_v1, CDR_DelayedWRCorrect_v1, 
                CDR_SimpleRTms_v1, CDR_ChoiceRTms_v1, CDR_NumericWMOverallAcc_v1, 
                CDR_DigitVigRTms_v1, SCB_ConStroopRTms_v1, SCB_InconStroopRTms_v1, 
                Jenson8ChoiceDecisionTimems_v1,Jenson4ChoiceDecisionTimems_v1,
                Jenson1ChoiceDecisionTimems_v1,WASImatrixreasoningrawscore, 
                WASI_vocabrawscore, MMSE, Serial3sAcc_v1, Serial7sAcc_v1, TMT_Ams_v1)


## Remove rows without cog results

cog_vars <- cog_vars[which(rowMeans(is.na(cog_vars)) < 0.80),]


# log transform RT variables

RT_vars <- cog_vars %>% select(
  CDR_RVIPRTms_v1,TMT_Bms_v1,CDR_SimpleRTms_v1, CDR_ChoiceRTms_v1, 
  CDR_DigitVigRTms_v1, SCB_ConStroopRTms_v1, SCB_InconStroopRTms_v1, 
  Jenson8ChoiceDecisionTimems_v1,Jenson4ChoiceDecisionTimems_v1,
  Jenson1ChoiceDecisionTimems_v1, TMT_Ams_v1) %>% 
  names() # variable names

cog_vars <- cog_vars %>% 
  mutate(across(all_of(RT_vars), ~ log2(.), .names = "log{.col}")) # log 2 transform


# check histograms #

cog_vars[,c(24:34)] %>% select(where(is.numeric)) %>% 
  multi.hist(global = F)

cog_vars %>% select(CDR_ImmWRCorrect_v1, CDR_DelayedWRCorrect_v1,
                               CDR_NumericWMOverallAcc_v1, MMSE,
                               WASImatrixreasoningrawscore, WASI_vocabrawscore,
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
                                 wWASImatrixreasoningrawscore, wWASI_vocabrawscore, 
                                 wCDR_ImmWRCorrect_v1,wCDR_DelayedWRCorrect_v1, 
                                 wCDR_NumericWMOverallAcc_v1, wMMSE,wSerial3sAcc_v1,
                                 wSerial7sAcc_v1) # select relevant variables

# histograms of winsorised variables

cog_vars2 %>% select(starts_with("wlog")) %>% 
  multi.hist(global = F)

cog_vars2 %>% select(wCDR_ImmWRCorrect_v1, wCDR_DelayedWRCorrect_v1,
                    wCDR_NumericWMOverallAcc_v1, wMMSE,
                    wWASImatrixreasoningrawscore, wWASI_vocabrawscore,
                    wSerial3sAcc_v1, wSerial7sAcc_v1) %>% 
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

# rescale WASI scores to per 10 points 

cog_vars2 <- cog_vars2 %>% 
  mutate_at(vars("wWASImatrixreasoningrawscore",
                 "wWASI_vocabrawscore"),
            .funs = ~ . / 10)

# select relevant variables #

cog_vars3 <- cog_vars2 %>% 
  select(IDno, Sex, Initialage, ends_with("_rev"), wCDR_ImmWRCorrect_v1, wWASI_vocabrawscore, wWASImatrixreasoningrawscore,
         wCDR_DelayedWRCorrect_v1, wCDR_NumericWMOverallAcc_v1, wMMSE,
         wSerial3sAcc_v1, wSerial7sAcc_v1)

## single imputation

# predictive mean matching

predictormatrix <- quickpred(cog_vars3, 
                           include=c("Sex", "Initialage"),
                           exclude="IDno",
                           mincor = 0.1)

cog_vars3_imp <- mice(cog_vars3, m = 1, method = "pmm",
                      predictorMatrix = predictormatrix, diagnostics = T,
                      seed = 10)

imp_dat <- complete(cog_vars3_imp) # imputed data 

# check imputation for each variable 

stripplot(cog_vars3_imp, wlogCDR_SimpleRTms_v1_rev +
                  wlogCDR_ChoiceRTms_v1_rev +
                  wlogCDR_DigitVigRTms_v1_rev +
                  wlogJenson8ChoiceDecisionTimems_v1_rev +
                  wlogJenson4ChoiceDecisionTimems_v1_rev +
                  wlogJenson1ChoiceDecisionTimems_v1_rev ~ .imp)

stripplot(cog_vars3_imp, wlogSCB_ConStroopRTms_v1_rev +
            wlogSCB_InconStroopRTms_v1_rev +
            wlogTMT_Ams_v1_rev +
            wlogTMT_Bms_v1_rev  ~ .imp)


stripplot(cog_vars3_imp, wCDR_ImmWRCorrect_v1 +
            wCDR_DelayedWRCorrect_v1 +
            wSerial3sAcc_v1 +
            wSerial7sAcc_v1 +
            wCDR_NumericWMOverallAcc_v1  ~ .imp)

stripplot(cog_vars3_imp, wWASImatrixreasoningrawscore + 
            wMMSE +
            wWASI_vocabrawscore ~ .imp)


# data for CFA 

varTable(imp_dat)

# define CFA model 

model <- 
  
' # model 

Gt =~  wlogCDR_SimpleRTms_v1_rev +
       wlogCDR_ChoiceRTms_v1_rev +
       wlogCDR_DigitVigRTms_v1_rev +
       wlogJenson8ChoiceDecisionTimems_v1_rev +
       wlogJenson4ChoiceDecisionTimems_v1_rev +
       wlogJenson1ChoiceDecisionTimems_v1_rev
       
Gs =~  wlogSCB_ConStroopRTms_v1_rev +
       wlogSCB_InconStroopRTms_v1_rev +
       wlogTMT_Ams_v1_rev +
       # cross loading
       wlogTMT_Bms_v1_rev +
       wlogCDR_ChoiceRTms_v1_rev 

       
Glr_m6 =~ wCDR_ImmWRCorrect_v1 +
       wCDR_DelayedWRCorrect_v1 
       

Gwm_ac =~ wSerial3sAcc_v1 +
        wSerial7sAcc_v1 +
        wCDR_NumericWMOverallAcc_v1 +
        wlogTMT_Bms_v1_rev

Gf =~ wWASImatrixreasoningrawscore

Gc =~ wMMSE +
      wWASI_vocabrawscore


# correlated residuals
wlogTMT_Ams_v1_rev ~~ wlogTMT_Bms_v1_rev
wSerial3sAcc_v1 ~~ wSerial7sAcc_v1
wlogJenson4ChoiceDecisionTimems_v1_rev ~~ wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ wlogJenson1ChoiceDecisionTimems_v1_rev
wlogJenson8ChoiceDecisionTimems_v1_rev ~~ wlogJenson4ChoiceDecisionTimems_v1_rev
'

# fit model

test1 <- cfa(model, data = imp_dat, estimator = "MLR", std.lv = F)
?multi.hist

# results

summary(test1, fit.measures=T, standardized = T)

## Save factor scores ##

preds <- as_tibble(lavPredict(test1, newdata = imp_dat, type = "lv",
                    method = "Bartlett"))


new <- cbind(imp_dat, preds)

## save only factor scores and IDvar to merge back into data

new <- new %>% select(IDno, Gt:Gc)

## merge factor scores into original data

arcli2 <- arcli %>% full_join(new, by = "IDno")

write_sav(arcli2, "arcli with CHC.sav")


