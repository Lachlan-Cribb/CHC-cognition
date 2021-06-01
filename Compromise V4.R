library(haven)
library(tidyverse)
library(lavaan)
library(psych)
library(semTools)
library(semPlot)


## Load data 

arcli <- read_sav("C:\\Users\\lachy\\Dropbox\\ARCLI data entry\\BLINDED DATA SETS\\BLINDED _DO NOT TOUCH\\ARCLI_V1_Baseline_ALL_29-APR-2021.sav")


# remove -999

arcli <- na_if(arcli, -999)

arcli <- na_if(arcli, -99)


# Rename problematic variable names


arcli <- arcli %>% 
  rename(CDR_ImmWRCorrect_v1 = `CDR_ImmWRCorrect#_v1`)

arcli <- arcli %>% 
  rename(CDR_DelayedWRCorrect_v1 = `CDR_DelayedWRCorrect#_v1`)



# create new overall scores where missing


arcli$CDR_PicRecOverallAcc_v1 <- (arcli$CDR_PicRecNovelStimAcc_v1 + arcli$CDR_PicRecOrigStimAcc_v1) / 2

arcli$CDR_Spatial_WM_acc_v1 <- (arcli$CDR_SpatialWMNovelStimAcc_v1 + arcli$CDR_SpatialWMOrigStimAcc_v1) / 2

arcli$CDR_NumericWMOverallAcc_v1 <- (arcli$CDR_NumericWMOrigStimAcc_v1 + arcli$CDR_NumericWMNewStimAcc_v1) / 2




# select cognitive important variables 

cog_vars <- arcli %>% 
  dplyr::select(Sex, IDno, CDR_RVIPRTms_v1, Initialage,CDR_RVIPacc_v1,TMT_Bms_v1, CDR_ImmWRCorrect_v1, CDR_DelayedWRCorrect_v1, 
                CDR_SimpleRTms_v1, CDR_ChoiceRTms_v1, CDR_NumericWMOverallAcc_v1, 
                CDR_DigitVigRTms_v1, SCB_ConStroopRTms_v1, SCB_InconStroopRTms_v1, 
                Jenson8ChoiceDecisionTimems_v1,Jenson4ChoiceDecisionTimems_v1,
                Jenson1ChoiceDecisionTimems_v1,WASImatrixreasoningrawscore, 
                CDR_RVIPacc_v1,WASI_vocabrawscore, MMSE, Serial3sAcc_v1, Serial7sAcc_v1, TMT_Ams_v1)



## Remove rows without cog results


cog_vars <- cog_vars[which(rowMeans(is.na(cog_vars)) < 0.80),]


# check histograms #

#multi.hist(cog_vars[,-1])


# winsorise extreme outliers 


winsorise <- function(x, sd){
  
  maxin <- max(x[scale(x) < sd],na.rm=T)
  
  minin <- min(x[scale(x) > -sd],na.rm=T)
  
  
  x <- ifelse(x > maxin & !is.na(x), maxin, x)
  
  x <- ifelse(x < minin & !is.na(x), minin, x)
  
  x
  
}


# Winsorise outliers 


cog_vars$CDR_SimpleRTms_v1c <- winsorise(cog_vars$CDR_SimpleRTms_v1, sd = 4)

cog_vars$CDR_ChoiceRTms_v1c <- winsorise(cog_vars$CDR_ChoiceRTms_v1, sd = 4)

cog_vars$Jenson4ChoiceDecisionTimems_v1c <- winsorise(cog_vars$Jenson4ChoiceDecisionTimems_v1, sd = 4)

cog_vars$Jenson1ChoiceDecisionTimems_v1c <- winsorise(cog_vars$Jenson1ChoiceDecisionTimems_v1, sd =4)

cog_vars$TMT_Ams_v1c <- winsorise(cog_vars$TMT_Ams_v1, sd = 4)

cog_vars$WASI_vocabrawscorec <- winsorise(cog_vars$WASI_vocabrawscore, sd = 4)

cog_vars$TMT_Bms_v1c <- winsorise(cog_vars$TMT_Bms_v1, sd = 4)

cog_vars$TMT_Bms_v1c <- winsorise(cog_vars$TMT_Bms_v1c, sd = 4)



## remove cases with more than half cog vars missing


cog_vars <- cog_vars[which(rowMeans(is.na(cog_vars)) < 0.50),]

  


# Reverse score reaction time measures

reverse_RT <- function(x){
  x <- x * -1
  x <- x + max(abs(x),na.rm=T) + 1
  x
}


cog_vars <- cog_vars %>% 
  mutate_at(vars("CDR_DigitVigRTms_v1","CDR_RVIPRTms_v1", "CDR_SimpleRTms_v1c",
                 "CDR_ChoiceRTms_v1c","Jenson8ChoiceDecisionTimems_v1",
                 "Jenson4ChoiceDecisionTimems_v1c", 
                 "Jenson1ChoiceDecisionTimems_v1c",
                 "SCB_ConStroopRTms_v1", 
                 "SCB_InconStroopRTms_v1",
                 "TMT_Ams_v1c", "TMT_Bms_v1c"),
            .funs = list(rev = ~ reverse_RT(.)))



## Rescale indicators 

# RT measures to per 25 miliseconds

cog_vars <- cog_vars %>% 
  mutate_at(vars("CDR_RVIPRTms_v1_rev","CDR_DigitVigRTms_v1_rev","CDR_SimpleRTms_v1c_rev",
                 "CDR_ChoiceRTms_v1c_rev","Jenson8ChoiceDecisionTimems_v1_rev",
                 "Jenson4ChoiceDecisionTimems_v1c_rev", "Jenson1ChoiceDecisionTimems_v1c_rev"),
            .funs = ~ . / 25)


# stroop vars to per 50 miliseconds


cog_vars <- cog_vars %>% 
  mutate_at(vars("SCB_ConStroopRTms_v1_rev", "SCB_InconStroopRTms_v1_rev"),
            .funs = ~ . / 50)


# Serial meausre per 10% 


cog_vars <- cog_vars %>% 
  mutate_at(vars("Serial3sAcc_v1","Serial7sAcc_v1"),
            .funs = ~ . / 10)


# Numeric WM acc per 5%

cog_vars$CDR_NumericWMOverallAcc_v1 <- cog_vars$CDR_NumericWMOverallAcc_v1 / 5


# TMT to 5 seconds 

cog_vars$TMT_Ams_v1c_rev <- cog_vars$TMT_Ams_v1c_rev / 5000 


cog_vars$TMT_Bms_v1c_rev <- cog_vars$TMT_Bms_v1c_rev / 5000 


# WASI scores to per 2.5 points 

cog_vars <- cog_vars %>% 
  mutate_at(vars("WASImatrixreasoningrawscore",
                 "WASI_vocabrawscore"),
            .funs = ~ . / 2.5)




# define CFA model 


model <- 
  
' # model 

Gt =~  CDR_SimpleRTms_v1c_rev +
       CDR_ChoiceRTms_v1c_rev +
       CDR_DigitVigRTms_v1_rev +
       Jenson8ChoiceDecisionTimems_v1_rev +
       Jenson4ChoiceDecisionTimems_v1c_rev +
       Jenson1ChoiceDecisionTimems_v1c_rev
       
Gs =~  SCB_ConStroopRTms_v1_rev +
       SCB_InconStroopRTms_v1_rev +
       TMT_Ams_v1c_rev +
       # cross loading
       TMT_Bms_v1c_rev +
       CDR_ChoiceRTms_v1c_rev 

       
Gl_m6 =~ CDR_ImmWRCorrect_v1 +
       CDR_DelayedWRCorrect_v1 
       

Gwm_ac =~ Serial3sAcc_v1 +
        Serial7sAcc_v1 +
        CDR_NumericWMOverallAcc_v1 +
        TMT_Bms_v1c_rev

Gf =~ WASImatrixreasoningrawscore

Gc =~ MMSE +
      WASI_vocabrawscorec


# correlated residuals
TMT_Ams_v1c_rev ~~ TMT_Bms_v1c_rev
Serial3sAcc_v1 ~~ Serial7sAcc_v1
Jenson4ChoiceDecisionTimems_v1c_rev ~~ Jenson1ChoiceDecisionTimems_v1c_rev
Jenson8ChoiceDecisionTimems_v1_rev ~~ Jenson1ChoiceDecisionTimems_v1c_rev
Jenson8ChoiceDecisionTimems_v1_rev ~~ Jenson4ChoiceDecisionTimems_v1c_rev
WASImatrixreasoningrawscore ~~ 0.10*WASImatrixreasoningrawscore
'


# fit model

test1 <- cfa(model, data = cog_vars, missing = "fiml", 
             estimator = "MLR", std.lv = F)


# results

summary(test1, fit.measures=T, standardized = T)



## SEM plot

semPaths(test1, intercepts = F, residuals = F,
         layout = "tree2", nCharNodes = 10, 
         sizeMan = 5, filetype = "png",
         filename = "sem plot CHC", width = 25,
         reorder = F, latents = c("Gt","Gs", "Gwm","Glr","Gf","Gc"))



## Save factor scores ##


preds <- as_tibble(lavPredict(test1, newdata = cog_vars, type = "lv",
                    method = "Bartlett"))


new <- cbind(cog_vars, preds)

## save only factor scores and IDvar to merge back into data

new <- new %>% select(IDno, Gt:Gc)


## merge factor scores into original data

arcli2 <- arcli %>% full_join(new, by = "IDno")




