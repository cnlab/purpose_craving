library(mediation)
library(dplyr)
library(lme4)
library(Hmisc)
library(lmerTest) # useful if you want to get the p values for the models of each path. 
#detach("package:lmerTest", unload = TRUE) # detach and re-run models prior to attempting mediation.
library(interactions)
library(jtools)
library(ggsci)
library(weights)
library(EMAtools)

rm(list = ls())

#load data
df_wide <- read.csv("~/df_wide_10_22_21.csv")
df_long <- read.csv("~/df_long_10_22_21.csv")

df_long$condition <- relevel(df_long$condition, ref='control')

length(unique(df_long$pID)) #279 who did covid ema
table(df_wide$did_ema) #double check that the N=279 in df_wide as well

length(unique(df_long$vs_alc_react)) #54 includes NA, so 53 with both the EMA and brain data
df_ema=subset(df_wide, did_ema==1)
df_ema_scan=subset(df_wide, did_ema==1&did_scan==1)
length(unique(df_ema_scan$pID)) #53=number of participants who did both fmri and covid ema

#define a function to get standardized beta from lmer
lmer.beta <- function(object) {
  sdy <- sd(attr(object, "resp")$y) # the y values are now in the 'y' slot 
  ###                                 of the resp attribute
  sdx <- apply(attr(object, "pp")$X, 2, sd) # And the X matriz is in the 'X' slot of the pp attr
  sc <- fixef(object)*sdx/sdy
  #mimic se.ranef from pacakge "arm"
  se.fixef <- function(obj) as.data.frame(summary(obj)[10])[,2] # last change - extracting 
  ##             the standard errors from the summary
  se <- se.fixef(object)*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

###################
##### METHODS #####
###################

#Abstract; number of data points

NROW(df_long$Purpose) - sum(is.na(df_long$Purpose))
#6451
NROW(df_long$drinks) - sum(is.na(df_long$drinks))
#12624
NROW(df_long$Craving_Alc) - sum(is.na(df_long$Craving_Alc))
#12592
6451 +12624 +12592
#31667

# ema subjects (n=279)
summary(df_ema$age)
sd(df_ema$age, na.rm=T)
table(df_ema$gender)
table(df_ema$race_numeric)
summary(df_ema$race_numeric) #25 NAs

#race_ema=data.frame(table(df_ema$race))
#write.csv(race_ema, "~/Desktop/race_ema.csv")

#fmri subjects who did ema (n=53)
summary(df_ema_scan$age)
sd(df_ema_scan$age)
table(df_ema_scan$gender)
table(df_ema_scan$race_numeric)
summary(df_ema_scan$race_numeric) #no NAs

#race_scan=data.frame(table(df_ema_scan$race))
#write.csv(race_scan, "~/Desktop/race_scan.csv")

#check is demographics are associated with main study variables
v1=df_ema_scan %>%
  dplyr::select(gender_numeric, age, race_numeric, rung_group, purpose_mean, vs_alc_react, drinks_mean)


M <- Hmisc::rcorr(as.matrix(v1))
library(corrplot)
corrplot(M$r, p.mat = M$P, insig = "label_sig",
         sig.level = c(.001, .05, .10), pch.cex=0.5, pch.col = "white",
         method="color", type="lower")
#gender is associated with vs (male with higher vs cue reactivity)
cor.test(df_ema_scan$gender_numeric, df_ema_scan$vs_alc_react)

###################
##### RESULTS #####
###################

#Alcohol Use Descriptives

#had alcohol
table(df_ema$hadalcohol_total) #out of 279, 90 never had alcohol
by(df_wide$hadalcohol_tota, df_wide$gender, table) #0 drinks for 59F 22M (9 no gender info) 
table(df_ema$gender) #190F 63M 
(279-90)/279*100 # total percentage who reported having alcohol at least once
(63-22)/63*100 #men percentage had alcohol
(190-59)/190*100 #women percentage had alcohol

#Of the ones who reported having drank alcohol at least once
drank_wide = subset(df_wide, drinks_total>0)
drank_long=subset(df_long, drinks_total>0)

summary(drank_wide$drinks_mean) #the average daily amount of consumption
sd(drank_wide$drinks_mean, na.rm=T)
table(drank_long$drinks) #daily drinks range to get the max (min is 1)


#############################################################
#### multilevel linear model with standardized variables ####
#############################################################
#standardize variables
df_long$purpose_dailyz <- with(df_long, ave(purpose_daily, pID, FUN=stdz))
df_long$purpose_daily_previousz <- with(df_long, ave(purpose_daily_previous, pID, FUN=stdz))
df_long$craving_previousz <- with(df_long, ave(craving_previous, pID, FUN=stdz))
df_long$drinksz <- with(df_long, ave(drinks, pID, FUN=stdz))

#Alcohol Craving and Subsequent Alcohol Consumption
#model 1
# We found that the more people craved alcohol, the larger amount of alcohol they subsequently ended up consuming
test <-lmer(drinksz~ 0 + craving_previousz +  (0+craving_previousz|groupID/pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summary(test)
summ(test)
confint(test)
lmer.beta(test)
lme.dscore(test,df_long,type="lme4")

#Interactions among Daily Purpose, Neural Alcohol Cue Reactivity, and Alcohol Craving on Subsequent Alcohol Consumption

#higher ventral striatum cue reactivity was associated with greater average alcohol craving throughout the EMA period 
cor.test(df_wide$craving_mean, df_wide$vs_alc_react)

#However, analysis of the Variance Inflation Factor (VIF) showed that the model did not have problems of multicollinearity
#check vif for the amount of multicollinearity

vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }

test <-lmer(drinksz ~ 0 + purpose_daily_previousz*vs_alc_react*craving_previousz + (0+purpose_daily_previousz|pID) + (0+craving_previousz|pID), df_long)

vif <-data.frame(vif.lme(test)) #low correlations of vif=~1
summary(vif$vif.lme.test.)


#model 2
#There was a significant two-way interaction between daily purpose in life and ventral striatum activity predicting subsequent alcohol use
test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|groupID/pID) + (0+purpose_daily_previousz|groupID/pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)
lme.dscore(test,df_long,type="lme4")


#Follow-up simple slopes analyses 

#within-person mean-center variables
df_long$purpose_cen =df_long$purpose_daily_previous - df_long$purpose_mean
df_long$craving_cen =df_long$craving_previous - df_long$craving_mean
df_long$drinks_cen =df_long$drinks - df_long$drinks_mean
df_long$vs_cen = df_long$vs_alc_react - mean(df_long$vs_alc_react,na.rm=T)
df_long$age_cen = df_long$age - mean(df_long$age,na.rm=T)
df_long$rung_group_cen = df_long$rung_group - mean(df_long$rung_group,na.rm=T) 


df_long=df_long %>%
  group_by(pID) %>%
  dplyr::mutate(purpose_cen_sd = sd(purpose_cen, na.rm=T))

                

df_long$purpose_low = df_long$purpose_cen + df_long$purpose_cen_sd
df_long$purpose_high = df_long$purpose_cen - df_long$purpose_cen_sd

df_long$vs_low <- df_long$vs_cen + sd(df_long$vs_cen, na.rm=T)
df_long$vs_high <- df_long$vs_cen - sd(df_long$vs_cen, na.rm=T)

#simple slopes analysis
## see stats for vs_cen for purpose*vs two way interaction (reported stats at purpose_low and purpose_high in text)
## see stats for vs_cen:craving_cen for purpose*vs*craving three-way interaction

test <-lmer(drinks_cen ~ 0 + purpose_low*vs_cen*craving_cen + (0+purpose_low|pID) + (0+craving_cen|pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_high*vs_cen*craving_cen + (0+purpose_high|pID) + (0+craving_cen|pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)


#a solution for failures to converge)
#test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_high*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
#lme4::allFit(test)  
#summary(lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa")))

 
#table 1
#M1
test <-lmer(drinksz~ 0 + craving_previousz +  (0+craving_previousz|pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summary(test)
summ(test)
confint(test)
lmer.beta(test)

#M2
test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|pID) + (0+purpose_daily_previousz|pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)

#Figure1
test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID)+age_cen+gender+race_numeric+rung_group_cen+condition,df_long,control = lmerControl(optimizer="bobyqa"))
interact_plot(test, pred = craving_cen,  modx=vs_cen, mod2=purpose_cen,interval = TRUE,
              x.label = "craving", y.label = "drinking",
              int.type = "confidence", int.width = 0.95)

#Discussion
#Although the average levels of daily purpose throughout the EMA period were not significantly associated with the neural reactivity within the ventral striatum in our dat

cor.test(df_wide$vs_alc_react, df_wide$purpose_mean)
cor.test(df_wide$vs_alc_react, df_wide$drinks_total)


#supplementary information (SI)

#SI1.  Race/ethnicity composition of the study samples

##sample 1
df_scan=subset(df_wide, did_scan==1)
summary(df_scan$race)

##sample 2
summary(df_ema$race) 

#sample 1 and 2
summary(df_ema_scan$race) 


#SI2. Results not controlling for potential covariates


#model 1
test <-lmer(drinksz~ 0 + craving_previousz +  (0+craving_previousz|groupID/pID),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)
lme.dscore(test,df_long,type="lme4")

#model 2
#test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|pID) + (0+purpose_daily_previousz|pID),df_long)
test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|groupID/pID) + (0+purpose_daily_previousz|groupID/pID),df_long)

summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)


test <-lmer(drinks_cen ~ 0 + purpose_low*vs_cen*craving_cen + (0+purpose_low|pID) + (0+craving_cen|pID),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_high*vs_cen*craving_cen + (0+purpose_high|pID) + (0+craving_cen|pID),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)


#SI4. Neurosynth map of ‘craving’


#model 2
test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*neurosynth_craving_alc_react + (0+craving_previousz|pID) + (0+purpose_daily_previousz|pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)


df_long$neurosynth_craving_alc_react_cen = df_long$neurosynth_craving_alc_react - mean(df_long$neurosynth_craving_alc_react,na.rm=T)

#followup simple slopes analyses
test <-lmer(drinks_cen ~ 0 + purpose_high*neurosynth_craving_alc_react_cen + (0+purpose_low|pID) +age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*neurosynth_craving_alc_react_cen + (0+purpose_low|pID) +age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_low*neurosynth_craving_alc_react_cen + (0+purpose_low|pID) +age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)


#Figure SI.
test <-lmer(drinks_cen ~ 0 + purpose_cen*neurosynth_craving_alc_react_cen*craving_cen + (0+purpose_cen|pID) + (0+craving_cen|pID)+age_cen+gender+race_numeric+rung_group_cen+condition,df_long,control = lmerControl(optimizer="bobyqa"))
interact_plot(test, pred = craving_cen,  modx=neurosynth_craving_alc_react_cen, mod2=purpose_cen,interval = TRUE,
              x.label = "craving", y.label = "drinking",
              int.type = "confidence", int.width = 0.95)

