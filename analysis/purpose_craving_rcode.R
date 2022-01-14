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
df_wide <- read.csv("~/df_wide.csv")
df_wide <- read.csv("~/Box Sync/CurrentProjects_Penn/MURI/Papers/Purpose/analysis/github/df_wide.csv")

df_long <- read.csv("~/df_long.csv")
df_long <- read.csv("~/Box Sync/CurrentProjects_Penn/MURI/Papers/Purpose/analysis/github/df_long.csv")


df_long$condition <- relevel(df_long$condition, ref='control')

df_ema=subset(df_wide, did_ema==1)

df_ema_scan=subset(df_wide, did_ema==1&did_scan==1)

length(unique(df_ema_scan$pID)) #54=number of participants who did both fmri and ema

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


#demography
summary(df_ema_scan$age)
sd(df_ema_scan$age)
table(df_ema_scan$gender)
table(df_ema_scan$race_numeric)
summary(df_ema_scan$race_numeric) #no NAs


###################
##### RESULTS #####
###################

#Alcohol Use Descriptives

#had alcohol
table(df_ema_scan$hadalcohol_total) #out of 54, 11 never had alcohol
(54-11)/54*100 # total percentage who reported having alcohol at least once
by(df_ema_scan$hadalcohol_total, df_ema_scan$gender, table) #0 drinks for 5F 6M 
table(df_ema_scan$gender) #37F 16M (exclude 1 other for percentage calculation to be consistent with previous report)
(16-6)/16*100 #men percentage had alcohol
(37-5)/37*100 #women percentage had alcohol

#Of the ones who reported having drank alcohol at least once
drank_wide = subset(df_wide, drinks_total>0 & did_scan==1)
drank_long=subset(df_long, drinks_total>0 & did_scan==1)
length(unique(drank_long$pID))

summary(drank_wide$drinks_mean) #the average daily amount of consumption
sd(drank_wide$drinks_mean, na.rm=T)
table(drank_long$drinks) #daily drinks range to get the max (min is 1)

drank_wide$drinks_weekly = drank_wide$drinks_mean*7
summary(drank_wide$drinks_weekly)
sd(drank_wide$drinks_weekly, na.rm=T)
table(drank_wide$drinks_weekly)


#############################################################
#### multilevel linear model with standardized variables ####
#############################################################
#standardize variables
df_long$purpose_dailyz <- with(df_long, ave(purpose_daily, pID, FUN=stdz))
df_long$purpose_daily_previousz <- with(df_long, ave(purpose_daily_previous, pID, FUN=stdz))
df_long$craving_previousz <- with(df_long, ave(craving_previous, pID, FUN=stdz))
df_long$drinksz <- with(df_long, ave(drinks, pID, FUN=stdz))

#table 1
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
## (SI6) see stats for vs_cen for purpose*vs two way interaction (reported stats at purpose_low and purpose_high in text) 
## see stats for vs_cen:craving_cen for purpose*vs*craving three-way interaction

test <-lmer(drinks_cen ~ 0 + purpose_low*vs_cen*craving_cen + (0+purpose_low|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_high*vs_cen*craving_cen + (0+purpose_high|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)


#solution for failures to converge
#test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_high*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long)
#lme4::allFit(test)  
#summary(lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa")))

 
#higher ventral striatum cue reactivity was associated with greater average alcohol craving throughout the EMA period 
cor.test(df_wide$craving_mean, df_wide$vs_alc_react)

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

test <-lmer(drinksz ~ 0 + purpose_daily_previousz*vs_alc_react*craving_previousz + (0+purpose_daily_previousz|groupID/pID) + (0+craving_previousz|groupID/pID), df_long)

vif <-data.frame(vif.lme(test)) #low correlations of vif=~1
summary(vif$vif.lme.test.)


#table 1

test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|groupID/pID) + (0+purpose_daily_previousz|groupID/pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)
lme.dscore(test,df_long,type="lme4")

#Figure1
test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+gender+race_numeric+rung_group_cen+condition,df_long,control = lmerControl(optimizer="bobyqa"))
interact_plot(test, pred = craving_cen,  modx=vs_cen, mod2=purpose_cen,interval = TRUE,
              x.label = "craving", y.label = "drinking",
              int.type = "confidence", int.width = 0.95)

#Discussion
#Although the average levels of daily purpose throughout the EMA period were not significantly associated with the neural reactivity within the ventral striatum in our dat

cor.test(df_wide$vs_alc_react, df_wide$purpose_mean)
#cor.test(df_wide$vs_alc_react, df_wide$drinks_total)


#supplementary information (SI)

#SI1.  Race/ethnicity composition of the study samples
summary(df_ema_scan$race) 
# Demographic variables were not associated with any of the main study variables (ps>.10), with an exception that men compared to women showed greater reactivity within the ventral striatum 
v1=df_ema_scan %>%
  dplyr::select(gender_numeric, age, race_numeric, rung_group, purpose_mean, vs_alc_react, drinks_mean)


M <- Hmisc::rcorr(as.matrix(v1))
library(corrplot)
corrplot(M$r, p.mat = M$P, insig = "label_sig",
         sig.level = c(.001, .05, .10), pch.cex=0.5, pch.col = "white",
         method="color", type="lower")
#gender is associated with vs (male with higher vs cue reactivity)
cor.test(df_ema_scan$gender_numeric, df_ema_scan$vs_alc_react)



#SI3. Neurosynth map of ‘craving’

test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*neurosynth_craving_alc_react + (0+craving_previousz|groupID/pID) + (0+purpose_daily_previousz|groupID/pID)+age+as.factor(gender)+as.factor(race_numeric)+rung_group+as.factor(condition),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)
lme.dscore(test,df_long,type="lme4")

df_long$neurosynth_craving_alc_react_cen = df_long$neurosynth_craving_alc_react - mean(df_long$neurosynth_craving_alc_react,na.rm=T)

#followup simple slopes analyses
test <-lmer(drinks_cen ~ 0 + purpose_low*neurosynth_craving_alc_react_cen*craving_cen  + (0+purpose_low|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*neurosynth_craving_alc_react_cen*craving_cen  + (0+purpose_low|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_high*neurosynth_craving_alc_react_cen*craving_cen  + (0+purpose_low|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+as.factor(gender)+as.factor(race_numeric)+rung_group_cen+as.factor(condition),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)


#Figure SI.
test <-lmer(drinks_cen ~ 0 + purpose_cen*neurosynth_craving_alc_react_cen*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID)+age_cen+gender+race_numeric+rung_group_cen+condition,df_long,control = lmerControl(optimizer="bobyqa"))
interact_plot(test, pred = craving_cen,  modx=neurosynth_craving_alc_react_cen, mod2=purpose_cen,interval = TRUE,
              x.label = "craving", y.label = "drinking",
              int.type = "confidence", int.width = 0.95)


#SI5. Results not controlling for potential covariates

test <-lmer(drinksz ~ 0 + craving_previousz*purpose_daily_previousz*vs_alc_react + (0+craving_previousz|groupID/pID) + (0+purpose_daily_previousz|groupID/pID),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)
lmer.beta(test)


test <-lmer(drinks_cen ~ 0 + purpose_low*vs_cen*craving_cen + (0+purpose_low|groupID/pID) + (0+craving_cen|groupID/pID),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_cen*vs_cen*craving_cen + (0+purpose_cen|groupID/pID) + (0+craving_cen|groupID/pID),df_long,control = lmerControl(optimizer="bobyqa"))
summ(test, digits = 3, confint = TRUE)
summary(test)

test <-lmer(drinks_cen ~ 0 + purpose_high*vs_cen*craving_cen + (0+purpose_high|groupID/pID) + (0+craving_cen|groupID/pID),df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)




