library(DOS2)
library(VIM)
library(mice)
library(MASS)
library(plyr)
library(dplyr)
library(tidyr)
library(pastecs)
library(mlbench)
library(RItools)
library(ggplot2)
library(optmatch)
library(rcbalance)
library(MatchIt)
library(naniar) # visualization of missing data
library(tidyverse)
library(sensitivitymult)
library(sensitivitymv)


################################################
##              helper functions              ##
################################################

# function to summarize matches 
summarize.match <- function(dat) {
  match.df <- matchit(z ~ Pressure.height + Wind.speed + Humidity + Tem.Sanburg + Base.height + Pressure.gradient+ 
                          Visibility,data = dat,, distance = "mahalanobis")

  c <- match.data(match.df,'control') 
  t <- match.data(match.df,'treated') 
  pair_df <- merge(c, t, by = "subclass", all.x = TRUE,all.y = TRUE) 
  pair_df <- pair_df %>% mutate(effect = maxO3.y - maxO3.x)
  DIM <- mean(pair_df$effect)
  VAR <- var(pair_df$effect)
  r = c(DIM, VAR)
  return(r)
}

summarize.simple <- function(df) {
  df_0 <- df[df$z == 0,]
  df_1 <- df[df$z == 1,]
  DIM <- mean(df_1$maxO3) - mean(df_0$maxO3)
  VAR <- var(df_1$maxO3)/nrow(df_1) + var(df_0$maxO3)/nrow(df_0)
  r = c(DIM, VAR)
  return(r)
}

get.result <- function(DIMs,VARs) {
  DIM_bar = mean(DIMs)
  VAR_bar = mean(VARs)
  VAR_var = var(DIMs)
  VAR_all <- VAR_bar*5/4*VAR_var
  r = c(DIM_bar, VAR_bar,VAR_var,VAR_all)
  return(r)
}

add_z <- function(oze,df){
  threshold = mean(oze$Base.tem, na.rm=TRUE)
  df$z <- rep(NA, nrow(df))
  df$z <- ifelse(df$Base.tem <= threshold, 0, 1)
  df <- select(df, select = -Base.tem)
  return(df)
}

#########################################
########Importing & Cleaning Data########
#########################################
data(Ozone)
attach(Ozone)
# Preprocess data set
oze_all <- dplyr::rename(Ozone,c(Month = V1,Day.month = V2,Day.week = V3,maxO3 = V4, Pressure.height = V5, Wind.speed = V6,Humidity = V7,Tem.Sanburg = V8,
         Tem.EM = V9, Base.height = V10, Pressure.gradient = V11, Base.tem = V12, Visibility = V13));
oze <- subset(oze_all, select = -c(Month,Day.month,Day.week) )
res<-summary(aggr(oze, sortVar=TRUE))$combinations # check the missing data distribution
oze = oze[!is.na(oze$maxO3),] # remove rows without maxO3
oze <- subset(oze, select = -c(Tem.EM) ) # remove Tem.EM because of 40% missing data

# Anlaysis on missing data pattern
summary(oze)
oze_state <- data.frame(stat.desc(oze))
mu <- colMeans(oze, na.rm = TRUE)
cv <- cov(oze, use = "pairwise")

md.pattern(oze)
gg_miss_var(oze)
vis_miss(oze, sort_miss = TRUE)

##########################
########Imputation########
##########################
set.seed(1)
# Stochastic Regression Imputation
sto_reg_imp <- mice(oze, method = "norm.nob", m = 5, maxit = 5, print = FALSE)
densityplot(sto_reg_imp)
sto_reg_imp$imp$Humidity


# MICE - pmm(predictive mean matching) a Markov Chain Monte Carlo (MCMC) method
pmm.imp.all <- mice(oze, m=5,seed=1)
pmm1 <- mice::complete(pmm.imp.all,1)
pmm2 <- mice::complete(pmm.imp.all,2)
pmm3 <- mice::complete(pmm.imp.all,3)
pmm4 <- mice::complete(pmm.imp.all,4)
pmm5 <- mice::complete(pmm.imp.all,5)

#pmm.imp.all <- vector("list", 5)
#for (i in 1:5){
  #imp <- mice(oze, maxit = 5, seed = i, print = FALSE)
  #pmm.imp.all[[i]] <- imp
#}
pmm1 <- mice::complete(pmm.imp.all[[1]])
pmm2 <- mice::complete(pmm.imp.all[[2]])
pmm3 <- mice::complete(pmm.imp.all[[3]])
pmm4 <- mice::complete(pmm.imp.all[[4]])
pmm5 <- mice::complete(pmm.imp.all[[5]])

summary(pmm.imp.all)
#pmm_oze$imp$Pressure.height
xyplot(pmm.imp.all, maxO3 ~ Tem.EM | .imp, pch = 20, cex = 1.4)
densityplot(pmm.imp.all)
pmm_fit <- with(data = pmm.imp.all, exp = lm(maxO3 ~ Pressure.height + Wind.speed + Humidity + Tem.Sanburg + Base.height + 
                                          Pressure.gradient+ Base.tem + Visibility ))
pmm_combine <- pool(pmm_fit)
summary(pmm_combine)


############################################
######## Calculate treatment effect ########
############################################

# get each imputated table
sr1 <- mice::complete(sto_reg_imp, 1)
sr2 <- mice::complete(sto_reg_imp, 2)
sr3 <- mice::complete(sto_reg_imp, 3)
sr4 <- mice::complete(sto_reg_imp, 4)
sr5 <- mice::complete(sto_reg_imp, 5)
pmm1 <- mice::complete(pmm.imp.all, 1)
pmm2 <- mice::complete(pmm.imp.all, 2)
pmm3 <- mice::complete(pmm.imp.all, 3)
pmm4 <- mice::complete(pmm.imp.all, 4)
pmm5 <- mice::complete(pmm.imp.all, 5)
# add Z column
sr1_z <- add_z(oze,sr1)
sr2_z <- add_z(oze,sr2)
sr3_z <- add_z(oze,sr3)
sr4_z <- add_z(oze,sr4)
sr5_z <- add_z(oze,sr5)
pmm1_z <- add_z(oze,pmm1)
pmm2_z <- add_z(oze,pmm2)
pmm3_z <- add_z(oze,pmm3)
pmm4_z <- add_z(oze,pmm4)
pmm5_z <- add_z(oze,pmm5)
summary(pmm1_z)
summary(pmm2_z)
summary(pmm3_z)
summary(oze)

##########################
######## Matching ########
##########################
dev.off()
### Stochastic Regression Imputation###
#ms.sr1 <- pairmatch(z ~ Pressure.height + Wind.speed + Humidity + Tem.Sanburg +
                   #Base.height + Pressure.gradient + Visibility, data=sr1_z, remove.unmatchables = TRUE)
r_sr1 <- summarize.match(sr1_z)
DIM_srm1 <-r_sr1[1]
VAR_srm1 <- r_sr1[2]
r_sr2 <- summarize.match(sr2_z)
DIM_srm2 <-r_sr2[1]
VAR_srm2 <- r_sr2[2]
r_sr3 <- summarize.match(sr3_z)
DIM_srm3 <-r_sr3[1]
VAR_srm3 <- r_sr3[2]
r_sr4 <- summarize.match(sr4_z)
DIM_srm4 <-r_sr4[1]
VAR_srm4 <- r_sr4[2]
r_sr5 <- summarize.match(sr5_z)
DIM_srm5 <-r_sr5[1]
VAR_srm5 <- r_sr5[2]
DIMs_sr_match <- c(DIM_srm1,DIM_srm2,DIM_srm3,DIM_srm4,DIM_srm5)
VARs_sr_match <- c(VAR_srm1,VAR_srm2,VAR_srm3,VAR_srm4,VAR_srm5)
sr_result_match <- get.result(DIMs_sr_match,VARs_sr_match)

# add propensity score
#sr1_z$prop <- glm(z ~  Wind.speed + Humidity + Pressure.gradient+ Visibility, family=binomial, data=sr1_z)$fitted.values
#ggplot(data=sr1_z, aes(x=prop, group=as.factor(z), fill=as.factor(z))) + geom_density(alpha=0.5) + theme_bw()
#plot(xBalance(z ~ Pressure.height + Wind.speed + Humidity + Tem.Sanburg + Base.height + Pressure.gradient+ 
                #Visibility + prop + strata(ms.sr2)-1, data = sr1_z))

#### MICE pmm ####
r_pmm1 <- summarize.match(pmm1_z)
DIM_pmm1 <-r_pmm1[1]
VAR_pmm1 <- r_pmm1[2]
r_pmm2 <- summarize.match(pmm2_z)
DIM_pmm2 <-r_pmm2[1]
VAR_pmm2 <- r_pmm2[2]
r_pmm3 <- summarize.match(pmm3_z)
DIM_pmm3 <-r_pmm3[1]
VAR_pmm3 <- r_pmm3[2]
r_pmm4 <- summarize.match(pmm4_z)
DIM_pmm4 <-r_pmm4[1]
VAR_pmm4 <- r_pmm4[2]
r_pmm5 <- summarize.match(pmm5_z)
DIM_pmm5 <-r_pmm5[1]
VAR_pmm5 <- r_pmm5[2]
DIMs_pmm_match <- c(DIM_pmm1,DIM_pmm2,DIM_pmm3,DIM_pmm4,DIM_pmm5)
VARs_pmm_match <- c(VAR_pmm1,VAR_pmm2,VAR_pmm3,VAR_pmm4,VAR_pmm5)
pmm_result_match <- get.result(DIMs_pmm_match,VARs_pmm_match)


##########################
######## Basic DIM #######
##########################
rs_sr1 <- summarize.simple(sr1_z)
DIM_sr1 <-rs_sr1[1]
VAR_sr1 <- rs_sr1[2]
rs_sr2 <- summarize.simple(sr2_z)
DIM_sr2 <-rs_sr2[1]
VAR_sr2 <- rs_sr2[2]
rs_sr3 <- summarize.simple(sr3_z)
DIM_sr3 <-rs_sr3[1]
VAR_sr3 <- rs_sr3[2]
rs_sr4 <- summarize.simple(sr4_z)
DIM_sr4 <-rs_sr4[1]
VAR_sr4 <- rs_sr4[2]
rs_sr5 <- summarize.simple(sr5_z)
DIM_sr5 <-rs_sr5[1]
VAR_sr5 <- rs_sr5[2]
DIMs_sr_basic <- c(DIM_sr1,DIM_sr2,DIM_sr3,DIM_sr4,DIM_sr5)
VARs_sr_basic <- c(VAR_sr1,VAR_sr2,VAR_sr3,VAR_sr4,VAR_sr5)
sr_result_basic <- get.result(DIMs_sr_basic,VARs_sr_basic)

rs_pmm1 <- summarize.simple(pmm1_z)
DIMs_pmm1 <-rs_pmm1[1]
VARs_pmm1 <- rs_pmm1[2]
rs_pmm2 <- summarize.simple(pmm2_z)
DIMs_pmm2 <-rs_pmm2[1]
VARs_pmm2 <- rs_pmm2[2]
rs_pmm3 <- summarize.simple(pmm3_z)
DIMs_pmm3 <-rs_pmm3[1]
VARs_pmm3 <- rs_pmm3[2]
rs_pmm4 <- summarize.simple(pmm4_z)
DIMs_pmm4 <-rs_pmm4[1]
VARs_pmm4 <- rs_pmm4[2]
rs_pmm5 <- summarize.simple(pmm5_z)
DIMs_pmm5 <-rs_pmm5[1]
VARs_pmm5 <- rs_pmm5[2]
DIMs_pmm_basic <- c(DIMs_pmm1,DIMs_pmm2,DIMs_pmm3,DIMs_pmm4,DIMs_pmm5)
VARs_pmm_basic <- c(VARs_pmm1,VARs_pmm2,VARs_pmm3,VARs_pmm4,VARs_pmm5)
pmm_result_basic <- get.result(DIMs_pmm_basic,VARs_pmm_basic)
