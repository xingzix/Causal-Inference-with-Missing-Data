library(mlbench)
library(VIM)
library(naniar)
library(missMDA)
library(Amelia)
library(mice)
library(missForest)
library(FactoMineR)
library(tidyverse)
library(corrplot)
library(DOS2)
library(optmatch)
library(RItools)
library(plyr)
library(rcbalance)

##############################
#Part I. Importing Data#######
##############################

data(Ozone)
attach(Ozone)
oze_all <- dplyr::rename(Ozone,c(Month = V1,Day.month = V2,Day.week = V3,maxO3 = V4, Pressure.height = V5, Wind.speed = V6,Humidity = V7,Tem.Sanburg = V8,
                                 Tem.EM = V9, Base.height = V10, Pressure.gradient = V11, Base.tem = V12, Visibility = V13))
oze <- subset(oze_all, select = -c(Month,Day.month,Day.week) )



dim(oze)
summary(oze)
head(oze)
mu <- colMeans(oze, na.rm = TRUE)
cv <- cov(oze, use = "pairwise")

#How is the missing data?
cat('Deleting observations with missing data for ozone data leads to a table with', dim(na.omit(oze))[1], 'rows')
cat('The percentage of the missing data in the original dataframe is', sum(is.na(oze))/prod(dim(oze)))
md.pattern(oze, rotate.names = TRUE)
gg_miss_var(oze)
res<-summary(aggr(oze, sortVar=TRUE))$combinations
head(res[rev(order(res[,2])),])

par(mfrow=c(1,3))
marginplot(oze[,c("Humidity","maxO3")])
marginplot(oze[,c("Base.height","maxO3")])
marginplot(oze[,c("Base.tem","maxO3")])
vis_miss(oze, sort_miss = TRUE) 

#What is the correlation between max03 and other variables?
res2 <- cor(oze, use = "complete.obs")
corrplot(res2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)



## Missing data exploratory analysis 

pct_miss(oze) # percentage of missing value in the data.
n_complete(oze) # without missing value
n_miss(oze$maxO3) # number of missing values in maxO3
# There are only 5 rows with maxO3 missing - delete these 5 rows


oze <- subset(oze, select = -c(Tem.EM) )
oze = oze[!is.na(oze$maxO3),]

dim(oze)
as_shadow(oze)
bind_shadow(oze)

library(dplyr)
oze %>%
  miss_var_summary()

oze %>%
  bind_shadow() %>%
  group_by(maxO3_NA) %>%
  summarise_at(.vars = "Tem.EM",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

oze %>%
  bind_shadow() %>%
  group_by(maxO3_NA) %>%
  summarise_at(.vars = "Humidity",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

oze %>%
  bind_shadow() %>%
  group_by(maxO3_NA) %>%
  summarise_at(.vars = "Base.height",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

oze %>%
  bind_shadow() %>%
  group_by(maxO3_NA) %>%
  summarise_at(.vars = "Base.tem",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)

oze %>%
  bind_shadow() %>%
  group_by(maxO3_NA) %>%
  summarise_at(.vars = "Pressure.height",
               .funs = c("mean", "sd", "var", "min", "max"),
               na.rm = TRUE)


#############################
#Part II. Imputing X ########
#############################

## Method 1: PCA
# estimate the number of components from incomplete data
nb <- estim_ncpPCA(oze ,method.cv = "Kfold", verbose = FALSE) 
nb$ncp #5
plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# iterativePCA algorithm
res.comp <- imputePCA(oze, ncp = nb$ncp) 
# a glimpse into the imputed dataset
res.comp$completeObs[1:5,] 

imp.pca <- res.comp$completeObs
res.pca <- PCA(imp.pca, ncp = nb$ncp, graph=TRUE)
# plot(res.pca, hab=12)
# plot(res.pca, choix="var")
imp.pca <- as.data.frame(imp.pca)
# Compare density of original and imputed Base.tem
par(mfrow=c(2,2))
plot(density(imp.pca$Pressure.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Pressure.height)),col='blue',lwd = 2, ann=FALSE)
mtext("Pressure.height", side=3, line=0)

plot(density(imp.pca$Humidity),col='red',ann=FALSE)
lines(density(na.omit(oze$Humidity)),col='blue',lwd = 2, ann=FALSE)
mtext("Humidity", side=3, line=0)

plot(density(imp.pca$Base.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.height)),col='blue',lwd = 2, ann=FALSE)
mtext("Base.height", side=3, line=0)

plot(density(imp.pca$Base.tem),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.tem)),col='blue',lwd = 2, ann=FALSE)
mtext("Base.tem", side=3, line=0)


## Method 2: Multiple Imputation
#assume a PCA based model - We use the R package missMDA
#Now generate 100 imputed data sets with the MIPCA method and 5 components.
#Store the result in a variable called res.MIPCA.
res.MIPCA <- MIPCA(oze, ncp = 5, nboot = 5) # MI with PCA using 5 dimensions 
res.MIPCA <- lapply(res.MIPCA$res.MI, as.data.frame)
# res.MIPCA[[1]] # to get the first imputed dataset

#Inspect imputed data
# plot(res.MIPCA,choice= "ind.supp")
# plot(res.MIPCA,choice= "var")

# Compare density of original and imputed Base.tem
par(mfrow=c(2,2))
plot(density(na.omit(oze$Pressure.height)),col='blue', lwd = 2, ann=FALSE)
lines(density(res.MIPCA[[1]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Pressure.height),col='red',ann=FALSE)
mtext("Pressure.height", side=3, line=0)

plot(density(na.omit(oze$Humidity)),col='blue', lwd = 2,ann=FALSE)
lines(density(res.MIPCA[[1]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Humidity),col='red',ann=FALSE)
mtext("Humidity", side=3, line=0)

plot(density(na.omit(oze$Base.height)),col='blue', lwd = 2,ann=FALSE)
lines(density(res.MIPCA[[1]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Base.height),col='red',ann=FALSE)
mtext("Base.height", side=3, line=0)

plot(density(na.omit(oze$Base.tem)),col='blue',lwd = 2, ann=FALSE)
lines(density(res.MIPCA[[1]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Base.tem),col='red',ann=FALSE)
mtext("Base.tem", side=3, line=0)




## Method 3: The naive approach: Delete all the missing data
oze_dropped <- oze %>% drop_na()

# Compare density of original and imputed Base.tem
par(mfrow=c(2,2))
plot(density(oze_dropped$Pressure.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Pressure.height)),col='blue',lwd = 2, ann=FALSE)
mtext("Pressure.height", side=3, line=0)

plot(density(oze_dropped$Humidity),col='red',ann=FALSE)
lines(density(na.omit(oze$Humidity)),col='blue', lwd = 2,ann=FALSE)
mtext("Humidity", side=3, line=0)

plot(density(oze_dropped$Base.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.height)),col='blue',lwd = 2, ann=FALSE)
mtext("Base.height", side=3, line=0)

plot(density(oze_dropped$Base.tem),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.tem)),col='blue',lwd = 2, ann=FALSE)
mtext("Base.tem", side=3, line=0)



## Method 4: Mean Imputation
mean_imp <- mice(oze, method = "mean", m = 1, maxit = 1)
densityplot(mean_imp)

## Method 4:Linear Regression
linear_oze <- data.frame(lapply(oze, function(X) na_interpolation(X)))
summary(linear_oze)
ggplot( bind_shadow(linear_oze),aes(x = Tem.EM,fill = oze)) + geom_density(alpha=0.5)

## Method 5:Stochastic Regression Imputation
sto_reg_imp <- mice(oze, method = "norm.nob", m = 5, maxit = 5, print = FALSE)
densityplot(sto_reg_imp)

## Method 6: Linear RegressionMICE - pmm
pmm_imp <- mice(oze, m=5)
summary(pmm_imp)
#pmm_oze$imp$Pressure.height
xyplot(pmm_imp, maxO3 ~ Tem.EM | .imp, pch = 20, cex = 1.4)
densityplot(pmm_imp)


pmm_fit <- with(data = pmm_imp, exp = lm(maxO3 ~ Pressure.height + Wind.speed + Humidity + Tem.Sanburg +Tem.EM + Base.height + 
                                           Pressure.gradient+ Base.tem + Visibility ))

pmm_combine <- pool(pmm_fit)
summary(pmm_combine)

## Method 7: Random Forest
# NRMSE: normalized mean squared error -> error derived from imputing continuous values
# PFC: proportion of falsely classified -> error derived from imputing categorical values
forest <- missForest(xmis = oze, maxiter = 20, ntree = 100)
oze.tree <-forest$ximp
oze.tree.error <- forest$OOBerror
oze.tree <- as.data.frame(oze.tree)

# Compare density of original and imputed Base.tem
par(mfrow=c(2,2))
plot(density(oze.tree$Pressure.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Pressure.height)),col='blue',lwd = 2, ann=FALSE)
mtext("Pressure.height", side=3, line=0)

plot(density(oze.tree$Humidity),col='red',ann=FALSE)
lines(density(na.omit(oze$Humidity)),col='blue', lwd = 2,ann=FALSE)
mtext("Humidity", side=3, line=0)

plot(density(oze.tree$Base.height),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.height)),col='blue', lwd = 2,ann=FALSE)
mtext("Base.height", side=3, line=0)

plot(density(oze.tree$Base.tem),col='red',ann=FALSE)
lines(density(na.omit(oze$Base.tem)),col='blue', lwd = 2,ann=FALSE)
mtext("Base.tem", side=3, line=0)


# Now we have all the imputed data sets. 
# We need to assign z on Base.tem for each data set and then delete Base.tem in the imputed dataset

#######################################################
############## function for adding z  #################

#df is an imputed dataset and oze is the original dataset
addz_IMdata <- function(df, oze){
  threshold = mean(oze$Base.tem, na.rm=TRUE)
  df$z <- rep(NA, nrow(df))
  df[df$Base.tem <= threshold, ][, "z"] <- 0
  df[df$Base.tem > threshold, ][, "z"] <- 1
}

################################################
##              helper functions              ##
################################################

# function to summarize matches 
summarize.match <- function(dat, ms, ps.name="prop") {
  adat <- dat
  adat$pair <- ms
  adat <- adat[!is.na(adat$pair),]
  adat.treat <- adat[adat$z==1, ]
  adat.ctrl <- adat[adat$z==0, ]
  
  adat.m <- merge(adat.treat, adat.ctrl, by="pair", suffixes=c(".1", ".0"))
  adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0", "pair"))]
  adat.m <- adat.m[, sort(names(adat.m), index.return=TRUE)$ix]
  
  p0.name <- paste0(ps.name,".", 0)
  p1.name <- paste0(ps.name,".",1)
  
  adat.m.tmp.1 <- adat.m[, -which(names(adat.m) %in% c(p0.name, p1.name))]
  adat.m.tmp.2 <- adat.m[, c(p0.name, p1.name)]
  
  adat.m <- cbind(adat.m.tmp.1, adat.m.tmp.2)
  
  return(adat.m)
}


ms.transform <- function(dat.arg, ms.rcbal) {
  ctrl <- seq(sum(dat.arg$z==0))
  matched.ctrl <- ms.rcbal
  unmatched.ctrl <- setdiff(ctrl,ms.rcbal)
  
  dat.tmp <- dat.arg
  dat.tmp$foo <- NA
  dat.tmp$foo[dat.tmp$z==1] <- matched.ctrl
  dat.tmp$foo[dat.tmp$z==0][matched.ctrl] <- matched.ctrl
  
  return(dat.tmp$foo)    
}

## preprocesses the results of pair matching for an analysis
## using `senm'.
cast.senm <- function(dat, ms.arg, two.outcomes=FALSE) {
  ms <- as.vector(ms.arg)
  
  y <- dat$y[!is.na(ms)]
  mset <- ms[!is.na(ms)]
  z <- dat$z[!is.na(ms)]
  
  dico.names <- unique(mset)
  dico <- seq(length(dico.names))
  names(dico) <- dico.names
  mset <- as.integer(dico[mset])
  
  if(two.outcomes==FALSE) {
    return(list(y=y, mset=mset, z=z))
  } else {
    y2 <- dat$y2[!is.na(ms)]
    return(list(y=y, y2=y2, mset=mset, z=z))
  }
}


##############################
#Part III. Calculate ATE #####
##############################

### For Method 3 The naive approach: Delete all the missing data
# Generate Z column for the imputed data set
threshold = mean(oze$Base.tem, na.rm=TRUE)
oze_dropped$z <- rep(NA, nrow(oze_dropped))
oze_dropped[oze_dropped$Base.tem <= threshold, ][, "z"] <- 0
oze_dropped[oze_dropped$Base.tem > threshold, ][, "z"] <- 1
oze_dropped <- dplyr::select(oze_dropped, -c(Base.tem))

## Estimator 1: DIM 
gp0_oze_dropped = filter(oze_dropped, z==0)
gp1_oze_dropped  = filter(oze_dropped, z==1)
N0_oze_dropped = nrow(gp0_oze_dropped)
N1_oze_dropped = nrow(gp1_oze_dropped)
N_oze_dropped = nrow(oze_dropped)
tao10_oze_dropped = sum(gp1_oze_dropped$maxO3)/N1_oze_dropped -sum(gp0_oze_dropped$maxO3)/N0_oze_dropped 
V0_oze_dropped = var(gp0_oze_dropped$maxO3)
V1_oze_dropped = var(gp1_oze_dropped$maxO3)
V_hat_oze_dropped= V1_oze_dropped/N1_oze_dropped + V0_oze_dropped/N0_oze_dropped
CI_left_oze_dropped=tao10_oze_dropped-qnorm(0.975)*sqrt(V_hat_oze_dropped)
CI_right_oze_dropped=tao10_oze_dropped+qnorm(0.975)*sqrt(V_hat_oze_dropped)
CI_oze_dropped = c(CI_left_oze_dropped, CI_right_oze_dropped)



########################################################################################


### For Method 1 PCA 
# Generate Z column for the imputed data set
threshold = mean(oze$Base.tem, na.rm=TRUE)
imp.pca$z <- rep(NA, nrow(imp.pca))
imp.pca[imp.pca$Base.tem <= threshold, ][, "z"] <- 0
imp.pca[imp.pca$Base.tem > threshold, ][, "z"] <- 1
imp.pca <- select(imp.pca, select = -Base.tem)

## Estimator: DIM 
gp0_imp.pca = filter(imp.pca, z==0)
gp1_imp.pca  = filter(imp.pca, z==1)
N0_imp.pca = nrow(gp0_imp.pca)
N1_imp.pca = nrow(gp1_imp.pca)
N_imp.pca = nrow(imp.pca)
tao10_imp.pca = sum(gp1_imp.pca$maxO3)/N1_imp.pca -sum(gp0_imp.pca$maxO3)/N0_imp.pca 
V0_imp.pca = var(gp0_imp.pca$maxO3)
V1_imp.pca = var(gp1_imp.pca$maxO3)
V_hat_imp.pca= V1_imp.pca/N1_imp.pca + V0_imp.pca/N0_imp.pca
CI_left_imp.pca=tao10_imp.pca-qnorm(0.975)*sqrt(V_hat_imp.pca)
CI_right_imp.pca=tao10_imp.pca+qnorm(0.975)*sqrt(V_hat_imp.pca)
CI_imp.pca = c(CI_left_imp.pca, CI_right_imp.pca)




########################################################################################
### For Method 4 Mean
## Generate Z column for the imputed data set
add_z <- function(oze,df){
  threshold = mean(oze$Base.tem, na.rm=TRUE)
  df$z <- rep(NA, nrow(df))
  df$z <- ifelse(df$Base.tem <= threshold, 0, 1)
  df <- select(df, select = -Base.tem)
  return(df)
}

## Estimator: DIM 
for(i in 1:ncol(oze)) {
  oze[ , i][is.na(oze[ , i])] <- mean(oze[ , i], na.rm = TRUE)
}
head(oze) # Check first 6 rows after substitution by mean
oze = add_z(oze,oze)
Y1 = filter(oze, z == 1)
Y0 = filter(oze, z == 0)
ate = sum(Y1$maxO3)/nrow(Y1) - sum(Y0$maxO3)/nrow(Y0)
ate

V1 = var(Y1$maxO3)
V0 = var(Y0$maxO3)
V_hat = V1 / nrow(Y1) + V0 / nrow(Y0) 
V_hat






### For Method 7 Random Forest  
## Generate Z column for the imputed data set
threshold = mean(oze$Base.tem, na.rm=TRUE)
oze.tree$z <- rep(NA, nrow(oze.tree))
oze.tree[oze.tree$Base.tem <= threshold, ][, "z"] <- 0
oze.tree[oze.tree$Base.tem > threshold, ][, "z"] <- 1
oze.tree <- select(oze.tree, select = -Base.tem)

# Three methods for ATE estimators 
## Estimator 1: DIM 
gp0_oze.tree = filter(oze.tree, z==0)
gp1_oze.tree  = filter(oze.tree, z==1)
N0_oze.tree = nrow(gp0_oze.tree)
N1_oze.tree = nrow(gp1_oze.tree)
N_oze.tree = nrow(oze.tree)
tao10_oze.tree = sum(gp1_oze.tree$maxO3)/N1_oze.tree -sum(gp0_oze.tree$maxO3)/N0_oze.tree 
V0_oze.tree = var(gp0_oze.tree$maxO3)
V1_oze.tree = var(gp1_oze.tree$maxO3)
V_hat_oze.tree= V1_oze.tree/N1_oze.tree + V0_oze.tree/N0_oze.tree
CI_left_oze.tree=tao10_oze.tree-qnorm(0.975)*sqrt(V_hat_oze.tree)
CI_right_oze.tree=tao10_oze.tree+qnorm(0.975)*sqrt(V_hat_oze.tree)
CI_oze.tree = c(CI_left_oze.tree, CI_right_oze.tree)


########################################################################################

### For Method 2 Multiple Imputation of PCA  
# Compare density of original and imputed Base.tem
par(mfrow=c(2,2))
plot(density(na.omit(oze$Pressure.height)),col='blue', ann=FALSE)
lines(density(res.MIPCA[[1]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Pressure.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Pressure.height),col='red',ann=FALSE)
mtext("Pressure.height", side=3, line=0)


plot(density(na.omit(oze$Humidity)),col='blue', ann=FALSE)
lines(density(res.MIPCA[[1]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Humidity),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Humidity),col='red',ann=FALSE)
mtext("Humidity", side=3, line=0)

plot(density(na.omit(oze$Base.height)),col='blue', ann=FALSE)
lines(density(res.MIPCA[[1]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Base.height),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Base.height),col='red',ann=FALSE)
mtext("Base.height", side=3, line=0)



plot(density(na.omit(oze$Base.tem)),col='blue', ann=FALSE)
lines(density(res.MIPCA[[1]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[2]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[3]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[4]]$Base.tem),col='red',ann=FALSE)
lines(density(res.MIPCA[[5]]$Base.tem),col='red',ann=FALSE)
mtext("Base.tem", side=3, line=0)


# Generate Z column for the imputed data set
threshold = mean(oze$Base.tem, na.rm=TRUE)
for (i in 1:5){
  res.MIPCA[[i]]$z <- rep(NA, nrow(res.MIPCA[[i]]))
  res.MIPCA[[i]][res.MIPCA[[i]]$Base.tem <= threshold, ][, "z"] <- 0
  res.MIPCA[[i]][res.MIPCA[[i]]$Base.tem > threshold, ][, "z"] <- 1
  res.MIPCA[[i]] <- dplyr::select(res.MIPCA[[i]], -Base.tem)
}

## Estimator 1: DIM 
tao10_mipca = rep(0, 5)
V_hat_mipca = rep(0, 5)
CI_left_mipca = rep(0, 5)
CI_right_mipca = rep(0, 5)

for (i in 1:5){
  gp0 = filter(res.MIPCA[[i]], z==0)
  gp1  = filter(res.MIPCA[[i]], z==1)
  N0 = nrow(gp0)
  N1 = nrow(gp1)
  N = nrow(res.MIPCA[[i]])
  tao10 = sum(gp1$maxO3)/N1 -sum(gp0$maxO3)/N0
  tao10_mipca[i] = tao10
  
  V0 = var(gp0$maxO3)
  V1 = var(gp1$maxO3)
  V_hat = V1/N1 + V0/N0
  V_hat_mipca[i] = V_hat
  
  CI_left = tao10-qnorm(0.975)*sqrt(V_hat)
  CI_right = tao10+qnorm(0.975)*sqrt(V_hat)
  CI = c(CI_left, CI_right)
  CI_left_mipca[i] = CI_left
  CI_right_mipca[i] = CI_right 
}
#mean_estimator
mean_tao10_mipca = mean(tao10_mipca)
#within-imputation variance
wi_var_tao10_mipca = var(V_hat_mipca)
#between-imputation variance
bi_var_tao10_mipca = sum((tao10_mipca-mean_tao10_mipca)**2)/4
#total_variance
t_var_tao10_mipca = wi_var_tao10_mipca + (1+1/5)*bi_var_tao10_mipca
