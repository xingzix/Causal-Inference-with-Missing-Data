# mean imputation 
mdf = oze

for(i in 1:ncol(mdf)) {
  mdf[ , i][is.na(mdf[ , i])] <- mean(mdf[ , i], na.rm = TRUE)
}
head(mdf) # Check first 6 rows after substitution by mean
mdf = add_z(oze,mdf)
Y1 = filter(mdf, z == 1)
Y0 = filter(mdf, z == 0)
ate = sum(Y1$maxO3)/nrow(Y1) - sum(Y0$maxO3)/nrow(Y0)
ate

V1 = var(Y1$maxO3)
V0 = var(Y0$maxO3)
V_hat = V1 / nrow(Y1) + V0 / nrow(Y0) 
V_hat

# linear regression
library(imputeTS)
linear_oze <- data.frame(lapply(oze, function(X) na_interpolation(X)))
lin = add_z(oze,linear_oze)
Y12 = filter(lin, z == 1)
Y02 = filter(lin, z == 0)
ate2 = sum(Y12$maxO3)/nrow(Y12) - sum(Y02$maxO3)/nrow(Y02)
ate2

V12 = var(Y12$maxO3)
V02 = var(Y02$maxO3)
V_hat2 = V12 / nrow(Y12) + V02 / nrow(Y02) 
V_hat2