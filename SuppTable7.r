library(performance)
library(glmnet)


###############################################################
## Supplementary Table 7

rca <- read.delim("pseudobulk_RCA_SuppTable7.txt", h=T)

summary(g0 <- glm(all_totMuts ~ scale(pat_express) + scale(mat_express) +
                      scale(phasedFracPatBinVal_patbins) +
                      scale(phasedFracMatBinVal_matbins) + 
                       scaleRepTime + scaleGC +  +
                       offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  rca))
cor.test(gonads$pat_express, gonads$mat_express)
check_collinearity(g0)


## -------------------------------------------------------
## check for coefficient stability with ridge
## ridge coefficient estimates similar

X <- model.matrix(~ scale(pat_express) + scale(mat_express) +
                      scale(phasedFracPatBinVal_patbins) +
                      scale(phasedFracMatBinVal_matbins) + 
                      scaleRepTime + scaleGC,
                  data = rca)[, -1]  # drop intercept
y <- rca$all_totMuts
# offset
off <- log(rca$denomBP)

# ridge (alpha = 0)
g0_ridge <- glmnet(X, y, 
                   family = "poisson",  # no quasipoisson in glmnet
                   offset = off,
                   alpha = 0, standardize=FALSE)

# cross-validate to find optimal lambda
cv_ridge <- cv.glmnet(X, y, family = "poisson", offset = off, alpha = 0)


## put it all into one dataframe
data.frame(ridge=coef(cv_ridge, s="lambda.min")[,1],
           pois=data.frame(coef(g0))[,1])
