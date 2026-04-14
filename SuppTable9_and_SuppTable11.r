library(performance)
library(glmnet)

rca <- read.delim("reproCellAtlas_CpGs_and_allOtherMutations.txt",
                   h=T)

###############################################################

## Supplementary Table 9
summary(sbs5 <- glm(mutsAllOthers~ scale(express)*sex*scale(phasedFrac) +
                scaleRepTime + scaleGC +
                offset(log(DenomAllOthers)),
            family=quasipoisson(link=log),
            rca))
check_collinearity(sbs5)

## high collinearity detected
## test for coefficient stability

# pre-scale predictors
rca$express_s    <- as.numeric(scale(rca$express))
rca$phasedFrac_s <- as.numeric(scale(rca$phasedFrac))

X <- model.matrix(~ express_s * sex * phasedFrac_s +
                      scaleRepTime + scaleGC,
                  data = rca)[, -1]

y   <- rca$mutsAllOthers
off <- log(rca$DenomAllOthers)

# ridge (alpha = 0)
cv_ridge <- cv.glmnet(X, y,
                      family  = "poisson",
                      offset  = off,
                      alpha   = 0,
                      standardize = FALSE)

# lasso (alpha = 1)
cv_lasso <- cv.glmnet(X, y,
                      family  = "poisson",
                      offset  = off,
                      alpha   = 1,
                      standardize = FALSE)

# coefficients at optimal lambda
coef(cv_ridge, s = "lambda.min")
coef(cv_lasso, s = "lambda.min")

# compare to unregularized glm
data.frame(ridge=coef(cv_ridge, s = "lambda.min")[,1],
           lasso=coef(cv_lasso, s = "lambda.min")[,1],
           quasipois = data.frame(coef(sbs5))[,1])


###############################################################

## Supplementary Table 11
summary(sbs1 <- glm(cpgMuts~ scale(express) + sex +
                      scale(avgMeth) +
                      scale(phasedFrac) +
                      scaleRepTime + scaleGC +
                      offset(log(cpgDenomBP)),
                  family=poisson(link=log),
                  rca))
check_collinearity(sbs1)

