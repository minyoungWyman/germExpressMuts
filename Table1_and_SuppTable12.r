
library (glmnet)
library(performance)

patstages <- read.delim("patMut_stages.txt", h=T)

###############################################################

## Table 1
summary(g1 <- glm(pat_totMuts ~ scalestage1 + scalestage2 +
                  scalestage3 +scalestage4 +
                  scalestage5+scalestage6 +
                  scalestage7+scalestage8 +
                  scalestage9+scalestage10+
                  scalestage11+scalestage12+
                  scalestage13 +
                  scaleRepTime + scaleGC + scalePatPhasedFrac +
                       offset(log(totalLen)),
                  family=poisson(link=log),
                  patstages))
check_collinearity(g1)

## keep only the germline stages and one somatic
summary(g2 <- glm(pat_totMuts ~ scalestage1 + scalestage2 +
                  scalestage3 + scalestage4 +
                  scalestage5 + scalestage6 +
                  scalestage7 + scalestage8 + scalestage12+
                  scaleRepTime + scaleGC + scalePatPhasedFrac +
                       offset(log(totalLen)),
                  family=poisson(link=log),
                  patstages))
check_collinearity(g2)

## final model
summary(g3 <- glm(pat_totMuts ~ scalestage1 + scalestage3 + 
                      scaleRepTime + scaleGC + scalePatPhasedFrac +
                  offset(log(totalLen)),
                  family=poisson(link=log),
                  patstages))
check_collinearity(g3)

########################################################################

## Supplementary Table 12

datMat <- data.matrix(patstages[,c("scalestage1", 
                                "scalestage3",
                                "scaleRepTime", "scaleGC",
                                "scalePatPhasedFrac")])
y <- patstages$pat_totMuts
logoff <- log(patstages$totalLen)

## LASSO
cvLasso1 <- cv.glmnet(x=datMat, y=y,
                      alpha=1,
                      offset=logoff,
                      family=poisson(link="log"),
                      standardize=FALSE)
bestLassoLam1 <- cvLasso1$lambda.min
bestLassoModel1 <- glmnet(datMat, y, alpha=1, lambda=bestLassoLam1,
                          offset=logoff, family=poisson("log"))
coef(bestLassoModel1)

## RIDGE
datMat <- data.matrix(patstages[,c("scalestage1",
                                   "scalestage3", 
                                   "scaleRepTime", "scaleGC",
                                   "scalePatPhasedFrac")])
y <- patstages$pat_totMuts
logoff <- log(patstages$totalLen)

ridge1 <- cv.glmnet(x=datMat, y=y,
                    offset=logoff,
                    alpha=0,
                    family=poisson(link="log"),
                    standardize=FALSE)
bestRidgeLam1 <- ridge1$lambda.min
bestRidgeModel1 <- glmnet(datMat, y, alpha=0, lambda=bestRidgeLam1,
                          offset=logoff, family=poisson("log"))
coef(bestRidgeModel1)
plot(ridge1)


## put it all into one dataframe
data.frame(ridge=coef(bestRidgeModel1)[,1],
           lasso=coef(bestLassoModel1)[,1],
           pois=data.frame(coef(g3))[,1])
