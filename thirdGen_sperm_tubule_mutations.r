
library(performance)

## consistency checks

##################################################################

## using only third generation phased mutations from Jonsson et al. 2017

## maternal
mat <- read.delim("mat_thirdGenMuts.txt", h=T)

## check with adult oocytes
summary(m1 <- glm(matThirdGenMuts ~ ArithMean_O +             
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=poisson(link=log), mat))
check_collinearity(m1)

## check with fetal germline expression
summary(m2 <- glm(matThirdGenMuts ~ atlasCPM_geoMean +
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=poisson(link=log), mat))
check_collinearity(m2)

## --------------------------------------------------------------

## paternal
pat <- read.delim("pat_thirdGenMuts.txt", h=T)

## check with adult germline expression
summary(p1 <- glm(patThirdGenMuts ~ scale(cpm_geoMean) +             
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=poisson(link=log), pat))
check_collinearity(p1)


## check with fetal germline expression
summary(p2 <- glm(patThirdGenMuts ~ scale(fetalmale_CPM_log2geoMean) +
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=poisson(link=log), pat))
check_collinearity(p2)

##################################################################

## using male-only mutations derived from
## sperm and seminiferous tubules

maleonly <- read.delim("spermSemTubuleMuts.txt", h=T)

## check with fetal germline expression
summary(sp1 <- glm(sperm_totMuts ~ scale(fetalGermExpress) +
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=quasipoisson(link=log),
                  maleonly))
check_collinearity(sp1)

## check with adult germline expression
summary(sp2 <- glm(sperm_totMuts ~ scale(adultGermExpress) +
                      scaleRepTime + scaleGC +
                      offset(log(totalLen)),
                  family=quasipoisson(link=log),
                  maleonly))
check_collinearity(sp2)
