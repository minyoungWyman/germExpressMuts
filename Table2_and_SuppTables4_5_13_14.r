library(MASS)
library(performance)

matstages <- read.delim("matMut_stages.txt", h=T)

########################################################################

## Table 2
summary(t2 <- glm.nb(mat_totMuts ~ scale(ArithMean_O) +
                          + scaleRepTime +
                         scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                     matstages))
check_collinearity(t2)

########################################################################

## Supplementary Table 4
summary(priO <- glm.nb(mat_totMuts ~ scale(ArithMean_primorO) +
                         scaleRepTime + scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                       matstages))
check_collinearity(priO)

########################################################################

## Supplementary Table 5
summary(st5 <- glm(all_totMuts ~ scale(ArithMean_O) +
                         scaleRepTime + scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                    family=quasipoisson(link=log), matstages))

check_collinearity(st5)

summary(st5 <- glm(all_totMuts ~ scale(ArithMean_primorO) +
                         scaleRepTime + scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                    family=quasipoisson(link=log), matstages))

check_collinearity(st5)

########################################################################

## Supplementary Table 13
summary(og <- glm.nb(mat_totMuts ~ scale(ArithMean_O) +
                         scale(ArithMean_G) + scaleRepTime +
                         scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                      matstages))
check_collinearity(og)

###############################################################

## Supplementary Table 14
summary(st4 <- glm.nb(mat_totMuts ~ scale(ArithMean_G) +
                          + scaleRepTime +
                         scaleGC + scale(phasedFracMatBinVal) +
                         offset(log(totalLen)),
                       matstages))
check_collinearity(st4)

########################################################################

## correlation between average of primordial oocytes and average of all oocyte types

cor.test(matstages$ArithMean_O,
         matstages$ArithMean_primorO,
         method="sp")
