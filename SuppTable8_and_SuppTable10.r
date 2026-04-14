library(performance)

gtex <- read.delim("renormedGtexData_CpGs_and_allOtherMutations.txt",
                   h=T)

###############################################################

## Supplementary Table 8
summary(sbs5 <- glm(mutsAllOthers~ scale(express)*sex*scale(phasedFrac) +
                scaleRepTime + scaleGC +
                offset(log(DenomAllOthers)),
            family=quasipoisson(link=log),
            gtex))
check_collinearity(sbs5)

###############################################################

## Supplementary Table 10
summary(sbs1 <- glm(cpgMuts~ scale(express) + sex +
                      scale(avgMeth) +
                      scale(phasedFrac) +
                      scaleRepTime + scaleGC +
                      offset(log(cpgDenomBP)),
                  family=poisson(link=log),
                  gtex))
check_collinearity(sbs1)

gtex$sex <- relevel(factor(gtex$sex), ref="mother")
