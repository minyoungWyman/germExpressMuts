
library(performance)

adultTestisGerm <- read.delim("pseudobulk_testisGermline.txt", h=T)

###############################################################

## Supplementary Table 3
summary(malephased <- glm(pat_totMuts ~ scale(pseudobulkTestisGerm)+
                  scale(meanTime) + scale(GCperc) + scale(phasedFracPatBinVal) +
                       offset(log(totalLen)),
                  family=poisson(link=log),
                  adultTestisGerm))
check_collinearity(malephased)

###############################################################

## Supplementary Table 6
summary(allMutations <- glm(all_totMuts ~ scale(pseudobulkTestisGerm) +
                      scale(meanTime) + scale(GCperc) +
                      scale(phasedFracPatBinVal) +
                      offset(log(totalLen)),
                  family=quasipoisson(link=log),
                  adultTestisGerm))
check_collinearity(allMutations)
