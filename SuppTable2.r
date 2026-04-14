library (performance)


gtex <- read.delim("renormedGtexData.txt", h=T)

## Supplementary Table 2

summary(g0 <- glm(totMuts ~ scale(express)*sex*scale(phasedFrac) + 
                       scaleRepTime + scaleGC +  +
                       offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  gtex))
summary(g1 <- glm(totMuts ~ scale(express)*sex +
                      scale(express)*scale(phasedFrac) + 
                      sex*scale(phasedFrac)+
                      scaleRepTime + scaleGC + 
                      offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  gtex))
## compare the 3- and three 2-ways: drop 3-ways
anova(g1, g0, test="F")

## this is the FINAL model
summary(g2a <- glm(totMuts ~ scale(express)*sex +
                       sex*scale(phasedFrac) +
                      scale(phasedFrac) +
                       scaleRepTime + scaleGC + 
                       offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  gtex))
anova(g2a, g1, test="F") ## suggests this is better (2 of the 2-ways) 
check_collinearity(g2a)


## compare to model without sex
summary(nosex <- glm(totMuts ~ scale(express) +
                      scale(phasedFrac) +
                       scaleRepTime + scaleGC + 
                       offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  gtex))
anova(nosex, g2a, test="F")
