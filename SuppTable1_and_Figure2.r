library(performance)
library(ggplot2)
library(dplyr)


rca <- read.delim("pseudobulk_reproCellAtlasData.txt", h=T)

###############################################################

## Supplementary Table 1

## the full model:
summary(g1 <- glm(muts ~ scale(logexpress)*sex +
                      sex*scale(phasedFrac)  +
                      scale(phasedFrac)*scale(logexpress) +
                      scale(meanTime) + scale(GCperc) +
                      offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  rca))

## keep 2 of the two-way interactions
summary(g2 <- glm(muts ~ scale(logexpress)*sex +
                      sex*scale(phasedFrac)  +
                      scale(meanTime) + scale(GCperc) +
                      offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  rca))
anova(g2, g1, test="F")
check_collinearity(g2)

## compare to a model without sex
summary(nosex <- glm(muts ~ scale(logexpress) +
                      scale(phasedFrac)  +
                      scale(meanTime) + scale(GCperc) +
                      offset(log(denomBP)),
                  family=quasipoisson(link=log),
                  rca))
anova(g2, nosex, test="F")


###############################################################

## Figure 2

n_trios <- 7596

## ----------------------------------------------------------
## 1. Scale covariates and compute covariate-corrected counts 
b <- coef(g2)

rca_corr <- rca %>%
  mutate(
    sc_logexpress = as.numeric(scale(logexpress)),
    sc_phasedFrac = as.numeric(scale(phasedFrac)),
    sc_meanTime   = as.numeric(scale(meanTime)),
    sc_GCperc     = as.numeric(scale(GCperc)),
    sex_mother    = as.numeric(sex == "mother")
  ) %>%
  mutate(
    nuisance = b["scale(meanTime)"]             * sc_meanTime   +
               b["scale(GCperc)"]               * sc_GCperc     +
               b["scale(phasedFrac)"]           * sc_phasedFrac +
               b["sexmother:scale(phasedFrac)"] * sex_mother * sc_phasedFrac,
    muts_corrected = muts * exp(-nuisance)
  )

## ----------------------------------------------------------
## 2. Individual gene-level corrected rates
zero_jitter <- c(father = 1e-10, mother = 1.1e-10)

## RAW
indiv_points <- rca_corr %>%
   group_by(sex) %>%
   mutate(expr_bin = ntile(logexpress, 4)) %>%
   ungroup() %>%
   mutate(
     raw_rate   = (muts / denomBP) / n_trios,        ## <-- muts, not muts_corrected
     indiv_rate = ifelse(raw_rate <= 0, zero_jitter[sex], raw_rate)
   )

## ----------------------------------------------------------
## 3. Binned summary
plot_corr <- rca_corr %>%
  group_by(sex) %>%
  mutate(expr_bin = ntile(logexpress, 4)) %>%
  group_by(sex, expr_bin) %>%
  summarise(
    corrRate   = (sum(muts_corrected, na.rm = TRUE) / sum(denomBP, na.rm = TRUE)) / n_trios,
    totCorMuts = sum(muts_corrected,  na.rm = TRUE),
    totDenom   = sum(denomBP,         na.rm = TRUE),
    ## use within-sex scaled expression for x position
    meanExpr   = mean(logexpress, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se  = sqrt(corrRate / (totDenom * n_trios)),
    lwr = pmax(corrRate - 1.96 * se, 1e-10),
    upr = corrRate + 1.96 * se
  )

## ----------------------------------------------------------
## 4. Plot
sex_colors <- c(father = "#2166AC", mother = "#D6604D")
sex_labels <- c(father = "Male", mother = "Female")

p <- ggplot() +
  ## ghosted individual points (using within-sex scaled x)
  geom_point(
    data = indiv_points,
    aes(x = logexpress, y = indiv_rate, color = sex),
    alpha = 0.2, size = 0.7, shape = 16,
    show.legend = FALSE
  ) +
  geom_line(
    data = plot_corr,
    aes(x = meanExpr, y = corrRate, color = sex, group = sex),
    linewidth = 0.9
  ) +
  geom_errorbar(
    data = plot_corr,
    aes(x = meanExpr, ymin = lwr, ymax = upr, color = sex),
    width = 0.08, linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = plot_corr,
    aes(x = meanExpr, y = corrRate, color = sex, fill = sex),
    size = 3.5, shape = 21, stroke = 1.5
  ) +
  scale_color_manual(
    name   = "Sex",
    values = sex_colors,
    labels = sex_labels
  ) +
  scale_fill_manual(
    name   = "Parent of origin",
    values = c(father = "white", mother = "white"),
    labels = sex_labels
  ) +
  guides(
    color = guide_legend(override.aes = list(
      shape  = 21,
      size   = 3.5,
      stroke = 1.5,
      fill   = "white",
      alpha  = 1
    )),
    fill = "none"
  ) +
  scale_x_continuous(
      name = "Expression",
      limits = c(-1,NA)
  ) +
 scale_y_log10(
    name   = "Mutation rate",
    labels = scales::trans_format("log10", scales::math_format()),
    limits = c(1e-10, NA)
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.85, 0.2),
    legend.background      = element_rect(fill = "white", color = "white"),
    legend.title           = element_text(size = 13),
    legend.text            = element_text(size = 12),
    axis.title             = element_text(size = 16),
    axis.text              = element_text(size = 14)
  ) 
 
## ----------------------------------------------------------
ggsave("Figure2.pdf",
       plot   = p,        
       width  = 7,
       height = 7,
       dpi    = 300,
       bg     = "white")

ggsave("Figure2_adjusted.png",
       plot=p,
       width=8,
       height=8,
       dpi=300,
       bg="white")
