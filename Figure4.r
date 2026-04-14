library(ggplot2)
library(dplyr)

nts_levels <- c("A>G", "C>G", "A>C", "A>T", "C>A", "nonCpG C>T", "CpG>TpG")
df_long <- read.delim("tasym_tests.txt", h = TRUE) %>%
  mutate(
    NTS = factor(NTS, levels = nts_levels),
    Sex = factor(Sex, levels = c("log_male", "log_fem"))
  )

pdf("Figure4.pdf", h = 8, w = 8)
ggplot(df_long, aes(x = NTS, y = log_value, color = Sex)) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "grey70", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position  = position_dodge(width = 0.2),
    width     = 0.15,
    linewidth = 0.5
  ) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_text(
    data        = df_long %>% distinct(NTS, .keep_all = TRUE),
    aes(x = NTS, y = y_pos, label = p_label),
    color       = "black",
    inherit.aes = FALSE,
    size        = 4,
    vjust       = 0,
    parse       = TRUE,
    nudge_x     = 0.1
  ) +
  scale_color_manual(
    values = c("log_male" = "#2166AC", "log_fem" = "#D6604D"),
    labels = c("log_male" = "Male",    "log_fem" = "Female")
  ) +
  scale_y_continuous(breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "", y = "T-asymmetry", color = "Sex") +
  theme_classic() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, margin = margin(t = 5),
                                      size = 16, color = "black"),
    plot.title         = element_text(hjust = 0, size = 20),
    axis.title.x       = element_text(margin = margin(t = 20), size = 18, vjust = 2),
    axis.text.y        = element_text(margin = margin(r = 5),  size = 16, color = "black"),
    axis.title.y       = element_text(margin = margin(r = 25), size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position    = c(0.92, 0.1),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.text        = element_text(size = 13),
    legend.title       = element_text(size = 14)
  )
dev.off()
