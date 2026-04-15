library(ggplot2)

quantplot <- read.delim("fetalGermExpressQuantiles_and_mutationRates.txt", h=T)


patQuant <- aggregate(quantplot$patMutRate,
                      by=list(pat_quant = quantplot$pat_quant),
                      mean)
matQuant <- aggregate(quantplot$matMutRate,
                      by=list(mat_quant = quantplot$mat_quant),
                      mean)

ratioMatrix <- outer(patQuant$x, matQuant$x, "/")
rownames(ratioMatrix) <- patQuant$pat_quant
colnames(ratioMatrix) <- matQuant$mat_quant


rng <- range(ratioMatrix)

df <- as.data.frame(as.table(ratioMatrix))
colnames(df) <- c("pat_quant", "mat_quant", "ratio")
df$pat_quant <- factor(df$pat_quant, levels = patQuant$pat_quant)
df$mat_quant <- factor(df$mat_quant, levels = matQuant$mat_quant)

midpoint <- mean(range(df$ratio))
df$text_color <- ifelse(df$ratio > midpoint, "white", "#08306B")

## png("Figure3_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pdf("Figure3.pdf", width = 8, height = 8)

ggplot(df, aes(mat_quant, pat_quant, fill = ratio)) +
  geom_tile(color = NA) +
  geom_text(aes(label = sprintf("%.2f", ratio), color = text_color), size = 5) +
  scale_color_identity() +                        
  scale_fill_gradient(
  high   = "#08306B",
  low    = "#F7FBFF",
  limits = range(df$ratio),
  breaks = seq(rng[1], rng[2], by = 0.25),
  labels = function(x) sprintf("%.1f", x)) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 15, color = "black"),
    axis.text.y  = element_text(size = 15, color = "black"),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10))) +
  labs(x = "female expression quantile", y = "male expression quantile")

dev.off()
