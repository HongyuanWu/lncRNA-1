library(stringi)
library(seqinr)
library(RMySQL)
library(keras)
library(tensorflow)
library(pROC)
library(ggplot2)
library(corrplot)
library(PKNCA)
library(visNetwork)
#=================================================
#lncRNA classification and subgroup analysis
load(file = "lncRNA_class.Rdata")

xx = which(as.character(lrna_class$class) == "Divergent")
xrna = lrna_class[-xx, ]

uclass = as.character(unique(xrna$class))
for (x in uclass)
{
  xx = which(as.character(xrna$class) == x)
  assign(x, xrna$mean[xx])
}

sense_inter_ks = ks.test(
  as.numeric(as.character(Sense)),
  as.numeric(as.character(Intergenic)),
  alternative = c("two.sided", "less", "greater")[1],
)

sense_intro_ks = ks.test(
  as.numeric(as.character(Sense)),
  as.numeric(as.character(Intronic)),
  alternative = c("two.sided", "less", "greater")[1],
)

sense_anti_ks = ks.test(
  as.numeric(as.character(Sense)),
  as.numeric(as.character(Antisense)),
  alternative = c("two.sided", "less", "greater")[1],
)

inter_intro_ks = ks.test(
  as.numeric(as.character(Intergenic)),
  as.numeric(as.character(Intronic)),
  alternative = c("two.sided", "less", "greater")[1],
)

inter_anti_ks = ks.test(
  as.numeric(as.character(Intergenic)),
  as.numeric(as.character(Antisense)),
  alternative = c("two.sided", "less", "greater")[1],
)

intro_anti_ks = ks.test(
  as.numeric(as.character(Intronic)),
  as.numeric(as.character(Antisense)),
  alternative = c("two.sided", "less", "greater")[1],
)

sen = ecdf(as.numeric(as.character(Sense)))
inter = ecdf(as.numeric(as.character(Intergenic)))
intro = ecdf(as.numeric(as.character(Intronic)))
anti = ecdf(as.numeric(as.character(Antisense)))

par(mfrow = c(1, 1))
b = 0.05
f = 1.2

plot(
  c(0),
  c(0),
  main = "",
  cex = f,
  col = "white",
  font = 1,
  xaxt = "n",
  xlim = c(0, 50),
  xlab = "Half-life (h)",
  yaxt = "n",
  ylim = c(0, 1),
  ylab = "Percentage of transcripts (%)",
  asp = c(grid(nx = 10, ny = 10))
)

axis(
  2,
  at = seq(0, 1, 0.1),
  las = 1,
  labels = seq(0, 1, 0.1) * 100,
  font.axis = 1,
  cex.axis = f
)

axis(
  1,
  at = seq(0, 50, 5),
  las = 1,
  labels = seq(0, 50, 5),
  font.axis = 1,
  cex.axis = f
)

plot(
  intro,
  col = "red",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)
plot(
  anti,
  col = "blue",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)
plot(
  sen,
  col = "black",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)
plot(
  inter,
  col = "green",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)

legend(
  x = 30,
  y = 20 * b,
  legend = c("Intronic"),
  col = c("red"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = 30,
  y = 19 * b,
  legend = c("Antisense"),
  col = c("blue"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = 30,
  y = 18 * b,
  legend = c("Sense"),
  col = c("black"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)
legend(
  x = 30,
  y = 17 * b,
  legend = c("Intergenic"),
  col = c("green"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

text(15,
     12 * b,
     labels = "Kolmogorov-Smirnov test",
     adj = 0,
     cex = f)
lines(c(15, 50), c(11 * b, 11 * b), lwd = 1)
text(15,
     10 * b,
     labels = "P1",
     adj = 0,
     cex = f)
text(25,
     10 * b,
     labels = "P2",
     adj = 0,
     cex = f)
text(35,
     10 * b,
     labels = "D",
     adj = 0,
     cex = f)
text(42,
     10 * b,
     labels = "p(same)",
     adj = 0,
     cex = f)
lines(c(15, 50), c(9 * b, 9 * b), lwd = 1)
text(15,
     8 * b,
     labels = "Sense",
     adj = 0,
     cex = f)
text(25,
     8 * b,
     labels = "Intergenic",
     adj = 0,
     cex = f)
text(
  35,
  8 * b,
  labels = format(sense_inter_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  8 * b,
  labels = format(sense_inter_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     7 * b,
     labels = "Sense",
     adj = 0,
     cex = f)
text(25,
     7 * b,
     labels = "Intronic",
     adj = 0,
     cex = f)
text(
  35,
  7 * b,
  labels = format(sense_intro_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  7 * b,
  labels = format(sense_intro_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     6 * b,
     labels = "Sense",
     adj = 0,
     cex = f)
text(25,
     6 * b,
     labels = "Antisense",
     adj = 0,
     cex = f)
text(
  35,
  6 * b,
  labels = format(sense_anti_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  6 * b,
  labels = format(sense_anti_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     5 * b,
     labels = "Intergenic",
     adj = 0,
     cex = f)
text(25,
     5 * b,
     labels = "Intronic",
     adj = 0,
     cex = f)
text(
  35,
  5 * b,
  labels = format(inter_intro_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  5 * b,
  labels = format(inter_intro_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     4 * b,
     labels = "Intergenic",
     adj = 0,
     cex = f)
text(25,
     4 * b,
     labels = "Antisense",
     adj = 0,
     cex = f)
text(
  35,
  4 * b,
  labels = format(inter_anti_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  4 * b,
  labels = format(inter_anti_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     3 * b,
     labels = "Intronic",
     adj = 0,
     cex = f)
text(25,
     3 * b,
     labels = "Antisense",
     adj = 0,
     cex = f)
text(
  35,
  3 * b,
  labels = format(intro_anti_ks$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  3 * b,
  labels = format(intro_anti_ks$p.value, digits = 3),
  adj = 0,
  cex = f
)

lines(c(15, 50), c(2 * b, 2 * b), lwd = 1)
