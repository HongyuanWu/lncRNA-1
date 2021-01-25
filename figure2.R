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
#====================================================
#ks.test for both lncRNAs and mRNAs
load(file = "lncRNA.Rdata")  #obtain lrna
load(file = "mRNA.Rdata")    #obtain mRNA

xx = ks.test(as.numeric(as.character(lrna$mean)),
             as.numeric(as.character(mrna$mean)),
             alternative = "two.sided")

mp = ecdf(as.numeric(as.character(mrna$mean)))
lp = ecdf(as.numeric(as.character(lrna$mean)))

par(mfrow = c(1, 1))

plot(
  c(0),
  c(0),
  main = "",
  cex = 1.2,
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
  cex.axis = 1
)

axis(
  1,
  at = seq(0, 50, 5),
  las = 1,
  labels = seq(0, 50, 5),
  font.axis = 1,
  cex.axis = 1
)

plot(mp, col = "red",  lwd = 2, add = TRUE)
plot(lp, col = "blue", lwd = 2, add = TRUE)

# Add a legend
legend(
  x = 20,
  y = 0.80,
  legend = c("lncRNA", "mRNA"),
  col = c("blue", "red"),
  lty = 1:1,
  cex = 0.8,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

text(20,
     0.40,
     labels = "Kolmogorov-Smirnov test",
     adj = 0,
     cex = 1.2)
lines(c(20, 45), c(0.35, 0.35), lwd = 1)
text(
  20,
  0.30,
  cex = 1.2,
  adj = 0,
  labels = stri_c("D value= ",
                  format(xx$statistic, digits = 3))
)

text(
  20,
  0.20,
  cex = 1.2,
  adj = 0,
  labels = stri_c("p(same)=",
                  format(xx$p.value, digits = 3))
)
