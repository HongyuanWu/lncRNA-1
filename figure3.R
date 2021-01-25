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
#============================================================
#for exon-based distribution:ks.test for both mRNA and lncRNA
#============================================================
load(file = "lncRNA_exon.Rdata") 
load(file = "mRNA_exon.Rdata")

p1=which(as.numeric(as.character(lrna_exon$exonnum))==1)
p2=which(as.numeric(as.character(lrna_exon$exonnum))>1)

lrna1=lrna_exon[p1,]
lrna2=lrna_exon[p2,]

p1=which(as.numeric(as.character(mrna_exon$exonnum))==1)
p2=which(as.numeric(as.character(mrna_exon$exonnum))>1)

mrna1=mrna_exon[p1,]
mrna2=mrna_exon[p2,]
 
lrna_ks = ks.test(
  as.numeric(as.character(lrna1$mean)),
  as.numeric(as.character(lrna2$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

mrna_ks = ks.test(
  as.numeric(as.character(mrna1$mean)),
  as.numeric(as.character(mrna2$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

lrna1_mrna1_ks = ks.test(
  as.numeric(as.character(lrna1$mean)),
  as.numeric(as.character(mrna1$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

lrna1_mrna2_ks = ks.test(
  as.numeric(as.character(lrna1$mean)),
  as.numeric(as.character(mrna2$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

lrna2_mrna1_ks = ks.test(
  as.numeric(as.character(lrna2$mean)),
  as.numeric(as.character(mrna1$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

lrna2_mrna2_ks = ks.test(
  as.numeric(as.character(lrna2$mean)),
  as.numeric(as.character(mrna2$mean)),
  alternative = c("two.sided", "less", "greater")[1],
)

lp1 = ecdf(as.numeric(as.character(lrna1$mean)))
lp2 = ecdf(as.numeric(as.character(lrna2$mean)))
mp1 = ecdf(as.numeric(as.character(mrna1$mean)))
mp2 = ecdf(as.numeric(as.character(mrna2$mean)))

par(mfrow=c(1,1))
b=0.05
f=1

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

plot(mp1, col = "black", lwd = 2, add=TRUE,do.points=FALSE)
plot(mp2, col = "green",  lwd = 2,add = TRUE)
plot(lp1, col = "red", add = TRUE, lwd = 2)
plot(lp2, col = "blue", add = TRUE, lwd = 2)

legend(
  x = 30,
  y = 20*b,
  legend = c("lnc-human1"),
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
  y = 19*b,
  legend = c("lnc-human2" ),
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
  y = 18*b,
  legend = c("m-human1"),
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
  y = 17*b,
  legend = c("m-human2"),
  col = c("green"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

text(15,12*b,labels = "Kolmogorov-Smirnov test",adj = 0,cex=f)
lines(c(15,50),c(11*b,11*b),lwd=1)
text(15,10*b,labels = "P1",adj = 0,cex=f)
text(25,10*b,labels = "P2",adj = 0,cex=f)
text(35,10*b,labels = "D",adj = 0,cex=f)
text(42,10*b,labels = "p(same)",adj = 0,cex=f)
lines(c(15,50),c(9*b,9*b),lwd=1)
text(15,8*b,labels = "lnc-human1",adj = 0,cex=f)
text(25,8*b,labels = "lnc-human2",adj = 0,cex=f)
text(35,8*b,labels = format(lrna_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,8*b,labels = format(lrna_ks$p.value, digits = 3 ),adj = 0,cex=f)

text(15,7*b,labels = "lnc-human1",adj = 0,cex=f)
text(25,7*b,labels = "m-human1",adj = 0,cex=f)
text(35,7*b,labels = format(lrna1_mrna1_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,7*b,labels = format(lrna1_mrna1_ks$p.value, digits = 3 ),adj = 0,cex=f)

text(15,6*b,labels = "lnc-human1",adj = 0,cex=f)
text(25,6*b,labels = "m-human2",adj = 0,cex=f)
text(35,6*b,labels = format(lrna1_mrna2_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,6*b,labels = format(lrna1_mrna2_ks$p.value, digits = 3 ),adj = 0,cex=f)

text(15,5*b,labels = "lnc-human2",adj = 0,cex=f)
text(25,5*b,labels = "m-human1",adj = 0,cex=f)
text(35,5*b,labels = format(lrna2_mrna1_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,5*b,labels = format(lrna2_mrna1_ks$p.value, digits = 3 ),adj = 0,cex=f)

text(15,4*b,labels = "lnc-human2",adj = 0,cex=f)
text(25,4*b,labels = "m-human2",adj = 0,cex=f)
text(35,4*b,labels = format(lrna2_mrna2_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,4*b,labels = format(lrna2_mrna2_ks$p.value, digits = 3 ),adj = 0,cex=f)

text(15,3*b,labels = "m-human1",adj = 0,cex=f)
text(25,3*b,labels = "m-human2",adj = 0,cex=f)
text(35,3*b,labels = format(mrna_ks$statistic, digits = 3 ),adj = 0,cex=f)
text(42,3*b,labels = format(mrna_ks$p.value, digits = 3 ),adj = 0,cex=f)

lines(c(15,50),c(2*b,2*b),lwd=1)
