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
#===============================================================
# nuclear or  cytoplasmic-based half-life comparison
#for lncRNAs
load(file = "lncRNA_nuclear.Rdata") 
load(file = "lncRNA_cytoplasm.Rdata")

cname = intersect(lrna_n$noncode, lrna_c$noncode)
lrna_n_com = lrna_n[match(cname, lrna_n$noncode), ]
lrna_c_com = lrna_c[match(cname, lrna_c$noncode), ]

lrna_n = lrna_n[match(setdiff(lrna_n$noncode, cname), lrna_n$noncode), ]
lrna_c = lrna_c[match(setdiff(lrna_c$noncode, cname), lrna_c$noncode), ]

#for mRNAs
load(file = "mRNA_nuclear.Rdata") 
load(file = "mRNA_cytoplasm.Rdata")
 
cname = intersect(mrna_n$transcript_id, mrna_c$transcript_id)
mrna_n_com = mrna_n[match(cname, mrna_n$transcript_id), ]
mrna_c_com = mrna_c[match(cname, mrna_c$transcript_id), ]

mrna_n = mrna_n[match(setdiff(mrna_n$transcript_id, cname),
                      mrna_n$transcript_id), ]
mrna_c = mrna_c[match(setdiff(mrna_c$transcript_id, cname),
                      mrna_c$transcript_id), ]

lrna_n_d = ecdf(as.numeric(as.character(lrna_n$half.life)))
lrna_c_d = ecdf(as.numeric(as.character(lrna_c$half.life)))
mrna_n_d = ecdf(as.numeric(as.character(mrna_n$half.life)))
mrna_c_d = ecdf(as.numeric(as.character(mrna_c$half.life)))

lrna_n_d_com = ecdf(as.numeric(as.character(lrna_n_com$half.life)))
lrna_c_d_com = ecdf(as.numeric(as.character(lrna_c_com$half.life)))
mrna_n_d_com = ecdf(as.numeric(as.character(mrna_n_com$half.life)))
mrna_c_d_com = ecdf(as.numeric(as.character(mrna_c_com$half.life)))


lrna_nc = ks.test(as.numeric(as.character(lrna_n$half.life)),
                  as.numeric(as.character(lrna_c$half.life)),
                  alternative = "two.sided")

mrna_nc = ks.test(as.numeric(as.character(mrna_n$half.life)),
                  as.numeric(as.character(mrna_c$half.life)),
                  alternative = "two.sided")

lrna_nc_com = ks.test(as.numeric(as.character(lrna_n_com$half.life)),
                      as.numeric(as.character(lrna_c_com$half.life)),
                      alternative = "two.sided")

mrna_nc_com = ks.test(as.numeric(as.character(mrna_n_com$half.life)),
                      as.numeric(as.character(mrna_c_com$half.life)),
                      alternative = "two.sided")


lrna_mrna_n = ks.test(as.numeric(as.character(lrna_n$half.life)),
                      as.numeric(as.character(mrna_n$half.life)),
                      alternative = "two.sided")

lrna_mrna_c = ks.test(as.numeric(as.character(lrna_c$half.life)),
                      as.numeric(as.character(mrna_c$half.life)),
                      alternative = "two.sided")

lrna_mrna_n_com = ks.test(as.numeric(as.character(lrna_n_com$half.life)),
                          as.numeric(as.character(mrna_n_com$half.life)),
                          alternative = "two.sided")

lrna_mrna_c_com = ks.test(as.numeric(as.character(lrna_c_com$half.life)),
                          as.numeric(as.character(mrna_c_com$half.life)),
                          alternative = "two.sided")


par(mfrow = c(1, 2))
plot(
  0,
  0,
  cex.axis = 1,
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
  lwd = 1,
  labels = seq(0, 1, 0.1) * 100,
  font.axis = 1,
  cex.axis = 1
)

axis(
  1,
  at = seq(0, 50, 5),
  las = 3,
  labels = seq(0, 50, 5),
  font.axis = 1,
  cex.axis = 1
)

title(main = "(A)")

b = 0.05
f = 1
plot(
  lrna_n_d,
  col = "red",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)
plot(
  lrna_c_d,
  col = "blue",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)
plot(
  mrna_n_d,
  col = "black",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)
plot(
  mrna_c_d,
  col = "green",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)

p1 = "lnc_nuc."
p2 = "lnc_cyt."
p3 = "m_nuc."
p4 = "m_cyt."

ystart = 0.95

#try to use legend
lx = 15
legend(
  x = lx,
  y = 0.9,
  legend = c("lnc-nuc."),
  col = c("red"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = lx,
  y = 0.85,
  legend = c("lnc-cyt."),
  col = c("blue"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = lx,
  y = 0.8,
  legend = c("m-nuc."),
  col = c("black"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)
legend(
  x = lx,
  y = 0.75,
  legend = c("m-cyt."),
  col = c("green"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)
#try end

text(
  1,
  ystart,
  labels = stri_c("lnc-nuc.=", nrow(lrna_n)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - b,
  labels = stri_c("lnc-cyt.=", nrow(lrna_c)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - 2 * b,
  labels = stri_c("m-nuc.=", nrow(mrna_n)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - 3 * b,
  labels = stri_c("m-cyt.=", nrow(mrna_c)),
  adj = 0,
  cex = f
)

 
ydown = 0.1
text(15,
     12 * b - ydown,
     labels = "Kolmogorov-Smirnov test",
     adj = 0,
     cex = f)
lines(c(15, 50), c(11 * b - ydown, 11 * b - ydown), lwd = 1)
text(15,
     10 * b - ydown,
     labels = "P1",
     adj = 0,
     cex = f)
text(25,
     10 * b - ydown,
     labels = "P2",
     adj = 0,
     cex = f)
text(35,
     10 * b - ydown,
     labels = "D",
     adj = 0,
     cex = f)
text(42,
     10 * b - ydown,
     labels = "p(same)",
     adj = 0,
     cex = f)
lines(c(15, 50), c(9 * b - ydown, 9 * b - ydown), lwd = 1)
text(15,
     8 * b - ydown,
     labels = p1,
     adj = 0,
     cex = f)
text(25,
     8 * b - ydown,
     labels = p2,
     adj = 0,
     cex = f)
text(
  35,
  8 * b - ydown,
  labels = format(lrna_nc$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  8 * b - ydown,
  labels = format(lrna_nc$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     7 * b - ydown,
     labels = p3,
     adj = 0,
     cex = f)
text(25,
     7 * b - ydown,
     labels = p4,
     adj = 0,
     cex = f)
text(
  35,
  7 * b - ydown,
  labels = format(mrna_nc$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  7 * b - ydown,
  labels = format(mrna_nc$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     6 * b - ydown,
     labels = p1,
     adj = 0,
     cex = f)
text(25,
     6 * b - ydown,
     labels = p3,
     adj = 0,
     cex = f)
text(
  35,
  6 * b - ydown,
  labels = format(lrna_mrna_n$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  6 * b - ydown,
  labels = format(lrna_mrna_n$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     5 * b - ydown,
     labels = p2,
     adj = 0,
     cex = f)
text(25,
     5 * b - ydown,
     labels = p4,
     adj = 0,
     cex = f)
text(
  35,
  5 * b - ydown,
  labels = format(lrna_mrna_c$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  5 * b - ydown,
  labels = format(lrna_mrna_c$p.value, digits = 3),
  adj = 0,
  cex = f
)
lines(c(15, 50), c(4 * b - ydown, 4 * b - ydown), lwd = 1)
#end for the first plot

plot(
  0,
  0,
  cex.axis = 1,
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
  lwd = 1,
  labels = seq(0, 1, 0.1) * 100,
  font.axis = 1,
  cex.axis = 1
)

axis(
  1,
  at = seq(0, 50, 5),
  las = 3,
  labels = seq(0, 50, 5),
  font.axis = 1,
  cex.axis = 1
)
title(main = "(B)")

b = 0.05
f = 1
plot(
  lrna_n_d_com,
  col = "red",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)
plot(
  lrna_c_d_com,
  col = "blue",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)
plot(
  mrna_n_d_com,
  col = "black",
  lwd = 2,
  add = TRUE,
  do.points = FALSE
)
plot(
  mrna_c_d_com,
  col = "green",
  add = TRUE,
  lwd = 2,
  do.points = FALSE
)

p1 = "lnc_nuc."
p2 = "lnc_cyt."
p3 = "m_nuc."
p4 = "m_cyt."

text(
  1,
  ystart,
  labels = stri_c("lnc-nuc.=", nrow(lrna_n_com)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - b,
  labels = stri_c("lnc-cyt.=", nrow(lrna_c_com)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - 2 * b,
  labels = stri_c("m-nuc.=", nrow(mrna_n_com)),
  adj = 0,
  cex = f
)
text(
  1,
  ystart - 3 * b,
  labels = stri_c("m-cyt.=", nrow(mrna_c_com)),
  adj = 0,
  cex = f
)

#try to use legend
lx = 15
legend(
  x = lx,
  y = 0.9,
  legend = c("lnc-nuc."),
  col = c("red"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = lx,
  y = 0.85,
  legend = c("lnc-cyt."),
  col = c("blue"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)

legend(
  x = lx,
  y = 0.8,
  legend = c("m-nuc."),
  col = c("black"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)
legend(
  x = lx,
  y = 0.75,
  legend = c("m-cyt."),
  col = c("green"),
  lty = 1:1,
  cex = f,
  lwd = 2,
  bty = "n",
  adj = 0,
  horiz = TRUE
)
#try end

# text(35,0.9,labels = "??lnc-nuc.",col = "red",cex=f,family="sans",font = 2,adj = 0)
# text(35,0.8,labels = "??lnc-cyt.",col = "blue",cex=f,family="sans",font = 2,adj = 0)
# text(35,0.7,labels = "??m-nuc.",col = "black",cex=f,family="sans",font = 2,adj = 0)
# text(35,0.6,labels = "??m-cyt.",col = "green",cex=f,family="sans",font = 2,adj = 0)
#

#lines(c(15,40),c(0.65,0.65),lwd=1)
ydown = 0.1
text(15,
     12 * b - ydown,
     labels = "Kolmogorov-Smirnov test",
     adj = 0,
     cex = f)
lines(c(15, 50), c(11 * b - ydown, 11 * b - ydown), lwd = 1)
text(15,
     10 * b - ydown,
     labels = "P1",
     adj = 0,
     cex = f)
text(25,
     10 * b - ydown,
     labels = "P2",
     adj = 0,
     cex = f)
text(35,
     10 * b - ydown,
     labels = "D",
     adj = 0,
     cex = f)
text(42,
     10 * b - ydown,
     labels = "p(same)",
     adj = 0,
     cex = f)
lines(c(15, 50), c(9 * b - ydown, 9 * b - ydown), lwd = 1)
text(15,
     8 * b - ydown,
     labels = p1,
     adj = 0,
     cex = f)
text(25,
     8 * b - ydown,
     labels = p2,
     adj = 0,
     cex = f)
text(
  35,
  8 * b - ydown,
  labels = format(lrna_nc_com$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  8 * b - ydown,
  labels = format(lrna_nc_com$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     7 * b - ydown,
     labels = p3,
     adj = 0,
     cex = f)
text(25,
     7 * b - ydown,
     labels = p4,
     adj = 0,
     cex = f)
text(
  35,
  7 * b - ydown,
  labels = format(mrna_nc_com$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  7 * b - ydown,
  labels = format(mrna_nc_com$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     6 * b - ydown,
     labels = p1,
     adj = 0,
     cex = f)
text(25,
     6 * b - ydown,
     labels = p3,
     adj = 0,
     cex = f)
text(
  35,
  6 * b - ydown,
  labels = format(lrna_mrna_n_com$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  6 * b - ydown,
  labels = format(lrna_mrna_n_com$p.value, digits = 3),
  adj = 0,
  cex = f
)

text(15,
     5 * b - ydown,
     labels = p2,
     adj = 0,
     cex = f)
text(25,
     5 * b - ydown,
     labels = p4,
     adj = 0,
     cex = f)
text(
  35,
  5 * b - ydown,
  labels = format(lrna_mrna_c_com$statistic, digits = 3),
  adj = 0,
  cex = f
)
text(
  42,
  5 * b - ydown,
  labels = format(lrna_mrna_c_com$p.value, digits = 3),
  adj = 0,
  cex = f
)
lines(c(15, 50), c(4 * b - ydown, 4 * b - ydown), lwd = 1)