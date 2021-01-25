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

load(file = "halflife_codon.Rdata")

codon_name = apply(as.data.frame(colnames(codon[2:ncol(codon)])),
                   1,
                   function(y)
                   {
                     return(stri_sub(y, 1, 3))
                   })
colnames(codon)[2:ncol(codon)] = codon_name


#for correlation analysis
gene_cor = c()
for (i in 2:ncol(codon))
{
  xx = cor.test(as.numeric(as.character(codon[,i])),
                as.numeric(as.character(codon$halflife)),
                alternative = "two.sided",exact = FALSE,
                method = c("pearson", "kendall", "spearman")[3])
  gene_cor = rbind(gene_cor, c(xx$estimate,xx$p.value))
}

gene_cor = as.data.frame(gene_cor)

gene_fdr=p.adjust(gene_cor$V2,
                  method ="bonferroni",
                  n=nrow(gene_cor))

gene_cor=data.frame(cor=gene_cor$rho,p=gene_cor$V2,FDR=gene_fdr)
rownames(gene_cor)=colnames(codon)[2:ncol(codon)]

gene_cor=gene_cor[order(as.numeric(as.character(gene_cor$cor)),
                        decreasing = TRUE),]

xx = intersect(which(gene_cor$cor > 0), which(gene_cor$FDR <= 0.01))
yy = intersect(which(gene_cor$cor < 0), which(gene_cor$FDR <= 0.01))

xx = which(as.numeric(as.character(gene_cor$cor)) >= 0)
barcolor=rep("grey",64)
cor=as.numeric(as.character(gene_cor$cor))
fdr=as.numeric(as.character(gene_cor$FDR))
for (i in 1:nrow(gene_cor))
{
  if ((cor[i]>0)  & (fdr[i]<0.01))  barcolor[i]="red"
  if ((cor[i]<0)  & (fdr[i]<0.01))  barcolor[i]="green"
}

par(mfrow=c(1,1))

plot(0,0,
     main = "correlation between mRNA codons frequencies and their half-lives",
     cex.axis =1,
     col = "white",
     font = 1,
     xaxt = "n",
     yaxt = "n",
     ylim = c(-0.1,0.1),
     xlim = c(1,80),
     ylab = "Correlation",
     asp = c(grid(nx=10,ny=10)
     ))


#plot
xb = barplot(
  cor,
  col = barcolor,
  xaxt = "n",
  yaxt = "n",
  add = TRUE
)


axis(
  2,
  at = seq(-0.1, 0.1, 0.025),
  las = 1,
  labels = seq(-0.1, 0.1, 0.025),
  font.axis = 1,
  cex.axis = 1
)



axis(
  1,
  tick = FALSE,
  at = xb,
  las = 3,
  cex.axis =0.8,
  #family="serif",
  font = 1,
  adj = 0,
  labels = rownames(gene_cor),
)

rect(
  xleft = xb[1] - 0.5,
  ybottom = 0.0,
  xright = (xb[21]+xb[22]) / 2,
  ytop = 0.095
)
text(
  x = xb[9],
  y = 0.08,
  labels = "FDR<0.01",
  las=2,
  cex = 1.5
)

rect(
  xleft = (xb[48] + xb[49]) / 2,
  ybottom = 0.0,
  xright = xb[64] + 0.5,
  ytop = -0.095
)

text(
  x = xb[55],
  y = -0.07,
  labels = "FDR<0.01",
  cex = 1.5
)
