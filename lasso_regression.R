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
library(lars)
 
#for deep learning-based regression for both mRNA and lncRNA
RNA = c("mRNA", "lncRNA")

for (xrna in RNA)
{
  if (xrna == "mRNA")
  {
    load(file = "mRNA_deeplearning.RData")
    load(file = "mRNA_deeplearning_cor_fdr.RData")
  } else
  {
    load(file = "lncRNA_deeplearning.RData")
    load(file = "lncRNA_deeplearning_cor_fdr.RData")
  }
  sim_cor = c()
  xvalue = c(0.1, 0.05, 0.01) #pvalue or FDR value
  TM = TN
  for (xp in xvalue)
  {
    if (xrna == "mRNA")
    {
      xx = which(as.numeric(as.character(TN_cor$fdr)) < xp)
    } else
    {
      xx = which(as.numeric(as.character(TN_cor$pvalue)) < xp)
    }
    TN_cor1 = TN_cor[xx,]
    TN = TM[, c(1, as.numeric(as.character(TN_cor1$no)))]
    
    
    rnum = sample(nrow(TN))
    TN =as.data.frame(TN[rnum, ])
    
    xTN = matrix(
      data = as.numeric(as.matrix(TN)),
      nrow = nrow(TN),
      ncol = ncol(TN),
      byrow = FALSE
    )
    
    xTN = as.data.frame(xTN)
    colnames(xTN) = colnames(TN)
    rownames(xTN) = rownames(TN)
    TN=xTN
    
    
    
    
    xmodel = lm(DF2formula(TN), data = TN)
    pred = predict(xmodel, TN)
    train_info = data.frame(real = TN[,1],
                       pred = pred)
     
    T=summary(xmodel)
    
    
    xx1 = cor.test(
      as.numeric(as.character(train_info$pred)),
      as.numeric(as.character(train_info$real)),
      alternative = c("two.sided", "less", "greater")[1],
      method = c("pearson", "kendall", "spearman")[1],
      exact = FALSE
    )
    xx2 = cor.test(
      as.numeric(as.character(train_info$pred)),
      as.numeric(as.character(train_info$real)),
      alternative = c("two.sided", "less", "greater")[1],
      method = c("pearson", "kendall", "spearman")[3],
      exact = FALSE
    )
    sim_cor = rbind(sim_cor, c(xp, length(xx),T$adj.r.squared,  xx1$estimate, xx1$p.value,
                               xx2$estimate, xx2$p.value))
  }
  
  sim_cor = as.data.frame(sim_cor)
  
  if (xrna == "mRNA")
  {
    colnames(sim_cor) = c("FDR","Features", "AdjustR2",  "pearson_cor", "pearson_pvalue",
                          "spearman_cor","spearman_p")
    write.csv(file = "mrna_training_lm_regressionmodel.csv", sim_cor)
  }else
  {
    colnames(sim_cor) = c("Pvalue", "Features","AdjustR2", "pearson_cor", "pearson_pvalue",
                          "spearman_cor","spearman_p")
    write.csv(file = "lncrna_training_lm_regressionmodel.csv", sim_cor)
  }
}

 