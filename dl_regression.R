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
    TN = TN[rnum, ]
    TN_train = TN
    
    
    train_x = TN_train[, 2:ncol(TN)]
    train_y = as.numeric(as.character(TN_train[, 1]))
    
    train_x = matrix(
      data = as.numeric(as.matrix(train_x)),
      nrow = nrow(train_x),
      ncol = ncol(train_x),
      byrow = FALSE
    )
    
    train_y = log(as.numeric(as.character(train_y)))
    
    train_x = scale(train_x)
    
    model = keras_model_sequential() %>%
      layer_dense(
        units = 500,
        activation = "relu",
        input_shape = dim(train_x)[2]
      ) %>%
      layer_dropout(rate = 0.2) %>%
      
      layer_dense(units = 250, activation = "relu") %>%
      layer_dropout(rate = 0.2) %>%
      
      layer_dense(units = 100, activation = "relu") %>%
      layer_dropout(rate = 0.2) %>%
      
      layer_dense(units = 50, activation = "relu") %>%
      layer_dropout(rate = 0.2) %>%

      layer_dense(units = 30, activation = "relu") %>%
      layer_dropout(rate = 0.2) %>%
      
      layer_dense(units = 10, activation = "relu") %>%
      layer_dense(units = 1)
    
    
    model %>% compile(optimizer = "rmsprop",
                      loss = "mse",
                      metrics = c("mae"))
    
    history = model %>% fit(
      train_x,
      train_y,
      shuffle = TRUE,
      epochs = 100,
      batch_size = 1000,
      validation_split = 0.2
    )
    
    train_info = data.frame(pred = as.vector(model %>% predict(train_x)),
                            real = train_y)
    
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
    sim_cor = rbind(sim_cor, c(xp, length(xx), xx1$estimate, xx1$p.value,
                               xx2$estimate, xx2$p.value))
  }
  
  sim_cor = as.data.frame(sim_cor)
  
  if (xrna == "mRNA")
  {
    colnames(sim_cor) = c("FDR","Features",  "pearson_cor", "pearson_pvalue",
                          "spearman_cor","spearman_p")
    write.csv(file = "mrna_training_deeplearning_regressionmodel.csv", sim_cor)
  }else
  {
    colnames(sim_cor) = c("Pvalue", "Features", "pearson_cor", "pearson_pvalue",
                          "spearman_cor","spearman_p")
    write.csv(file = "lncrna_training_deeplearning_regressionmodel.csv", sim_cor)
  }
}

 