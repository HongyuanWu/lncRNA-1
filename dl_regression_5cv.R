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
  xvalue = c(0.1, 0.05, 0.01)  #p or FDR value
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
    
    sam_index = rep(1:5, round(nrow(TN) / 5) + 1)[1:nrow(TN)]
    
    test_cv = c()
    for (i in 1:5)
    {
      tn_test = which(as.numeric(as.character(sam_index)) == i)
      TN_test = TN[tn_test,]
      TN_train = TN[-tn_test,]
      
      train_x = TN_train[, 2:ncol(TN)]
      train_y = as.numeric(as.character(TN_train[, 1]))
      test_x = TN_test[, 2:ncol(TN)]
      test_y = as.numeric(as.character(TN_test[, 1]))
      
      train_x = matrix(
        data = as.numeric(as.matrix(train_x)),
        nrow = nrow(train_x),
        ncol = ncol(train_x),
        byrow = FALSE
      )
      
      test_x = matrix(
        data = as.numeric(as.matrix(test_x)),
        nrow = nrow(test_x),
        ncol = ncol(test_x),
        byrow = FALSE
      )
      
      train_center = apply(train_x, 2, mean)
      train_sd = apply(train_x, 2, sd)
      
      train_y = as.numeric(as.character(train_y))
      
      train_x = scale(train_x)
      
      test_x = scale(
        matrix(
          data = as.numeric(as.matrix(test_x)),
          nrow = nrow(test_x),
          ncol = ncol(test_x),
          byrow = FALSE
        ),
        center = train_center,
        scale = train_sd
      )
      
      test_y = as.numeric(as.character(test_y))
      
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
        epochs = 200,
        batch_size = 1000,
        validation_split = 0.2
      )
      
      #  train_info = data.frame(pred = as.vector(model %>% predict(train_x)),
      #                          real = train_y)
      
      test_info = data.frame(pred = as.vector(model %>% predict(test_x)),
                             real = test_y)
      
      test_cv = rbind(test_cv, test_info)
    }
    
    test_cv = as.data.frame(test_cv)
    xx1 = cor.test(
      as.numeric(as.character(test_cv$pred)),
      as.numeric(as.character(test_cv$real)),
      alternative = c("two.sided", "less", "greater")[1],
      method = c("pearson", "kendall", "spearman")[1],
      exact = FALSE
    )
    
    xx2 = cor.test(
      as.numeric(as.character(test_cv$pred)),
      as.numeric(as.character(test_cv$real)),
      alternative = c("two.sided", "less", "greater")[1],
      method = c("pearson", "kendall", "spearman")[3],
      exact = FALSE
    )
    sim_cor = rbind(sim_cor,
                    c(
                      xp,
                      ncol(TN)-1,
                      xx1$estimate,
                      xx1$p.value,
                      xx2$estimate,
                      xx2$p.value
                    ))
  }
  
  sim_cor = as.data.frame(sim_cor)
  
  if (xrna == "mRNA")
  {
    colnames(sim_cor) = c(
      "FDR",
      "Features",
      "pearson_cor",
      "pearson_pvalue",
      "spearman_cor",
      "spearman_p"
    )
    write.csv(file = "mrna_5cv_deeplearning_regressionmodel.csv", sim_cor)
  } else
  {
    colnames(sim_cor) = c(
      "Pvalue",
      "Features",
      "pearson_cor",
      "pearson_pvalue",
      "spearman_cor",
      "spearman_p"
    )
    write.csv(file = "lncrna_5cv_deeplearning_regressionmodel.csv", sim_cor)
  }
}
