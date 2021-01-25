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
#===============================================
#network-based comprehensive analysis of lncRNAs
#===============================================
#loading dataset TN
load(file = "lncrna_network.Rdata")

#calculating spearman correlation
TN_cor=c()
for (i in 1:(ncol(TN)-1))
{
  for (j in (i+1):ncol(TN))
  {
    xx=cor.test(as.numeric(as.character(TN[,i])),
                as.numeric(as.character(TN[,j])),
                alternative = "two.sided",
                method = c("pearson", "kendall", "spearman")[3],
                exact = FALSE)
    TN_cor=rbind(TN_cor,c(i,colnames(TN)[i],
                          j,colnames(TN)[j],
                          xx$estimate,
                          xx$p.value))
  }
}
TN_cor=as.data.frame(TN_cor)
colnames(TN_cor)=c("from_node","first_name",
                   "to_node","second_name",
                   "rho","p")
fdr=p.adjust(as.numeric(as.character(TN_cor$p)),
             method = "bonferroni",
             n=nrow(TN_cor))
TN_cor=data.frame(TN_cor,fdr=fdr)


#Display network
xx=which(as.numeric(as.character(TN_cor$fdr))<0.01)
TN_cor=TN_cor[xx,]

#set edge property
half_edge=data.frame(from=as.numeric(as.character(TN_cor$from_node)),
                     to=as.numeric(as.character(TN_cor$to_node))
)

edge_color=c()
edge_width=c()
for (i in 1:nrow(half_edge))
{
  if (as.numeric(as.character(TN_cor$rho))[i] >=0)
  {
    edge_color=c(edge_color,"red")
  }else
  {
    edge_color=c(edge_color,"black")
  }
  
  if (as.numeric(as.character(TN_cor$p))[i]<=1.0e-33)
  {
    edge_width=c(edge_width,3)
  }
  
  if ((as.numeric(as.character(TN_cor$p))[i]<=1.0e-26)  &
      (as.numeric(as.character(TN_cor$p))[i]>=1.0e-33))
  {
    edge_width=c(edge_width,2)
  }
  
  if ((as.numeric(as.character(TN_cor$p))[i]<=1.0e-2)  &
      (as.numeric(as.character(TN_cor$p))[i]>=1.0e-26))
  {
    edge_width=c(edge_width,1)
  }
}



half_node=data.frame(id=1:ncol(TN),
                     label=colnames(TN),
                     shape = rep("ellipse",ncol(TN)),
                     level = c(1,2,2,3,3,4),
                     borderWidth=2,
                     physics=rep(FALSE,ncol(TN)),
                     x=c(300,rep(c(100,450),2),300),  # x coordinate
                     y=c(50,rep(100,2),rep(150,2),200),  # y coordinate
                     color=c("red",rep("lightgreen",ncol(TN)-1)),
                     shadow = rep(TRUE,ncol(TN))
)

half_edge=data.frame(half_edge,
                     color=edge_color,
                     width=edge_width)

visNetwork(half_node,
           half_edge,
           main = "half-life regulation network for lnc-human1")%>%
  visHierarchicalLayout(levelSeparation = 80,
                        nodeSpacing = 200)
