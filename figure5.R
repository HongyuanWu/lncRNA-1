library(RColorBrewer)
library(corrplot)

# Correlation matrix
load(file = "length_gc_rna2d_cor.Rdata")

# p-value matrix
load(file = "length_gc_rna2d_p.Rdata")

# Color palette
col <- colorRampPalette(c("red", "white", "blue"))(8)
 
corrplot(
  r.mat,
  method = "circle",
  type = "full",
  #add=TRUE,
  col = col,
  title = "Correlation between RNAs half-lives and their lengths, GC contents, and RNA2Ds",
  bg = "white",
  outline = FALSE,
  p.mat = p.mat,
  mar=c(0, 0, 1, 0),
  
  sig.level = 0.01,
  insig = "pch", # "blank" or "label_sig" or "pch" or "p-value"
  pch = 4, pch.col = "black", pch.cex = 2,
  
  
  
  tl.col = "black",
  cl.pos = "r",
  cl.length = 6,
 
  tl.offset = 1.3,
  
  win.asp = 1.1,
  
  tl.cex = 1,
  tl.srt = 90,
  rect.lwd = 1,
  na.label = "na",
  is.corr = FALSE
)
text(2,3.9,"lnc-human1",cex = 1.2)
text(6,3.9,"lnc-human2",cex=1.2)
text(9,3.9,"lnc-human",cex=1.2)

text(12,3.9,"m-human1",cex=1.2)
text(16,3.9,"m-human2",cex=1.2)
text(20,3.9,"m-human",cex=1.2)

 