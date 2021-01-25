x=c((1:4)*100,434,(5:10)*100, 100*seq(15,40,5))
y=c(8,7,6,5,5,4.3,3.7,3.2,2.6,2.3,1.8,0.6,0.2,0.04,0.02,0.01,0)*100

sam=data.frame(number=x,time=y)

par(mfrow=c(1,1))
plot(x,y,type = "o",
     pch=20,
     xaxt = "n",
     yaxt = "n",
     cex.axis =1,
     xlab = "The number of lncRNAs sampled from lnc-human1",
     ylab = "The number of times with p>0.1",
     asp = c(xlim = c(0,max(x)),
             ylim = c(min(y),max(y)),
             grid(nx=10,ny=10)
     ))

axis(2,
     tick = TRUE,
     at = seq(0, 800,200),
     las = 1,
     cex.axis =1,
     labels = seq(0, 800,200))

axis(
  1,
  tick = TRUE,
  at = seq(0,4000,400),
  las = 2,
  cex.axis =1,
  labels = seq(0,4000,400)
)
