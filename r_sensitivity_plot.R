require("lattice")
require("MASS")
require("sensitivity")
require("plyr")
fits <- read.csv("summ.fits.csv",header=T)
names(fits)

fits <- fits[fits$p2 > 0.005,] # take out 
# effect of density dependent mortality
plot(fits$mud.p,fits$loglik,bty="n",pch=1)

minll <- ddply(fits,.(fits$mud.p),summarize,minll=min(loglik)) 

# partial correlation coefficient
pcorrel <- pcc(fits[,c(5:10)],fits[,3],nboot=100,rank=T)

pcorrelvals <-  c(-0.41763190, -0.50512787,0.09241388, -0.97787489,0.88625805,0.01287321)
pcorrelcimin <- c(-0.481859620,-0.600357682,0.001782897, -0.986831647, 0.858332224,-0.071951688)
pcorrelcimax <- c(-0.34688267, -0.39904573,0.17321569, -0.96974943,0.91107564,0.09708653)

tiff("fig_6_sensitvity_plot.tiff", height = 3, width = 4.5, units = 'in', compression="lzw", res=400)
par(cex=0.6,mar=c(4,4,1,1))

plot(c(1:6),pcorrelvals,bty="n",ylim=c(-1,1),xaxt="n",xlab="Parameter",ylab="Partial rank correlation coefficient")
abline(h=0)
segments(x0=c(1,2,3,4,5,6),x1=c(1,2,3,4,5,6),y0=pcorrelcimin,y1=pcorrelcimax)
mtext(c("Larviposition \n rate"
      ,"Adult emergence \n rate"
      ,"Pupal \n density-independent \n mortality rate"
      ,"Days to \n starvation"
      ,"Pupal \n density-dependent \n mortality coef."
      ,"Number of \n starting pupae")
      ,side=1,
      at=c(1.1,1.9,3,4,5,5.9),cex=0.5)

dev.off()
