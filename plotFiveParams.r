Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.26/bin/gswin32c.exe")
library(extrafont)
library(Cairo)
graphics.off()
baseFolder="D:\\Downloads\\10_15_2020\\DataFolder_pVavTests_0\\" 
precursor="NoDynamicsLater"

outNoise=paste(baseFolder,precursor,"Noise.pdf",sep="")
#pdf(outNoise,onefile=FALSE)
CairoPDF(outNoise,6.5,7.5)
plot.new()
par(family="CMU Serif")
plotParameters=matrix(1:12,nrow=6,ncol=2)
noiseString=paste(baseFolder,precursor,"_bestParticles.txt",sep="")
trueString=paste(baseFolder,precursor,"_outRunge_noNoise.txt",sep="")
outMeans=numeric(6)
trueValues=numeric(6)
inData=read.table(noiseString,header=FALSE)
inData=inData[,2:ncol(inData)]
trueData=read.table(trueString,header=FALSE)
par(family="CMU Serif")
for(param in 1:6){
	xMin=(min(inData[,param]))*(.9)
	plotParameters[param,1]=xMin
	xMax=max(inData[,param])*1.1
	plotParameters[param,2]=xMax
	roundNumber=3
	trueValues[param]=trueData[1,param]
	if(param==1){
		printLabels=formatC(seq(xMin,xMax,length.out=5),format="e")
	}
	else{
		printLabels=round(seq(xMin,xMax,length.out=5),printAccuracy)
	}
	
	plot(trueData[1,param],param+.12,xlim=c(xMin,xMax),ylim=c(.8,6.2),axes=FALSE,xlab="Parameter Values",ylab="",col='red',cex=2)
	
	points(mean(inData[2:length(inData[,param]),param]),param+.12,col='black',cex=2)
	outMeans[param]=mean(inData[2:length(inData[,param]),param])
	par(cex=.8)
	printAccuracy=3
	axis(pos=(param-.12),side=1,at=round(seq(xMin,xMax,length.out=5),5),labels=printLabels)
	par(cex=1)
	
	points(inData[2:length(inData[,param]),param],replicate(length(inData[,param])-1,param),col='blue')
	shift=.4
	polygon(c(xMin,xMin,xMax,xMax), c(param,param+shift,param+shift,param),
    col = rgb(param/10,0,1-param/10,.25), border = NA)
	if(param==1){
		text(xMax,param+.25,paste("True Value is: ",formatC(trueValues[param],format='e'),sep=""),pos=2)
	}
	else{
		text(xMax,param+.25,paste("True Value is: ",round(trueValues[param],printAccuracy),sep=""),pos=2)
	}
	
	par(new=TRUE)
}
axis(side=2, labels=c("k0","k1","k2","k3","k4","k5"),at=c(1:6),pos=xMin*.7)
title(main="Nonlinear pVav Response Regime wrt Shp1\nOnly pVav Mean in Cost Function")
dev.off()
embed_fonts(outNoise)
write.table(outMeans,paste(baseFolder,"outFoundParameters.txt",sep=""),sep=" ",row.names=FALSE,col.names=FALSE)
write.table(trueValues,paste(baseFolder,"outTrueParameters.txt",sep=""),sep=" ",row.names=FALSE,col.names=FALSE)


if(file.exists(paste(baseFolder,"trueDists.txt",sep=""))){
	labels=c("T","I","V","R")
	inData=read.table(paste(baseFolder,"trueDists.txt",sep=""),sep=",")
	inData2=read.table(paste(baseFolder,"testDists.txt",sep=""),sep=",")
	inData=inData[,-c(ncol(inData))]
	inData2=inData2[,-c(ncol(inData2))]
	pdf(paste(baseFolder,"trueDists.pdf",sep=""),6.5,14,onefile=TRUE)
	times=inData[1,]
	inData=inData[-c(1),]
	inData2=inData2[-c(1),]
	plot.new()
	par(family="CMU Serif")
	par(mfrow=c(4,1))
	for(times in 1:4){
		hist(inData[,times],breaks=40,col=rgb(0,0,1,.65),xlab=paste("P(",labels[times],")",sep=""),main=paste("Histogram for ",labels[times]," at 10 Minutes",sep=""))
		hist(inData2[,times],breaks=40,col=rgb(1,0,0,.5),add=TRUE)
		legend(x="topright",c("True Distribution","Test Distribution"),col=c(rgb(0,0,1,.65),rgb(1,0,0,.5)),pch=4)
	}
	plot.new()
	par(family="CMU Serif")
	par(mfrow=c(4,1))
	for(times in 1:4){
		hist(inData[,7*4+times],breaks=40,col=rgb(0,0,1,.65),xlab=paste("P(",labels[times],")",sep=""),main=paste("Histogram for ",labels[times]," at 80 Minutes",sep=""))
		hist(inData2[,7*4+times],breaks=40,col=rgb(1,0,0,.5),add=TRUE)
		legend(x="topright",c("True Distribution","Test Distribution"),col=c(rgb(0,0,1,.65),rgb(1,0,0,.5)),pch=4)
	}
	
	dev.off()
	embed_fonts(paste(baseFolder,"trueDists.pdf",sep=""))
}

