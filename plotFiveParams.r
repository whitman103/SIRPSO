Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.27/bin/gswin64c.exe")
library(extrafont)
library(Cairo)
graphics.off()



deterSwitch=FALSE
if(deterSwitch){
	baseFolder="D:\\Downloads\\12_02_2020\\DataFolder_ODEMeansDeter_2\\" 
	precursor="noExtrinsicNoise"
	trueString=paste(baseFolder,precursor,"_outRunge_noNoise.txt",sep="")
}
if(!deterSwitch){
	baseFolder="D:\\Downloads\\12_04_2020\\DataFolder_ODEMeansOnly_1\\" 
	precursor="extrinsicNoise"
	trueString=paste(baseFolder,precursor,"_outRunge_testNoise.txt",sep="")
}

bounds=c(1.38e-6,5e-3,8.3e-4,3,8.3e-4,3,1.38e-5,5e-2,2.7e-3,10,2.7e-3,10)

noiseString=paste(baseFolder,precursor,"_bestParticles.txt",sep="")
trueString=paste(baseFolder,precursor,"_outRunge_testNoise.txt",sep="")
inData=read.table(noiseString,header=FALSE)
trueParameters=read.table(trueString,header=FALSE)
colorVector=colormap()
minFitness=min(inData[,1])
maxFitness=max(inData[,1])



plot(inData[,2],rep(c(trueParameters[1]),times=length(inData[,1])),col=plotColors,pch=2)

points(inData[,2],rep(c(trueParameters[1]*1.1),times=length(inData[,1])),col=c(0,1))

outNoise=paste(baseFolder,precursor,"Noise.pdf",sep="")
#pdf(outNoise,onefile=FALSE)
CairoPDF(outNoise,6.5,7.5)
par(family="CMU Serif")
plotParameters=matrix(1:12,nrow=6,ncol=2)
noiseString=paste(baseFolder,precursor,"_bestParticles.txt",sep="")

outMeans=numeric(6)
trueValues=numeric(6)
inData=read.table(noiseString,header=FALSE)
inData=inData[,2:ncol(inData)]
trueData=read.table(trueString,header=FALSE)
par(family="CMU Serif")
for(param in 1:6){
	trueValues[param]=trueData[1,param]
	xMin=(min(median(inData[,param]),trueValues[param]))*(.9)
	plotParameters[param,1]=xMin
	xMax=max(median(inData[,param]),trueValues[param])*1.1
	xMin=bounds[2*param-1]
	xMax=bounds[2*param]
	plotParameters[param,2]=xMax
	roundNumber=3
	plotColors=numeric(length(inData[,1]))
	for(i in 1:length(inData[,1])){
		plotColors[i]=colorVector[floor(length(colorVector)*inData[i,1]/(maxFitness-minFitness))+1]
	}

	if(param==1||param==4){
		printLabels=formatC(xMin*10^((0:5)/5*log10(xMax/xMin)),format="e",digits=3)
	}
	else{
		printLabels=round(xMin*10^((0:5)/5*log10(xMax/xMin)),printAccuracy)
	}
	plot(median(inData[2:length(inData[,param]),param]),param+.12,col='black',cex=2,axes=FALSE,ylab="",ylim=c(.8,6.2),xlab="Parameter Values",xlim=c(xMin,xMax),log="x")
	points(trueValues[param],param+.12,col='red',cex=2)
	
	outMeans[param]=mean(inData[2:length(inData[,param]),param])
	par(cex=.8)
	printAccuracy=3
	labelPoints=xMin*10^((0:5)/5*log10(xMax/xMin))
	axis(pos=(param-.12),side=1,at=labelPoints,labels=printLabels)
	par(cex=1)
	
	points(inData[2:length(inData[,param]),param],replicate(length(inData[,param])-1,param),col=colorVector)
	shift=.4
	polygon(c(xMin,xMin,xMax,xMax), c(param,param+shift,param+shift,param),
    col = rgb(param/10,0,1-param/10,.25), border = NA)
	
	par(new=TRUE)
}
axis(side=2, labels=c("k0","k1","k2","k3","k4","k5"),at=c(1:6),pos=0.0015)
title(main="PSO Results, Extrinsic Noise Included\nMeans Used Only")
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

