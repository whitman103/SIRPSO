baseFolder="D:\\SIRPSO\\interactionFolder\\"
currentBase="bindingTwoStateGenes_base.txt"
library(extrafont)
par(family="CMU Serif")

currentInteraction="bindingTwoStateGenes.txt"
inData2=read.table(paste(baseFolder,currentInteraction,sep=""))
hist(inData2[,8],breaks=40,col=rgb(1,0,0,.5),xlab="Protein Number",ylab="Counts")

inData=read.table(paste(baseFolder,currentBase,sep=""))
hist(inData[,8],breaks=40,add=T,col=rgb(0,0,1,.5))

legend(x="topright",c("Base Three State Gene","Three State With Interaction"),col=c(rgb(1,0,0,.5),rgb(0,0,1,.5)),pch=7)