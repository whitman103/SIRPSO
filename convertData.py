baseFileName="NKData//NKData_"
times=[8,32,64,128,256]
ending="min.csv"
inData=open(baseFileName+str(times[0])+ending,'r')

#
for line in inData:
