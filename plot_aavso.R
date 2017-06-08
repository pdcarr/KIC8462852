# this is designed to plot AAVSO magnitudes after some data cleaning, with decreasing magnitude

# housekeeping
source("data_funcs.R") # a file of supporting functions
rm(allFits)
rm(cleanBand)
options(digits=12) # hard to read JDs without this setting

allFits <- list()


##################
# input parameters
llightcurve_name <- "aavsodata_28Mar2017.csv"
maxairmass = 4.0 # air mass values above this will be filtered, as well as missing air masses set >= 100 to turn this off
maxuncertainty = 0.1  # maximum AAVSO uncertainty estimate
wildsd = 3.0 # worst number of standard deviations from mean allowed

earliestJD = 2457300 # only data on or after this JD will be used
ExclCodes <- c("None") # observers not be used in the fit
allBands <- data.frame(bandinQ=c("I","R","V","B"),plotColor=c("#FF00FF","red","green","blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)

plotExcluded <- TRUE
####################

# load light curve
lightcurve <- read.csv(file=llightcurve_name,header=TRUE,check.names=TRUE)
# lightcurve, which is the AAVSO data read in via read.table with header=TRUE
totalRec = length(lightcurve$JD)
cat("\n\nTotal records read from file: ",totalRec,"\n\n")

#  myPlotTitle:  string for plot title
myPlotTitle = paste("AAVSO Data","Airmass <= ",toString(maxairmass),ExclCode," excluded",sep=" ")


# make a list of every record with a Magnitude reported
goodMags <- !is.na(lightcurve$Magnitude)


# loop over the desired bands
numBands = length(allBands$bandinQ)
cleanBand <- matrix(nrow=totalRec,ncol=numBands)
allFits <- list()

index = 1


#clean and separate and the bands in question
for (thisBand in allBands$bandinQ) {
	print(thisBand)	
#	print(index)
	cleanBand[,index] <- cleanAAVSO(lightcurve,thisBand,ExclCodes,maxairmass,maxuncertainty,wildsd,earliestJD)
	allFits <- rbind(allFits,lm(lightcurve[cleanBand[,index],"Magnitude"] ~ lightcurve[cleanBand[,index],"JD"], weights= 1/lightcurve[cleanBand[,index],"Uncertainty"]))

	if(index==1) {
		allClean <- cleanBand[,index]
	} else {
		allClean <- cleanBand[,index] | allClean
	}
	index = index + 1	
}


# axis limits

myxlims = c(min(lightcurve[allClean,"JD"],na.rm=TRUE),max(lightcurve[allClean,"JD"],na.rm=TRUE))

# calculate pretty y limits

myYlims = c(ceiling(max(lightcurve[allClean,"Magnitude"],na.rm=TRUE)*10)/10,floor(10*min(lightcurve[allClean,"Magnitude"],na.rm=TRUE))/10) # set up Y limits for reversed Y axis

# plot the cleaned data, the fit lines and the excluded points
icol=1
for (thisBand in allBands$bandinQ) {
	if (icol > 1){par(new=TRUE)}
	ourCleanData <- cleanBand[,icol]
	print(allBands$plotColor[icol])
	plot(lightcurve[ourCleanData,"JD"],lightcurve[ourCleanData,"Magnitude"],col=allBands$plotColor[icol],xlab="Julian Date",ylab="Magnitude",xlim= myxlims,ylim = myYlims,main=myPlotTitle,pch=3) # plot band in prple
	grid()
	
	# plot line fit
	myslope <- coefficients(allFits[icol,])[2]
	basemag <- coefficients(allFits[icol,])[1]
	
	par(new=TRUE)
	curve(basemag + myslope*x,from=min(lightcurve[ourCleanData,"JD"],na.rm=TRUE),to=max(lightcurve[ourCleanData,"JD"],na.rm=TRUE), add=TRUE)
	
	if(plotExcluded) {
		# plot excluded points in black
		IResid <- !ourCleanData & goodMags & (lightcurve$Band == allBands$bandinQ[icol])
		par(new=TRUE)
		
		plot(lightcurve[IResid,"JD"],lightcurve[IResid,"Magnitude"],col="black",xlab="Julian Date",ylab="Magnitude",xlim=myxlims,ylim = myYlims,pch=20) # plot I removed points in black
	}	
	icol = icol + 1
}



# who the observers were for cleaned data:
#buncha <- AAVSOObsStats(lightcurve,cleanBand,allBands)
#cat("\n\nObserver Summary\n")
#print(buncha)

# summaries of the fits
cat("\n\n Observer Code: ", ExclCode," observations not included\n")
cat(length(lightcurve$JD),"total observations loaded\n")
cat(length(lightcurve[cleanI | cleanR | cleanV | cleanB,"JD"]),"observations after cleaning\n")

cat("\n\nObserver Stats for this set of fits")
AAVSOObsStats(lightcurve,cleanBand,allBands)

for (index in 1:numBands) {
	cat("\n\n",allBands$bandinQ[index]," Summary")
	print(summary(allFits[index,]))
}


