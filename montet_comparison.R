# this is designed to plot AAVSO magnitudes after some data cleaning, with decreasing magnitude, and to fit the data to a straight line.

# housekeeping
#compile the function libraries
source("data_funcs.R") # a file of supporting functions
source("astro_funcs.R")
# remove old garbage if it's there
rm(allFits)
rm(cleanBand)
rm(binCurve)
rm(mars)
# options
options(digits=12) # hard to read JDs without this setting
# load the required packages
library("MASS") # for rlm() and lqs()
library("smooth") # for smoothing
library("earth") # for MARS
library("crayon") # to add color to text

allFits <- list()


##################
# input parameters
llightcurve_name <- "aavso_05Jun2017.csv"
mands.curve.name <- "Montet_SIMON_FFI_table.csv"
maxairmass <- 1.5 # air mass values above this will be filtered out, as well as missing air masses. Set >= 100 to turn this off
maxuncertainty <- 0.02  # maximum AAVSO uncertainty estimate
maxBinUncertainty <- 0.2 # worst standard deviation to accept for a binned set of observations
wildsd <- 10.0 # worst number of standard deviations from mean allowed

earliestJD = 2457294 # only data on or after this JD will be used
#earliestJD <- 2457700
#earliestJD <- 2457800
startPlot <- earliestJD
#startPlot <- 2457800
plotRelTimes <- TRUE
#####################
mands.JD.base <- 2454833
mands.shift <- earliestJD - mands.JD.base
extra.shift <- -640 # days
##########
includeExclude <- FALSE # TRUE if your list of observer codes is to to be included, FALSE if excluded or not used
ExclCodes <- "None"
#ExclCodes <- c("DUBF","SGEA","ELYA","OAR","AAM","ATE","HJW","BALB","BMAK","BPAD","BSM","LPB","GKA","HBB")
#ExclCodes <- c("LDJ","SGEA","ELYA","OAR","AAM","ATE","HJW","BALB","BMAK","BPAD","BSM","LPB","GKA","HBB")
#ExclCodes <- c("SGEA","ELYA","OAR","AAM","ATE","HJW","BALB","BMAK","BPAD","BSM","LPB","GKA","HBB")
#ExclCodes <- c("LDJ","DUBF","SGEA","ELYA","OAR","AAM","ATE","HJW","BALB","BMAK","BPAD","BSM","LPB","GKA","HBB")
#ExclCodes <- c("LDJ","DUBF","ELYA","HJW","JM")
#ExclCodes <- c("DUBF","LDJ","OAR","SGEA","DKS","OJJ","LPB","BPAD")
#ExclCodes <- c("DUBF","LDJ","OAR","SGEA","DKS")
#ExclCodes <- c("DUBF","ELYA","LPB","OJJ","HJW","OAR")
#ExclCodes <- c("JM","LDJ","ELYA","DKS","OJJ","OAR","ATE","BPAD","HJW")
#ExclCodes <- c("LDJ","UJHA","DKS","OJJ","JM","DUBF","ELYA","HJW")
#ExclCodes <- c("JM","LDJ","OAR","LPB")
#ExclCodes <- c("DUBF","JM","LDJ","LPB")
#ExclCodes <- c("DUBF","DKS","ELYA","OAR","NRNA","ATE","HJW","BPAD","OJJ","LBG","LDJ","UJHA","OYE","GFRB","OAS","MJB","EEY") # V ensemble
#ExclCodes <- c("DUBF","GKA","BPAD","LPB","SJAR","LBG","LDJ","LWHA") # R ensemble
#ExclCodes <- c("OAR","OJJ","GKA","MJB","SJAR","LWHA","LBG","LPB","LDJ","CMP") # I ensemble
#ExclCodes <- "JM"
#ExclCodes <- "LDJ"
########
plotMee <- NA # do not highlight any particular observer code
#plotMee <- "HJW" # observer code to plot with special character
#plotMee <- "JM"
#plotMee <- "LDJ"
#plotMee <- "DUBF"
#plotMee <- "ELYA"
meeColor <- "darkviolet"
########
allBands <- data.frame(bandinQ=c("I","R","V","B"),plotColor=c("darkviolet","red","green","blue"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("V"),plotColor=c("darkgreen"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("V","B"),plotColor=c("green","blue"), stringsAsFactors=FALSE)
allBands <- data.frame(bandinQ=c("B"),plotColor=c("blue"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I"),plotColor=c("darkviolet"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("R"),plotColor=c("red"), stringsAsFactors=FALSE)
#allBands <- data.frame(bandinQ=c("I","R","B"),plotColor=c("darkviolet","red","blue"), stringsAsFactors=FALSE)

########################
deltaJD <- 4.0 # bin width in days
########################
plotExcluded <- FALSE # set to TRUE to plot the points in the lightcurve not used in the fit.
plotQuadratic <- FALSE # set to TRUE to plot a quadratic fit
generateTS <- TRUE # set to TRUE to creat a time series from the data
tsBinWidth <- 10.0 # time series bin width in days. Important if generateTS is TRUE
smoothTS <- TRUE # set to TRUE to smooth the times series
tsSmoothOrder <- 10 # the order for the moving average to smooth the time series
tryLQS <- FALSE # set to TRUE is you want to try resistant regression.
userlm <- TRUE # set to TRUE to use robust lm, or rlm() if not using MARS.
plotMARS <-  TRUE # set to TRUE to try a MARS fit instead of lm() or rlm()
plotResiduals <- TRUE # set to true to plot the residuals vs. time
plot2Lines <-  FALSE  # two line feature doesn't work well
lqsColor <- "darkgreen"
weightedBins <- FALSE # set to TRUE to weight lower uncertainty bins more.

####### MARS
marsOrder <- 9
marsPenalty <- 2 # set to 0 to avoid penalizing knots in pruning pass
marsPMethod <- "none" # set to "none" to avoid pruning#marsPMethod <- "backward" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning#marsPMethod <- "backward" # set to "none" to avoid pruning
splineRaw <-  FALSE # do the spline on the raw lightcuve, not binned.
##############################
okComparison <- "(000-?BLS-?556)|(000-?BLS-?551)|(000-?BLS-?553)|000-?BLS-?552)|(000-?BLS-?554)|(000-?BLS-?549)|(000-?BLS-?555)|(108)|(113)|(116)|(118)|(121)|(124)|(128)|(ENSEMBLE)|(APASS20062365[+-]442738)" # regular expression from AAVSO photometry table

editUser <- data.frame(obsCode= "UJHA",startJD=2457650,endJD=2457760,band="V",stringsAsFactors=FALSE)
# fill in some missing airmass values
lasCruces <- c(32.31994,-106.763654) # center of Las Cruces, NM in decimal degrees latitude, longitude.
tabbysLoc <- c("+44d 27m 24.61s","20h 06m 15.457s") # right ascension and declination of the star.
missingAirmass <- "JM"	# observer code
### option to draw a vertical date line
drawDateLine <-  TRUE
jdLine <- 2457892.0
jdLineColor <- "red"
jdLineText <- "18May17"
##########################################################################

if (includeExclude) {
	inclWord <- "used"
} else {
	inclWord <- "not used"
}

# load light curve
lightcurve <- read.csv(file=llightcurve_name,header=TRUE,check.names=TRUE)

# lightcurve, which is the AAVSO data read in via read.csv with header=TRUE
totalRec = length(lightcurve$JD)
cat("\n\nTotal records read from file: ",totalRec,"\n\n")

# read in Montet and Simon FFI data:
MandSFFI <- read.csv(file= mands.curve.name,header=TRUE)


#  myPlotTitle:  string for plot title


# make a list of every record with a Magnitude reported
goodMags <- !is.na(lightcurve$Magnitude)

# replace missing airmass for JM
jmtest = lightcurve$Observer_Code == missingAirmass & is.na(lightcurve$Airmass)
lightcurve$Airmass[jmtest] <- AirMass(lightcurve$JD[jmtest],lasCruces, tabbysLoc)

# loop over the desired bands
numBands = length(allBands$bandinQ)
cleanBand <- matrix(nrow=totalRec,ncol=numBands)
allFits <- list()
allQFits <- list()
allResistFit <- list()
tmin <- head(lightcurve$JD,n=1)

# create a vector of shifted times
MandS.shift.times <- MandSFFI$Time + mands.JD.base + mands.shift + extra.shift
fit.times = MandS.shift.times - tmin
MandS.Fit <- earth(x= fit.times,y= MandSFFI$Normalized.Flux,nk= 7,pmethod= "backward",penalty = 2)

index = 1

# edit specific observers over time range(s) and specific passband
editCurve <- ObserverJDEdit(editUser,lightcurve)

#clean and separate and the bands in question
for (thisBand in allBands$bandinQ) {
#	print(thisBand)	
#	print(index)
	# clean the data for this passband
	cleanBand[,index] <- cleanAAVSO3(lightcurve,thisBand,ExclCodes,includeExclude,maxairmass,maxuncertainty,wildsd,earliestJD,okComparison) & editCurve
	index = index + 1
}
	
	#bin the data for each Observer
binCurve <- binAAVSO	(lightcurve,cleanBand,allBands,deltaJD)

# test for less than max uncertainty
uncertaintyTest <- binCurve$Uncertainty <= maxBinUncertainty

index = 1

# loop over the passbands and do regression for each one

for (thisBand in allBands$bandinQ) {
	# fold in test for passband
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
	desmat <- binCurve$JD[btest] - tmin # subtract off the earliest time to get a meaningful intercept
	if(length(desmat) < 2) {
			cat("\n\nWarning: fewer than 2 points in the band:",thisBand,"\n")
			next
		}
	
	# determine the bin weights, if any
	if (userlm & !plotMARS) {
		# use robust algorithm
		thisFit <- rlm(binCurve$Magnitude[btest] ~ desmat,na.action="na.omit",psi=psi.bisquare)
	} else if(!plotMARS){
		if (weightedBins) {
			binWeights <- 1/binCurve[btest,"Uncertainty"]
		} else {
			binWeights <- NULL
		}
		# do the linear regression in the binned data for this band, starting at the earliest time in the file
		thisFit <- lm(binCurve[btest,"Magnitude"] ~ desmat, weights= binWeights)
	}	
	
	# if MARS is selected (plotMARS == TRUE), use it to do the regression
	if(plotMARS) {
		# use MARS algorithm in earth()
		cat("\n Using MARS algorithm\n")
		if(splineRaw) {
			# do the spline fit on the unbinned data
			desmat <- lightcurve$JD[cleanBand[,index]] - tmin
			thisFit <- earth(x=desmat,y=lightcurve$Magnitude[cleanBand[,index]],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)
		} else {
			thisFit <- earth(x=desmat,y=binCurve$Magnitude[btest],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)
	#		print(class(thisFit))
		}
	} 
	
	# output a summary to the console 	
	cat("\n\n Band",thisBand,"summary")
	print(summary(thisFit))
	
	# try resistant regression if selected
	if (tryLQS) {
		resistFit <- lqs(formula = binCurve$Magnitude[btest] ~ desmat)
		cat("\n\n Band",thisBand," resistant fit coefficients")
		cat(resistFit$coefficients,"\n")
		#build up the matrix of resistant fits
		allResistFit <- rbind(allResistFit,resistFit)
	}

	if (plotQuadratic) {
		# fit a quadratic just for fun
		desmat <- outer(binCurve$JD[btest]-tmin,1:2,"^")
		qfit <- lm(binCurve$Magnitude[btest] ~ desmat, weights= binWeights)
		allQFits <- rbind(allQFits,qfit)
	}
	#package the fits into one matrix
	allFits <- rbind(allFits,thisFit)

	if(index==1) {
		allClean <- cleanBand[,index]
	} else {
		allClean <- cleanBand[,index] | allClean
	}
	index = index + 1	
}

##################################### Plot This Stuff ##################################
# axis limits

myxlims = c(startPlot,max(binCurve$JD,na.rm=TRUE))

# calculate pretty y limits

myYlims = c(ceiling(max(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE)*10)/10,floor(10*min(binCurve$Magnitude[uncertaintyTest],na.rm=TRUE))/10) # set up Y limits for reversed Y axis

# set up plot title text
ocodesInTitle <- paste("Observer",head(ExclCodes,n=3),inclWord,sep=" ")
howManyObs = length(unique(binCurve$Observer_Code[uncertaintyTest]))
if (length(ExclCodes) > 3){ocodesInTitle <- c(ocodesInTitle,paste("and",howManyObs - 3,"more observer code(s)",sep=" "))}

myBands = paste(allBands$bandinQ,collapse=" ")
titleString <- c(paste("AAVSO",myBands,"Data with",deltaJD,"Day Bins",sep=" "), ocodesInTitle)
myPlotTitle <- paste(titleString,collapse="\n")


# plot the cleaned and binned data, the fit lines and the excluded points
icol=1
if(plotRelTimes) {
	myTimes <- binCurve$JD - tmin
	jdLine <- jdLine - tmin
	myxlims <- myxlims - tmin
	lcTimes = lightcurve$JD - tmin
	myXLabel <- paste("Julian Date -",as.character(tmin),sep=" ")
} else {
	myTimes <- binCurve$JD
	myXLabel <- "Julian Date"
	}

for (thisBand in allBands$bandinQ) {
#	if (icol > 1){par(new=TRUE)}
	ourCleanData <- cleanBand[,icol]
	btest <- (binCurve$Band == thisBand) & uncertaintyTest
		
	if(icol==1) {
		plot(myTimes[btest],binCurve[btest,"Magnitude"],col=allBands$plotColor[icol],xlab= myXLabel,ylab="Magnitude",xlim= myxlims,ylim = myYlims,main=myPlotTitle,pch=3,cex.main=0.7)
	} else {
		points(myTimes[btest],binCurve[btest,"Magnitude"],col=allBands$plotColor[icol],pch=3)
	}
	
	# plot line fit
	if(plot2Lines) {
		cat(red("plt2Lines doesn't really work. Use earth() instead"))
	} else if(!plotMARS) {
		# plot the lm() or rlm() fit
		myslope <- coefficients(allFits[icol,])[2]
		basemag <- coefficients(allFits[icol,])[1]
		curve(basemag + myslope*(x-tmin),from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col="black")
	} 
	
	########################## plot the MARS fit if that is selected
	
	if (plotMARS) {
		mars <- allFits[icol,]
		intercept <- mars$fitted.values[1] 
		closest.0time <- min(abs(fit.times))
		use.me <- abs(fit.times) <= closest.0time
		if(length(fit.times[use.me]) != 1) {
			cat("\nThere is a problem with finding the unique nearest time to zero for the M&S fit.")
			break
			}
		normalize.fit <- MandS.Fit$fitted.values[use.me]
		MandS.plot.Mags <- -2.5*log10(MandS.Fit$fitted.values/normalize.fit) + intercept
		
		if (splineRaw) {
			lines(x=lcTimes[cleanBand[,icol]],y=mars$fitted.values,col= allBands$plotColor[icol],lwd=2)
		} else {
			lines(x=myTimes[btest],y=mars$fitted.values,col= "black",lwd=2)
			lines(x= fit.times,y= MandS.plot.Mags,col="purple",lwd=2,lty="dotted")
		}
	}
	#optionally plot the LQS fit
	if (tryLQS) {
		myslope <- coefficients(allResistFit[icol,])[2]
		basemag <- coefficients(allResistFit[icol,])[1]
		curve(basemag + myslope*(x-tmin),from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col=lqsColor)
	}
	
	# an option to plot the quadratic fit
	if (plotQuadratic) {
		#quadratic fit
		basemag <- coefficients(allQFits[icol,])[1]
		linterm <- coefficients(allQFits[icol,])[2]
		qterm <- coefficients(allQFits[icol,])[3]
		curve(basemag + linterm*(x-tmin) + qterm*(x -tmin)^2,from=min(binCurve$JD,na.rm=TRUE),to=max(binCurve$JD,na.rm=TRUE), add=TRUE,col="red")
	}
	
#plot the excluded points, if option selected
	if(plotExcluded) {
		# plot excluded points in black
		IResid <- !ourCleanData & goodMags & (lightcurve$Band == allBands$bandinQ[icol])		
		points(lcTimes[IResid],lightcurve[IResid,"Magnitude"],col="black",pch=20) # plot removed points in black
	}	
	icol = icol + 1
}

if (!is.na(plotMee)) {
	imSpecial <- grep(plotMee,binCurve$Observer_Code,ignore.case=TRUE)
	points(binCurve$JD[imSpecial],binCurve$Magnitude[imSpecial],col=meeColor,pch=20,cex=1.5)
}

grid(col="black")

if(drawDateLine) {
	# draw a vertical line for a date of interest
	lines(x=c(jdLine,jdLine),y=myYlims,col=jdLineColor,lwd=1,lty="dashed")
	text(x= jdLine,y=myYlims[2],labels=jdLineText,pos=3,cex=0.5)
}


# generate a time series and plot if called for

if (generateTS) {
	quartz("Time Series")
	tsMain <- "Time Series"	
	myts <- binAAVSO_ts(lightcurve, cleanBand,allBands, tsBinWidth)
	if (smoothTS) {
		tsMain <- paste("Smoothed",tsMain," - Order = ", tsSmoothOrder)
		tsMain 
		okts <- !sapply(X=myts,FUN=is.na,simplify="logical")
		for(iband in 1:numBands) {
			tsScratch <- myts[,iband]
			tsScratch[okts[,iband]] <- sma(data = myts[okts[,iband],iband],order=tsSmoothOrder)$fitted
			if (iband == 1) {
					bts <- tsScratch
				} else {
					bts  <- cbind(bts,tsScratch,deparse.level=0)
				}
		}
	} else {
		bts <- myts[,1:numBands] # not smoothed
	}
	if(numBands > 1 ) {
		colnames(bts) <- sapply(allBands$bandinQ,FUN=paste,"Mag",sep=" ")
		plot.ts(bts,plot.type="multiple",xlim= myxlims,ylim = myYlims,main=tsMain,cex.axis=0.9,lwd=2)
	} else {
		plot(bts,ylim=myYlims,main=tsMain,lwd=2,col=allBands$plotColor,ylab=paste(allBands$bandinQ,"Mag"))
		points(myts,col="grey",pch=20,cex=0.5)
	}
	grid(col="black")
}

# plot the residuals if desired:
if (plotResiduals) {
	quartz("Residuals")
	irow <- 1
	
	for (thisBand in allBands$bandinQ) {
		btest <- (binCurve$Band == thisBand) & uncertaintyTest
		if(irow == 1) {
			plot(myTimes[btest],allFits[irow,]$residuals,col= allBands$plotColor[irow],xlab="Julian Date",ylab="Magnitude",xlim= myxlims,main="Residuals",pch=20,cex.main=1.0)
		} else {
			points(myTimes[btest],allFits[irow,]$residuals,col=allBands$plotColor[irow],pch=20)
		}
		irow <- irow + 1
	}
	grid(col="black")
}

# who the observers were for cleaned data:
buncha <- AAVSOObsStats(lightcurve,cleanBand,allBands)
cat("\n\nObserver Summary - Raw Observations\n")
print(buncha)

probOK <- logical()

for (thisBand in allBands$bandinQ) {probOK <- cbind(probOK,!is.na(binCurve$JD) & uncertaintyTest)}
binda <-AAVSOObsStats(binCurve,probOK,allBands)
cat("\n\n    Observer Summary - Binned Observations with acceptable scatter\n")
print(binda)


# summaries of the fits 

cat("\n\n    Observer Code: ", ExclCodes," observations",inclWord,"\n")
cat("   ",length(lightcurve$JD),"total observations loaded\n")
#cat("    ",length(lightcurve[cleanI | cleanR | cleanV | cleanB,"JD"]),"raw observations after cleaning\n")
cat("    ",length(binCurve$JD[uncertaintyTest]),"binned observations with",deltaJD,"day bins")

#cat("\n\nObserver Stats for this set of fits")
#AAVSOObsStats(lightcurve,cleanBand,allBands)

for (index in 1:numBands) {
	if(!plotMARS) {
		cat("\n\n",allBands$bandinQ[index],"    Summary\n")
		pctPerYear <- (10^(-1*coefficients(allFits[index,])[2]*365.24/2.5)-1)*100
		cat("    ",coefficients(allFits[index,])[2]*365.24*100,"magnitudes per century","or",pctPerYear,"% per year\n")
	}
}


