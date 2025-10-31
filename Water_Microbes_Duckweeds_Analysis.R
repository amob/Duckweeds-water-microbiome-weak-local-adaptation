#############
#Analyze data from local adaptation experiment factorially recombining:
## duckweeds, microbiomes from each duckweed site, water from each duckweed site

#data analysis includes
#effects of inoculation with microbes vs uninoculated plants
#effects of matched source site materials on duckweed growth and traits, microbiome growth
#basic model of source effects for main effects, regardless of local adaptation 


# ASSUMES THE WORKING DIRECTORY CONTAINS THE INPUT FILES
##################

#####
#libraries and functions
library(MCMCglmm) #package for linear models, others appropriate, but functions are different (i.e. code will break without this package)

std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))} #the standard error function

range01=function(x){ # rescale a vector to fall within 0 and 1, useful for plotting with color assignments based on a variable
	newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
	return(newnums)
}

ssbyvar <- function(response,category.vec){ #sums of squares function
                means <- tapply(response,category.vec,mean,na.rm=T) #take the means by category
                ssresid <- sum(sapply(sort(unique(category.vec)), function(z) sum( (response[category.vec==z] - means[names(means)==z])^2,na.rm=T ))) #square of difference of each datapoint from its associated treatment mean (residual variation)
                sstot <- sum((response-mean(response,na.rm=T))^2,na.rm=T) #square of difference of each datapoint from the grand mean (total variation)
                sst <- (sstot-ssresid) # total variation - residual variation = treatment variation
                return(sst/sstot) # treatment variance as a fraction of total variation
                }



#####
#read in data files
finaldata <- read.csv("FinalTimepoint.csv",header=T,stringsAsFactors=F)	
startdata <- read.csv("FirstPhotoTimepoint.csv",header=T,stringsAsFactors=F)	
treatmentdata <- read.csv("Treatments.csv",header=T,stringsAsFactors=F)#includes area but not traits
OD600 <- read.csv("OD Data - 600 Readings.csv",header=T,stringsAsFactors=F)
OD600sortvector <- sapply(1:nrow(treatmentdata), function(z) which(OD600$Plate==treatmentdata$Plate[z] & OD600$Column==treatmentdata$Column[z] & OD600$Row==treatmentdata$Row[z]  ) )
OD600sort <- OD600[OD600sortvector,]
colnames(OD600sort)[7:9] <- c("OD600", "OD600.BlankAverage", "OD600.Corrected") #further columns are not useful
OD420 <- read.csv("OD Data - 420 Readings.csv",header=T,stringsAsFactors=F)
OD420sortvector <- sapply(1:nrow(treatmentdata), function(z) which(OD420$Plate==treatmentdata$Plate[z] & OD420$Column==treatmentdata$Column[z] & OD420$Row==treatmentdata$Row[z]  ) )
OD420sort <- OD420[OD420sortvector,]
colnames(OD420sort)[7:9] <- c("OD420", "OD420.BlankAverage", "OD420.Corrected") #further columns are not useful

#check that all are organized the same (sum to 960, the length of each file in rows)
sum(paste(OD600sort$Plate, OD600sort$Row, OD600sort$Column) == paste(treatmentdata$Plate,treatmentdata$Row,treatmentdata$Column))
sum(paste(OD420sort$Plate, OD420sort$Row, OD420sort$Column) == paste(treatmentdata$Plate,treatmentdata$Row,treatmentdata$Column))
sum(finaldata$area == treatmentdata$Area)
 
#convert area to sqmm.
#as measured in pictures, a plate width is 2035 pixels, plate width in mm is 85.4
plateratio <- 2035/85.4
finaldata$sqmm <- finaldata$area/(plateratio^2)
startdata$sqmm <- startdata$area/(plateratio^2)
dsqmm <- finaldata$sqmm - startdata$sqmm
rsqmm <- (finaldata$sqmm - startdata$sqmm)/startdata$sqmm
rsqmm[is.infinite(rsqmm)] <- NA
cor(cbind(dsqmm,rsqmm,startdata$sqmm,finaldata$sqmm),use="complete.obs") #correlations among options for measuring growth
#note that "start" was not actually day 0, and some photos are missing for this data, we have therefore not used it in analyses

#put data and treatments together in one table, add useful vectors describing treatments
treatfinal <- as.data.frame(cbind(finaldata,treatmentdata,OD600sort[,7:9],OD420sort[,7:9]))
treatfinal$dsqmm <- dsqmm
treatfinal$rsqmm <- rsqmm
treatfinal$initialsqmm <- startdata$sqmm
treatfinal$isinoc <- ifelse(treatfinal$Microbe != "None", "yes","no")
#three optical density tests are outside the reasonable range for both wavelengths, and suggest a plant or debris in well
#OD >1 indicates >90% of light blocked, which does not track with images of final day of experiment
treatfinal$OD420[treatfinal$OD420>1] <- NA
treatfinal$OD600[treatfinal$OD600>1] <- NA

#normality tests show optical density data should be logged in models, but plant data is not improved by taking the log
shapiro.test(treatfinal$sqmm); shapiro.test(log(treatfinal$sqmm[-which(treatfinal$sqmm==0)]))
shapiro.test(treatfinal$dsqmm)
shapiro.test(treatfinal$rsqmm)
shapiro.test(treatfinal$greenness); shapiro.test(log(treatfinal$greenness[-which(treatfinal$greenness==0)]))
shapiro.test(treatfinal$OD420); shapiro.test(log(treatfinal$OD420))
shapiro.test(treatfinal$OD600); shapiro.test(log(treatfinal$OD600))
treatfinal$lnOD600 <- log(treatfinal$OD600)
treatfinal$lnOD420 <- log(treatfinal$OD420)

#vectors for indicating which and how many source sites are the same
treatfinal$Allsame <- ifelse(treatfinal$Plant==treatfinal$Water & treatfinal$Plant==treatfinal$Microbe,1,0)
treatfinal$PWsame <- ifelse(treatfinal$Plant==treatfinal$Water,1,0)
treatfinal$PMsame <- ifelse(treatfinal$Plant==treatfinal$Microbe,1,0)
treatfinal$MWsame <- ifelse(treatfinal$Microbe==treatfinal$Water,1,0)
treatfinal$Alldiff <- ifelse(rowSums(cbind(treatfinal$PWsame,treatfinal$PMsame,treatfinal$MWsame))==0,1,0)
treatfinal$NumSamePlant <- treatfinal$PWsame + treatfinal$PMsame
treatfinal$NumSameMicr <- treatfinal$MWsame + treatfinal$PMsame
#vectors for indicating which and how many source watersheds are the same
treatfinal$watershedP <- ifelse(treatfinal$Plant%in%c("WDM","DR"),"Pettee","Oyster")
treatfinal$watershedM <- ifelse(treatfinal$Microbe%in%c("WDM","DR"),"Pettee","Oyster")
treatfinal$watershedW <- ifelse(treatfinal$Water%in%c("WDM","DR"),"Pettee","Oyster")
treatfinal$AllWShedSame <- ifelse(treatfinal$watershedP==treatfinal$watershedM & treatfinal$watershedP==treatfinal$watershedW,1,0)
treatfinal$PW_WShedsame <- ifelse(treatfinal$watershedP==treatfinal$watershedW,1,0)
treatfinal$PM_WShedsame <- ifelse(treatfinal$watershedP==treatfinal$watershedM,1,0)
treatfinal$MW_WShedsame <- ifelse(treatfinal$watershedM==treatfinal$watershedW,1,0)
treatfinal$NumSamePlantWShed <- treatfinal$PW_WShedsame + treatfinal$PM_WShedsame
treatfinal$NumSameMicrWShed <- treatfinal$PM_WShedsame + treatfinal$MW_WShedsame
#we can view response measures as different, as they are not tightly correlated
cor(cbind(treatfinal$sqmm, treatfinal$greenness, treatfinal$lnOD0420,treatfinal$lnOD600),use="complete.obs")

#make subset of inoculated only data
inocOnly <- treatfinal[treatfinal$Microbe != "None",]  


#####
#Produce a full plot of all treatment means
#for a supplemental figure

#take means and standard errors
area_allmn 	<- tapply(treatfinal$sqmm, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
area_allse 	<- tapply(treatfinal$sqmm, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),std.error)
green_allmn <- 100*tapply(treatfinal$greenness, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
green_allse <- 100*tapply(treatfinal$greenness, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),std.error)
od420_allmn <- tapply(log(treatfinal$OD420), paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
od420_allse <- tapply(log(treatfinal$OD420), paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),std.error)
od600_allmn <- tapply(log(treatfinal$OD600), paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
od600_allse <- tapply(log(treatfinal$OD600), paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),std.error)


#set colors and plotting symbols
sourcepch <- c(22,24,25,23) #DR, M, UM, WDM, oyster gets triangles, pettee gets squares/diamonds
sourcecols <- c(rgb(1,0,0),rgb(0,0.5,1),rgb(0,0,1),rgb(1,0.5,0))
sourcecols_micr <- c(rgb(1,0,0),rgb(0,0.5,1),rgb(1,1,1),rgb(0,0,1),rgb(1,0.5,0)) #DR M None UM WDM
sitenames <- c("Durham Reservoir","Woodman","Upper Mill","Mill Pond")
# red for pettee watershed, bright for DR, dark for WDM; blue for oyster, bright for UM, dark for M 
# black for None
water_pch <- rep( sourcepch, each =  5*4 )
plant_col <- rep ( rep( sourcecols, each =  5), times = 4)
micr_col <- rep( sourcecols_micr, times = 4*4 )  

#set order
area_order <- order(area_allmn)
green_order <- order(green_allmn)
od420_order <- order(od420_allmn)
od600_order <- order(od600_allmn)

##get points tagged with a number based on how their means sort above.
pointfac <- paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr)
pointloc_A <- as.vector(sapply(names(area_allmn)[area_order], function(z) which(pointfac==z)))
pointloc_G <- as.vector(sapply(names(green_allmn)[green_order], function(z) which(pointfac==z)))
pointloc_O420 <- as.vector(sapply(names(od420_allmn)[od420_order], function(z) which(pointfac==z)))
pointloc_O600 <- as.vector(sapply(names(od600_allmn)[od600_order], function(z) which(pointfac==z)))

#plot
pdf("alltrtmns.pdf", height = 7, width=5)
par(mfrow=c(4,1))
par(oma=c(1,0,1,1))
par(mar=c(0,4,0,0))
plot(area_allmn[area_order] ~ c(1:length(area_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(treatfinal$sqmm) )
	arrows(c(1:length(area_allmn)), (area_allmn - area_allse)[area_order], y1 = (area_allmn + area_allse)[area_order],
		length=0, lwd=1, col=plant_col[area_order])
	points(area_allmn[area_order] ~ c(1:length(area_allmn)), pch = water_pch[area_order], 
		bg = micr_col[area_order] , col = plant_col[area_order])
	points(treatfinal$sqmm[pointloc_A]~rep(c(1:length(area_allmn)),each=12), pch = rep(water_pch[area_order],each=12), bg = rep(micr_col[area_order],each=12) , col = rep(plant_col[area_order],each=12), cex=0.1)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
plot(green_allmn[green_order] ~ c(1:length(green_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim=c(27, 100*max(treatfinal$greenness)) )
	arrows(c(1:length(green_allmn)), (green_allmn - green_allse)[green_order], y1 = (green_allmn + green_allse)[green_order],
		length=0, lwd=1, col=plant_col[green_order])
	points(green_allmn[green_order] ~ c(1:length(green_allmn)), pch = water_pch[green_order], 
		bg = micr_col[green_order] , col = plant_col[green_order])
	points(100*treatfinal$greenness[pointloc_G]~rep(c(1:length(green_allmn)),each=12), pch = rep(water_pch[green_order],each=12), bg = rep(micr_col[green_order],each=12) , col = rep(plant_col[green_order],each=12), cex=0.1)
	mtext("%Greenness", side=2,line=2)
	#seven points below 25% not plotted
plot(od420_allmn[od420_order] ~ c(1:length(od420_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(log(treatfinal$OD420),na.rm=T) )
	arrows(c(1:length(od420_allmn)), (od420_allmn - od420_allse)[od420_order], y1 = (od420_allmn + od420_allse)[od420_order],
		length=0, lwd=1, col=plant_col[od420_order])
	points(od420_allmn[od420_order] ~ c(1:length(od420_allmn)), pch = water_pch[od420_order], 
		bg = micr_col[od420_order] , col = plant_col[od420_order])
	abline(h=mean(unique(log(treatfinal$OD420.BlankAverage))),lty=3)
	points(log(treatfinal$OD420)[pointloc_O420]~rep(c(1:length(od420_allmn)),each=12), pch = rep(water_pch[od420_order],each=12), bg = rep(micr_col[od420_order],each=12) , col = rep(plant_col[od420_order],each=12), cex=0.1)
	mtext("ln(OD 420 nm)", side=2,line=2)
plot(od600_allmn[od600_order] ~ c(1:length(od600_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim= c(min(log(treatfinal$OD600),na.rm=T),0) ) # range(c(od600_allmn-od600_allse,od600_allmn+od600_allse)) )
	arrows(c(1:length(od600_allmn)), (od600_allmn - od600_allse)[od600_order], y1 = (od600_allmn + od600_allse)[od600_order],
		length=0, lwd=1, col=plant_col[od600_order])
	points(od600_allmn[od600_order] ~ c(1:length(od600_allmn)), pch = water_pch[od600_order], 
		bg = micr_col[od600_order] , col = plant_col[od600_order])
	abline(h=mean(unique(log(treatfinal$OD600.BlankAverage))),lty=3)
	points(log(treatfinal$OD600)[pointloc_O600]~rep(c(1:length(od600_allmn)),each=12), pch = rep(water_pch[od600_order],each=12), bg = rep(micr_col[od600_order],each=12) , col = rep(plant_col[od600_order],each=12), cex=0.1)
	mtext("ln(OD 600 nm)", side=2,line=2)
	text(49,-0,"Water source", adj=0)
	legend(47,-0,sitenames, adj = c(0,0.5),pch=c(22,23,24,25),bty="n")
	text(24,-0,"Duckweed source", adj=0)
	legend(22,-0,sitenames,adj = c(0,0.5), lty=1, col=sourcecols[c(1,4,2,3)],bty="n")
	text(-1,-0,"Microbe source", adj=0)
	legend(-2,-0,c(sitenames,"None"),adj = c(0,0.5),fill=sourcecols_micr[c(1,5,2,4,3)],bty="n")
dev.off()


#same figure, but re-color by which sources of treatments are the same
PWallmn <- tapply(treatfinal$PWsame, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
PMallmn <- tapply(treatfinal$PMsame, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
MWallmn <- tapply(treatfinal$MWsame, paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
isinocallmn <- tapply(as.numeric(as.factor(treatfinal$isinoc)), paste(treatfinal$Water, treatfinal$Plant, treatfinal$Micr),mean,na.rm=T)
#2 = yes, 1 = no
samecols <- rgb(0.8-0.6*PWallmn,0.8-0.6*PMallmn,0.8-0.6*MWallmn) 
samecols_allpts <-  rgb(0.8-0.6*treatfinal$PWsame,0.8-0.6*treatfinal$PMsame,0.8-0.6*treatfinal$MWsame)
#produces ----- dark grey: all same, cyano: just PW same, purple:just PM same, yellow:just MWsame, light gray: none same


pdf("alltrtmns_matchcol.pdf", height = 7, width=5)
par(mfrow=c(4,1))
par(oma=c(1,0,1,1))
par(mar=c(0,4,0,0))
plot(area_allmn[area_order] ~ c(1:length(area_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(treatfinal$sqmm) )
	arrows(c(1:length(area_allmn)), (area_allmn - area_allse)[area_order], y1 = (area_allmn + area_allse)[area_order],
		length=0, lwd=1, col=samecols[area_order])
	points(area_allmn[area_order] ~ c(1:length(area_allmn)), pch = 1+23*(isinocallmn[area_order]-1), 
		 bg = samecols[area_order], col=rgb(0,0,0,alpha=0.25))
	points(treatfinal$sqmm[pointloc_A]~rep(c(1:length(area_allmn)),each=12), pch = 1+23*(as.numeric(as.factor(treatfinal$isinoc[pointloc_A]))-1), bg = samecols_allpts[pointloc_A] , col = rgb(0,0,0,alpha=0.25), cex=0.1)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
plot(green_allmn[green_order] ~ c(1:length(green_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim = c(27, 100*max(treatfinal$greenness)) )
	arrows(c(1:length(green_allmn)), (green_allmn - green_allse)[green_order], y1 = (green_allmn + green_allse)[green_order],
		length=0, lwd=1, col=samecols[green_order])
	points(green_allmn[green_order] ~ c(1:length(green_allmn)), pch = 1+23*(isinocallmn[green_order]-1), 
		 bg = samecols[green_order], col=rgb(0,0,0,alpha=0.25))
	points(100*treatfinal$greenness[pointloc_G]~rep(c(1:length(green_allmn)),each=12), pch = 1+23*(as.numeric(as.factor(treatfinal$isinoc[pointloc_G]))-1), bg = samecols_allpts[pointloc_G] , col = rgb(0,0,0,alpha=0.25), cex=0.1)
	mtext("%Greenness", side=2,line=2)
###eliminates 7 points below greenness = 25%
plot(od420_allmn[od420_order] ~ c(1:length(od420_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim= range(log(treatfinal$OD420),na.rm=T) ) #range(c(od420_allmn-od420_allse,od420_allmn+od420_allse)) )
	arrows(c(1:length(od420_allmn)), (od420_allmn - od420_allse)[od420_order], y1 = (od420_allmn + od420_allse)[od420_order],
		length=0, lwd=1, col=samecols[od420_order])
	points(od420_allmn[od420_order] ~ c(1:length(od420_allmn)), pch = 1+23*(isinocallmn[od420_order]-1), 
		 bg = samecols[od420_order], col=rgb(0,0,0,alpha=0.25))
	abline(h=mean(unique(log(treatfinal$OD420.BlankAverage))),lty=3)
	points(log(treatfinal$OD420)[pointloc_O420]~rep(c(1:length(od420_allmn)),each=12), pch = 1+23*(as.numeric(as.factor(treatfinal$isinoc[pointloc_O420]))-1), bg = samecols_allpts[pointloc_O420] , col = rgb(0,0,0,alpha=0.25), cex=0.1)
	mtext("ln(OD 420 nm)", side=2,line=2)
plot(od600_allmn[od600_order] ~ c(1:length(od600_allmn)), pch =NA, xlab="",ylab="",xaxt="n", ylim= c(min(log(treatfinal$OD600),na.rm=T),0) )
	arrows(c(1:length(od600_allmn)), (od600_allmn - od600_allse)[od600_order], y1 = (od600_allmn + od600_allse)[od600_order],
		length=0, lwd=1, col=samecols[od600_order])
	points(od600_allmn[od600_order] ~ c(1:length(od600_allmn)), pch = 1+23*(isinocallmn[od600_order]-1), 
		 bg = samecols[od600_order], col=rgb(0,0,0,alpha=0.25))
	abline(h=mean(unique(log(treatfinal$OD600.BlankAverage))),lty=3)
	points(log(treatfinal$OD600)[pointloc_O600]~rep(c(1:length(od600_allmn)),each=12), pch = 1+23*(as.numeric(as.factor(treatfinal$isinoc[pointloc_O600]))-1), bg = samecols_allpts[pointloc_O600] , col = rgb(0,0,0,alpha=0.25), cex=0.1)
	mtext("ln(OD 600 nm)", side=2,line=2)
	text(-1,-0.1,"Source sites matching", adj=0)
 	legend(-2,-0.1,c("All", "Duckweeds and Water", "Microbes and Water","Duckweeds and Microbes", "None"), adj = c(0,0.5),fill=unique(samecols),bty="n")
 	legend(32,-0.1,c("Uninoculated","Inoculated"),adj = c(0,0.5),pch=c(1,17),bty="n")
dev.off()


######
### tests for fixed effects of plants, water, microbes
### supplemental figure

#estimate fixed effects of plant source site
area_fixp <- MCMCglmm(sqmm~Plant-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(area_fixp$Sol)  # DR ns> WDM, DR > M & UM, M > UM.   so in DR WDM UM M order that is: 
areafixpt <- c("a", "ab", "c", "b") #order of sitenames
green_fixp <- MCMCglmm(greenness~Plant-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(green_fixp$Sol) # all overlap
od420_fixp <- MCMCglmm(lnOD420~Plant-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od420_fixp$Sol) # UM and WDM lower than DR, but same as each other, and M same as all
od420fixpt <- c("a", "b", "b", "ab")
od600_fixp <- MCMCglmm(lnOD600~Plant-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od600_fixp$Sol) #again UM and WDM lower than DR, but same as each other, and M same as all
od600fixpt <- c("a", "b", "b", "ab")

#estimate fixed effects of microbe source site
area_fixm <- MCMCglmm(sqmm~Microbe-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(area_fixm$Sol) # all but None greater than UM, None sig lower than M
areafixmt <- c("ab","ab","c","a","bc")
green_fixm <- MCMCglmm(greenness~Microbe-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(green_fixm$Sol) # none lower than all others, and DR higher than only UM
greenfixmt <- c("a", "ab", "b", "ab", "c")
od420_fixm <- MCMCglmm(lnOD420~Microbe-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od420_fixm$Sol) #none lower than all others, no other differences sig
od420fixmt <- c("a", "a", "a", "a", "b")
od600_fixm <- MCMCglmm(lnOD600~Microbe-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od600_fixm$Sol) # all inoc same and higher than none
od600fixmt <- c("a", "a", "a", "a", "b") 

#estimate fixed effects of water source site
area_fixw <- MCMCglmm(sqmm~Water-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(area_fixw$Sol) # M higher than DR, all else same
areafixwt <- c("b","ab","ab","a")
green_fixw <- MCMCglmm(greenness~Water-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(green_fixw$Sol) # all ovelap
od420_fixw <- MCMCglmm(lnOD420~Water-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od420_fixw$Sol) #WDM lower than all, no other sig diffs
od420fixwt <- c("a","b","a","a")
od600_fixw <- MCMCglmm(lnOD600~Water-1, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
HPDinterval(od600_fixw$Sol) #WDM lower than M and UM, but no other diffs
od600fixwt <- c("ab","b","a","a")

#calculate means and standard errors for all traits, for plant, microbe, and water source sites
area_dmn <- tapply(treatfinal$sqmm, treatfinal$Plant,mean,na.rm=T)
area_dse <- tapply(treatfinal$sqmm, treatfinal$Plant,std.error)
green_dmn <- 100*tapply(treatfinal$greenness, treatfinal$Plant,mean,na.rm=T)
green_dse <- 100*tapply(treatfinal$greenness, treatfinal$Plant,std.error)
od420_dmn <- tapply(log(treatfinal$OD420), treatfinal$Plant,mean,na.rm=T)
od420_dse <- tapply(log(treatfinal$OD420), treatfinal$Plant,std.error)
od600_dmn <- tapply(log(treatfinal$OD600), treatfinal$Plant,mean,na.rm=T)
od600_dse <- tapply(log(treatfinal$OD600), treatfinal$Plant,std.error)

area_mmn <- tapply(treatfinal$sqmm, treatfinal$Microbe,mean,na.rm=T)
area_mse <- tapply(treatfinal$sqmm, treatfinal$Microbe,std.error)
green_mmn <- 100*tapply(treatfinal$greenness, treatfinal$Microbe,mean,na.rm=T)
green_mse <- 100*tapply(treatfinal$greenness, treatfinal$Microbe,std.error)
od420_mmn <- tapply(log(treatfinal$OD420), treatfinal$Microbe,mean,na.rm=T)
od420_mse <- tapply(log(treatfinal$OD420), treatfinal$Microbe,std.error)
od600_mmn <- tapply(log(treatfinal$OD600), treatfinal$Microbe,mean,na.rm=T)
od600_mse <- tapply(log(treatfinal$OD600), treatfinal$Microbe,std.error)

area_wmn <- tapply(treatfinal$sqmm, treatfinal$Water,mean,na.rm=T)
area_wse <- tapply(treatfinal$sqmm, treatfinal$Water,std.error)
green_wmn <- 100*tapply(treatfinal$greenness, treatfinal$Water,mean,na.rm=T)
green_wse <- 100*tapply(treatfinal$greenness, treatfinal$Water,std.error)
od420_wmn <- tapply(log(treatfinal$OD420), treatfinal$Water,mean,na.rm=T)
od420_wse <- tapply(log(treatfinal$OD420), treatfinal$Water,std.error)
od600_wmn <- tapply(log(treatfinal$OD600), treatfinal$Water,mean,na.rm=T)
od600_wse <- tapply(log(treatfinal$OD600), treatfinal$Water,std.error)

setnumplant <- as.factor(treatfinal$Plant)
levels(setnumplant) <- c(1, 4, 3, 2) #DR, M, UM, WDM become 1, 4, 3, and 2
setnummicr <- as.factor(treatfinal$Microbe)
levels(setnummicr) <- c(1, 4, 5, 3, 2) # "DR"   "M"    "None" "UM"   "WDM" become 1, 4, 5, 3, and 2
setnumwat <- as.factor(treatfinal$Water)
levels(setnumwat) <- c(1, 4, 3, 2) #DR, M, UM, WDM become 1, 4, 3, and 2

#set colors
sourcecols_micr2 <- c(rgb(1,0,0),rgb(0,0.5,1),rgb(0,0,0),rgb(0,0,1),rgb(1,0.5,0)) #DR M None UM WDM

#plot
pdf("MainEffectsMeans.pdf", height = 7, width=4)
layout(matrix(1:12,ncol=3,byrow=T),widths=c(4,5,4))
par(oma=c(8,4,2,1))
par(mar=c(0,0,0,0))
plot(area_dmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(treatfinal$sqmm),xlim=c(0.5,4.5),cex=1.5 )#ylim=c(9,13.5)
	points(treatfinal$sqmm~jitter(as.numeric(as.character(setnumplant))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(area_dmn ~ c(1,4,3,2), pch =1,col=sourcecols,cex=1)
	arrows(c(1,4,3,2),area_dmn-area_dse, y1=area_dmn+area_dse, length=0,col=sourcecols,lwd=2)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
	text(1:4,25,areafixpt)
	mtext("Plant Source",side=3,line=0)
plot(area_mmn ~ c(1,4,5,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(treatfinal$sqmm),xlim=c(0.5,5.5) ,cex=1.5 )
	points(treatfinal$sqmm~jitter(as.numeric(as.character(setnummicr))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(area_mmn ~ c(1,4,5,3,2), pch =1,col=sourcecols_micr2 )
	arrows(c(1,4,5,3,2),area_mmn-area_mse, y1=area_mmn+area_mse, length=0,col=sourcecols_micr2,lwd=2)
	text(1:5,25,areafixmt)
	mtext("Microbe Source",side=3,line=0)
plot(area_wmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(treatfinal$sqmm) ,xlim=c(0.5,4.5) ,cex=1.5 )
	points(treatfinal$sqmm~jitter(as.numeric(as.character(setnumwat))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(area_wmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),area_wmn-area_wse, y1=area_wmn+area_wse, length=0,col=sourcecols,lwd=2)
	text(1:4,25,areafixwt)
	mtext("Water Source",side=3,line=0)
plot(green_dmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", ylim=c(25,100*max(treatfinal$greenness)),xlim=c(0.5,4.5) ,cex=1.5 )
	points(100*treatfinal$greenness~jitter(as.numeric(as.character(setnumplant))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(green_dmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),green_dmn-green_dse, y1=green_dmn+green_dse, length=0,col=sourcecols,lwd=2)
	mtext("%Greenness", side=2,line=2)
	text(2.5,55,"n.s.")
	#cut off 7 points below 25% greennness
plot(green_mmn ~ c(1,4,5,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n",ylim=c(25,100*max(treatfinal$greenness)) ,xlim=c(0.5,5.5) ,cex=1.5 )
	points(100*treatfinal$greenness~jitter(as.numeric(as.character(setnummicr))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(green_mmn ~ c(1,4,5,3,2), pch =1,col=sourcecols_micr2 )
	arrows(c(1,4,5,3,2),green_mmn-green_mse, y1=green_mmn+green_mse, length=0,col=sourcecols_micr2,lwd=2)
	text(1:5,55,greenfixmt)
plot(green_wmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=c(25,100*max(treatfinal$greenness)) ,xlim=c(0.5,4.5) ,cex=1.5 )
	points(100*treatfinal$greenness~jitter(as.numeric(as.character(setnumwat))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(green_wmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),green_wmn-green_wse, y1=green_wmn+green_wse, length=0,col=sourcecols,lwd=2)
	text(2.5,55,"n.s.")
plot(od420_dmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(log(treatfinal$OD420),na.rm=T),xlim=c(0.5,4.5) ,cex=1.5 )
	points(log(treatfinal$OD420)~jitter(as.numeric(as.character(setnumplant))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od420_dmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),od420_dmn-od420_dse, y1=od420_dmn+od420_dse, length=0,col=sourcecols,lwd=2)
	mtext("ln(OD 420 nm)", side=2,line=2)
	text(1:4,-0.25,od420fixpt)
plot(od420_mmn ~ c(1,4,5,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(log(treatfinal$OD420),na.rm=T) ,xlim=c(0.5,5.5) ,cex=1.5 )
	points(log(treatfinal$OD420)~jitter(as.numeric(as.character(setnummicr))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od420_mmn ~ c(1,4,5,3,2), pch =1,col=sourcecols_micr2 )
	arrows(c(1,4,5,3,2),od420_mmn-od420_mse, y1=od420_mmn+od420_mse, length=0,col=sourcecols_micr2,lwd=2)
	text(1:5,-0.25,od420fixmt)
plot(od420_wmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(log(treatfinal$OD420),na.rm=T) ,xlim=c(0.5,4.5) ,cex=1.5 )
	points(log(treatfinal$OD420)~jitter(as.numeric(as.character(setnumwat))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od420_wmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),od420_wmn-od420_wse, y1=od420_wmn+od420_wse, length=0,col=sourcecols,lwd=2)
	text(1:4,-0.25,od420fixwt)
plot(od600_dmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", ylim=range(log(treatfinal$OD600),na.rm=T), xlim=c(0.5,4.5) ,cex=1.5 )
	points(log(treatfinal$OD600)~jitter(as.numeric(as.character(setnumplant))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od600_dmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),od600_dmn-od600_dse, y1=od600_dmn+od600_dse, length=0,col=sourcecols,lwd=2)
	mtext("ln(OD 600 nm)", side=2,line=2)
	axis(side=1,at=c(1:4), sitenames,las=2)
	text(1:4,-0.5,od600fixpt)
plot(od600_mmn ~ c(1,4,5,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(log(treatfinal$OD600),na.rm=T), ,xlim=c(0.5,5.5) ,cex=1.5 )
	points(log(treatfinal$OD600)~jitter(as.numeric(as.character(setnummicr))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od600_mmn ~ c(1,4,5,3,2), pch =1,col=sourcecols_micr2 )
	arrows(c(1,4,5,3,2),od600_mmn-od600_mse, y1=od600_mmn+od600_mse, length=0,col=sourcecols_micr2,lwd=2)
	axis(side=1,at=c(1:5), c(sitenames,"None"),las=2)
	text(1:5,-0.50,od600fixmt)
plot(od600_wmn ~ c(1,4,3,2), pch =NA, xlab="",ylab="",xaxt="n", yaxt="n", ylim=range(log(treatfinal$OD600),na.rm=T),xlim=c(0.5,4.5) ,cex=1.5 )
	points(log(treatfinal$OD600)~jitter(as.numeric(as.character(setnumwat))),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.15))
	points(od600_wmn ~ c(1,4,3,2), pch =1,col=sourcecols )
	arrows(c(1,4,3,2),od600_wmn-od600_wse, y1=od600_wmn+od600_wse, length=0,col=sourcecols,lwd=2)
	axis(side=1,at=c(1:4), sitenames,las=2)
	text(1:4,-0.5,od600fixwt)
dev.off()


#####
###Test effects of inoculation on response variables, create figure
### supplemental figure

#models --- previous analysis shows we should have plant and water as random effects
Inocarea <- MCMCglmm(sqmm~isinoc, random = ~ Plant + Water, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)# 
Inocgreen <- MCMCglmm(greenness~isinoc, random = ~ Plant + Water, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)# 
InocOD600 <- MCMCglmm(lnOD600~isinoc, random = ~ Plant + Water, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)# 
InocOD420 <- MCMCglmm(lnOD420~isinoc, random = ~ Plant + Water, data=treatfinal,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)# 
#means and standard errors
inocareamn 	<- tapply(treatfinal$sqmm, treatfinal$isinoc,mean,na.rm=T)
inocarease 	<- tapply(treatfinal$sqmm, treatfinal$isinoc,std.error)
inocgreenmn <- tapply(treatfinal$greenness, treatfinal$isinoc,mean,na.rm=T)
inocgreense <- tapply(treatfinal$greenness, treatfinal$isinoc,std.error)
inocOD420mn <- tapply(log(treatfinal$OD420), treatfinal$isinoc,mean,na.rm=T)
inocOD420se <- tapply(log(treatfinal$OD420), treatfinal$isinoc,std.error)
inocOD600mn <- tapply(log(treatfinal$OD600), treatfinal$isinoc,mean,na.rm=T)
inocOD600se <- tapply(log(treatfinal$OD600), treatfinal$isinoc,std.error)
#plot figure
pdf("Inoculation.pdf",height=6,width=2)
layout(matrix(1:4,ncol=1))
par(mar=c(0,0,0,0))
par(oma=c(9,4,1,1))
plot(inocareamn~c(1:2),pch=1, cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,2.75),ylim= range(treatfinal$sqmm)) #bg=rgb(0.25,0.25,0.25),#c(10.25,12))#range(c(inocareamn-inocarease, inocareamn+inocarease)))
	points(treatfinal$sqmm~jitter(as.numeric(as.factor(treatfinal$isinoc))),cex=0.25,col=rgb(0,0,0,alpha=0.05))
	arrows(c(1:2),inocareamn-inocarease,y1=inocareamn+inocarease,length=0,lwd=2)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
	text(1.5,25,"p < 0.1")
plot(100*inocgreenmn~c(1:2),pch=1,  cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,2.75),ylim=c(25,100*max(treatfinal$greenness)) ) #100*c(0.405,0.455)) #range(c(inocgreenmn-inocgreense, inocgreenmn+inocgreense)))
	points(100*treatfinal$greenness~jitter(as.numeric(as.factor(treatfinal$isinoc))),cex=0.25,col=rgb(0,0,0,alpha=0.05))
	arrows(c(1:2),100*(inocgreenmn-inocgreense),y1=100*(inocgreenmn+inocgreense),length=0,lwd=2)
	mtext("%Greenness", side=2,line=2)
	text(1.5,55,"*",cex=2)
	#7 points below 25% are not plotted
plot(inocOD420mn~c(1:2),pch=1,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,2.75),ylim= range(log(treatfinal$OD420),na.rm=T)) #c(0.05,0.13))#range(c(inocOD420mn-inocOD420se, inocOD420mn+inocOD420se)))
	points(log(treatfinal$OD420)~jitter(as.numeric(as.factor(treatfinal$isinoc))),cex=0.25,col=rgb(0,0,0,alpha=0.05))
	arrows(c(1:2),inocOD420mn-inocOD420se,y1=inocOD420mn+inocOD420se,length=0,lwd=2)
	mtext("ln(OD 420 nm)", side=2,line=2)
	text(1.5,-0.25,"*",cex=2)
	abline(h=mean(log(unique(treatfinal$OD420.BlankAverage))),lty=3) #even uninoculated microcosms > blanks; host cells may leave pigment in water
plot(inocOD600mn~c(1:2),pch=1,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,2.75),ylim= range(log(treatfinal$OD600),na.rm=T))#c(0.04,0.09))#range(c(inocOD600mn-inocOD600se, inocOD600mn+inocOD600se)))
	points(log(treatfinal$OD600)~jitter(as.numeric(as.factor(treatfinal$isinoc))),cex=0.25,col=rgb(0,0,0,alpha=0.05))
	arrows(c(1:2),inocOD600mn-inocOD600se,y1=inocOD600mn+inocOD600se,length=0,lwd=2)
	mtext("ln(OD 600 nm)", side=2,line=2)
	axis(side=1,at=c(1:2),labels=c("Uninoculated","Inoculated"),las=2,cex=2)
	text(1.5,-0.5,"*",cex=2)
	abline(h=mean(log(unique(treatfinal$OD600.BlankAverage))),lty=3)
dev.off()



#####
###Test for signatures of local adaptation in fitness and trait metrics
#using match to plant for plant growth/traits, match to microbe for microbe density
#use inoculated data when microbe effects included, uninoculated data for only test of plant adaptation to water source
#previous analyses (as well as blanquart et al, other references in main text) shows we should have plant, microbe, and water as random effects
### main text figure

# primary analysis for effects of treatment matching focal (plant or microbe) on response variables
LAarea <- MCMCglmm(sqmm~NumSamePlant, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)#marg pos
LAgreen <- MCMCglmm(greenness~NumSamePlant, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)#n
LAod420p <- MCMCglmm(lnOD420~NumSamePlant, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T) # n
LAod600p <- MCMCglmm(lnOD600~NumSamePlant, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T) # n
LAaream <- MCMCglmm(sqmm~NumSameMicr, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)#marg pos
LAgreenm <- MCMCglmm(greenness~NumSameMicr, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)#n
LAod420 <- MCMCglmm(lnOD420~NumSameMicr, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T) # n
LAod600 <- MCMCglmm(lnOD600~NumSameMicr, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T) # n

#means and standard errors
LAareamn <- tapply(inocOnly$sqmm, inocOnly$NumSamePlant,mean,na.rm=T)
LAarease <- tapply(inocOnly$sqmm, inocOnly$NumSamePlant,std.error)
LAgreenmn <- tapply(inocOnly$greenness, inocOnly$NumSamePlant,mean,na.rm=T)
LAgreense <- tapply(inocOnly$greenness, inocOnly$NumSamePlant,std.error)
LAOD420mnp <- tapply(log(inocOnly$OD420), inocOnly$NumSamePlant,mean,na.rm=T)
LAOD420sep <- tapply(log(inocOnly$OD420), inocOnly$NumSamePlant,std.error)
LAOD600mnp <- tapply(log(inocOnly$OD600), inocOnly$NumSamePlant,mean,na.rm=T)
LAOD600sep <- tapply(log(inocOnly$OD600), inocOnly$NumSamePlant,std.error)

LAareamnm <- tapply(inocOnly$sqmm, inocOnly$NumSameMicr,mean,na.rm=T)
LAareasem <- tapply(inocOnly$sqmm, inocOnly$NumSameMicr,std.error)
LAgreenmnm <- tapply(inocOnly$greenness, inocOnly$NumSameMicr,mean,na.rm=T)
LAgreensem <- tapply(inocOnly$greenness, inocOnly$NumSameMicr,std.error)
LAOD420mn <- tapply(log(inocOnly$OD420), inocOnly$NumSameMicr,mean,na.rm=T)
LAOD420se <- tapply(log(inocOnly$OD420), inocOnly$NumSameMicr,std.error)
LAOD600mn <- tapply(log(inocOnly$OD600), inocOnly$NumSameMicr,mean,na.rm=T)
LAOD600se <- tapply(log(inocOnly$OD600), inocOnly$NumSameMicr,std.error)

#plot figure
pdf("LAtests_full.pdf",height=7,width=6)
layout(matrix(1:16,ncol=4, byrow=F))
par(mar=c(0,2.5,0,0))
par(oma=c(7,2,1,1))
plot(LAareamn~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= range(inocOnly$sqmm))#c(10.25,13.5))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$sqmm~jitter(inocOnly$NumSamePlant+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAareamn-LAarease,y1=LAareamn+LAarease,length=0,lwd=2)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
	text(2,26,"p < 0.1")
plot(100*LAgreenmn~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(25,100*max(inocOnly$greenness)+5) )#c(0.433,0.467))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(100*inocOnly$greenness~jitter(inocOnly$NumSamePlant+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),100*(LAgreenmn-LAgreense),y1=100*(LAgreenmn+LAgreense),length=0,lwd=2)
	mtext("%Greenness", side=2,line=2)
	text(2,60,"n.s.")
	#cut off 7 points with greenness less than 25%
plot(LAOD420mnp~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= range(log(inocOnly$OD420),na.rm=T) ) #c(0.08,0.15))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(log(inocOnly$OD420)~jitter(inocOnly$NumSamePlant+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAOD420mnp-LAOD420sep,y1=LAOD420mnp+LAOD420sep,length=0,lwd=2)
	mtext("ln(OD 420 nm)", side=2,line=2)
	text(2,-0.3,"p < 0.1")
plot(LAOD600mnp~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= range(log(inocOnly$OD600),na.rm=T)) #c(0.06,0.11))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(log(inocOnly$OD600)~jitter(inocOnly$NumSamePlant+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAOD600mnp-LAOD600sep,y1=LAOD600mnp+LAOD600sep,length=0,lwd=2)
	text(2,-0.4,"p < 0.1")
	mtext("ln(OD 600 nm)", side=2,line=2)
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	mtext("Treatments match plants",side=1,line=5,at=4)
#
plot(LAareamn~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(10.25,13.5))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAareamn-LAarease,y1=LAareamn+LAarease,length=0,lwd=2)
	text(2,13.35,"p < 0.1")
plot(100*LAgreenmn~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= 100*c(0.433,0.467))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),100*(LAgreenmn-LAgreense),y1=100*(LAgreenmn+LAgreense),length=0,lwd=2)
	text(2,46.5,"n.s.")
	#cut off 7 points with greenness less than 25%
plot(LAOD420mnp~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(-2.53,-2.17))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAOD420mnp-LAOD420sep,y1=LAOD420mnp+LAOD420sep,length=0,lwd=2)
	text(2,-2.19,"p < 0.1")
plot(LAOD600mnp~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(-2.8,-2.57))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAOD600mnp-LAOD600sep,y1=LAOD600mnp+LAOD600sep,length=0,lwd=2)
	text(2,-2.57,"p < 0.1")
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	#
plot(LAareamnm~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),
		ylim= range(inocOnly$sqmm))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$sqmm~jitter(inocOnly$NumSameMicr+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAareamnm-LAareasem,y1=LAareamnm+LAareasem,length=0,lwd=2)
	text(2,26,"n.s.")
plot(100*LAgreenmnm~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75), ylim= c(25,100*max(inocOnly$greenness)+5) )
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(100*inocOnly$greenness~jitter(inocOnly$NumSameMicr+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),100*(LAgreenmnm-LAgreensem),y1=100*(LAgreenmnm+LAgreensem),length=0,lwd=2)
	text(2,60,"n.s.")
	#cut off 7 points with greenness less than 25%
plot(LAOD420mn~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= range(log(inocOnly$OD420),na.rm=T))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(log(inocOnly$OD420)~jitter(inocOnly$NumSameMicr+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAOD420mn-LAOD420se,y1=LAOD420mn+LAOD420se,length=0,lwd=2)
	text(2,-0.3,"n.s.")
plot(LAOD600mn~c(1:3),pch=1,cex=1,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= range(log(inocOnly$OD600),na.rm=T))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(log(inocOnly$OD600)~jitter(inocOnly$NumSameMicr+1),cex=0.15,col=rgb(0,0,0,alpha=0.15))
	arrows(c(1:3),LAOD600mn-LAOD600se,y1=LAOD600mn+LAOD600se,length=0,lwd=2)
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	text(2,-0.4,"n.s.")
	mtext("Treatments match microbes",side=1,line=5,at=4)
#
plot(LAareamnm~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(10.25,13.5))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAareamnm-LAareasem,y1=LAareamnm+LAareasem,length=0,lwd=2)
	text(2,13.35,"n.s.")
plot(100*LAgreenmnm~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= 100*c(0.433,0.467))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),100*(LAgreenmnm-LAgreensem),y1=100*(LAgreenmnm+LAgreensem),length=0,lwd=2)
	text(2,46.5,"n.s.")
	#cut off 7 points with greenness less than 25%
plot(LAOD420mn~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),ylim= c(-2.53,-2.17))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAOD420mn-LAOD420se,y1=LAOD420mn+LAOD420se,length=0,lwd=2)
	text(2,-2.19,"n.s.")
plot(LAOD600mn~c(1:3),pch=16,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.25,3.75),	ylim= c(-2.8,-2.57))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(c(1:3),LAOD600mn-LAOD600se,y1=LAOD600mn+LAOD600se,length=0,lwd=2)
	text(2,-2.57,"n.s.")
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	#
dev.off()


#####
#Test for signatures of local adaptation in fitness and trait metrics
#using match to plant for plant growth/traits, match to microbe for microbe density
#use inoculated data when microbe effects included, uninoculated data for only test of plant adaptation to water source
# This is NEARLY IDENTICAL to the previous section
# EXCEPT that this analysis accounts for whether the source sites are in the same watershed 
### supplemental figure

#models, as in previous section, but with matching watershed term
WLAarea_p <- MCMCglmm(sqmm~NumSamePlant + NumSamePlantWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAgreen_p <- MCMCglmm(greenness~NumSamePlant + NumSamePlantWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAod420_p <- MCMCglmm(lnOD420~NumSamePlant + NumSamePlantWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAod600_p <- MCMCglmm(lnOD600~NumSamePlant + NumSamePlantWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)

WLAarea_m <- MCMCglmm(sqmm~NumSameMicr + NumSameMicrWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAgreen_m <- MCMCglmm(greenness~NumSameMicr + NumSameMicrWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAod420_m <- MCMCglmm(lnOD420~NumSameMicr + NumSameMicrWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)
WLAod600_m <- MCMCglmm(lnOD600~NumSameMicr + NumSameMicrWShed, random = ~ Plant + Microbe + Water, data=inocOnly,nitt=100000,burnin=10000,thin=10,verbose=F,pr=T)

#means and standard errors
WLAareamn <- tapply(inocOnly$sqmm, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),mean,na.rm=T)
WLAarease <- tapply(inocOnly$sqmm, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),std.error)
WLAgreenmn <- tapply(inocOnly$greenness, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),mean,na.rm=T)
WLAgreense <- tapply(inocOnly$greenness, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),std.error)
WLAod420mn <- tapply(inocOnly$lnOD420, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),mean,na.rm=T)
WLAod420se <- tapply(inocOnly$lnOD420, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),std.error)
WLAod600mn <- tapply(inocOnly$lnOD600, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),mean,na.rm=T)
WLAod600se <- tapply(inocOnly$lnOD600, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),std.error)

WLAareamnm <- tapply(inocOnly$sqmm, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),mean,na.rm=T)
WLAareasem <- tapply(inocOnly$sqmm, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),std.error)
WLAgreenmnm <- tapply(inocOnly$greenness, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),mean,na.rm=T)
WLAgreensem <- tapply(inocOnly$greenness, paste(inocOnly$NumSameMicr,inocOnly$NumSameMicrWShed),std.error)
WLAod420mnp <- tapply(inocOnly$lnOD420, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),mean,na.rm=T)
WLAod420sep <- tapply(inocOnly$lnOD420, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),std.error)
WLAod600mnp <- tapply(inocOnly$lnOD600, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),mean,na.rm=T)
WLAod600sep <- tapply(inocOnly$lnOD600, paste(inocOnly$NumSamePlant,inocOnly$NumSamePlantWShed),std.error)

#set plotting locations for points
pwshed_add <- rep(0,times=nrow(inocOnly))
pwshed_add[inocOnly$NumSamePlant==1 & inocOnly$NumSamePlantWShed==1] <- -0.2
pwshed_add[inocOnly$NumSamePlant==1 & inocOnly$NumSamePlantWShed==2] <- 0.2
pwshed_add[inocOnly$NumSamePlant==0 & inocOnly$NumSamePlantWShed==0] <- -0.2
pwshed_add[inocOnly$NumSamePlant==0 & inocOnly$NumSamePlantWShed==2] <- 0.2
plotnumP <- inocOnly$NumSamePlant + 1 + pwshed_add
mwshed_add <- rep(0,times=nrow(inocOnly))
mwshed_add[inocOnly$NumSameMicr==1 & inocOnly$NumSameMicrWShed==1] <- -0.2
mwshed_add[inocOnly$NumSameMicr==1 & inocOnly$NumSameMicrWShed==2] <- 0.2
mwshed_add[inocOnly$NumSameMicr==0 & inocOnly$NumSameMicrWShed==0] <- -0.2
mwshed_add[inocOnly$NumSameMicr==0 & inocOnly$NumSameMicrWShed==2] <- 0.2
plotnumM <- inocOnly$NumSameMicr + 1 + mwshed_add

#plot
pdf("LA_watershed_tests.pdf",height=7,width=6.5)
WLAxaxt <- c(0.8,1,1.2,1.8,2.2,3)
WLAcolv <- c(1,0.5,0,0.5,0,0)
layout(matrix(1:16,ncol=4, byrow=F))
par(mar=c(0,1.5,0,1))
par(oma=c(7,2.5,1,1))
plot(WLAareamn~WLAxaxt,pch=NA,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25), ylim= range(inocOnly$sqmm))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$sqmm~jitter(plotnumP),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAareamn~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAareamn-WLAarease,y1=WLAareamn+WLAarease,length=0,lwd=2)
	mtext(expression("Frond area, mm"^2), side=2,line=2)
	text(2,27,"local adaptation: *")
plot(100*WLAgreenmn~WLAxaxt,pch=NA,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25), ylim= c(25,100*max(inocOnly$greenness)))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(100*inocOnly$greenness~jitter(plotnumP),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(100*WLAgreenmn~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,100*(WLAgreenmn-WLAgreense),y1=100*(WLAgreenmn+WLAgreense),length=0,lwd=2)
	mtext("%Greenness", side=2,line=2)
	text(2,56,"local adaptation: n.s.")
	#7 greenness points below 25% not included
plot(WLAod420mnp~WLAxaxt,pch=NA,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25), ylim= range(inocOnly$lnOD420,na.rm=T))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$lnOD420~jitter(plotnumP),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAod420mnp~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAod420mnp-WLAod420sep,y1=WLAod420mnp+WLAod420sep,length=0,lwd=2)
	mtext("OD 420 nm", side=2,line=2)
	text(2,-0.25,"local adaptation: n.s.")
plot(WLAod600mnp~WLAxaxt,pch=NA,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25), ylim= range(inocOnly$lnOD600,na.rm=T))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$lnOD600~jitter(plotnumP),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAod600mnp~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAod600mnp-WLAod600sep,y1=WLAod600mnp+WLAod600sep,length=0,lwd=2)
	mtext("OD 600 nm", side=2,line=2)
	text(2,-0.4,"local adaptation: n.s.")
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	mtext("Treatments match plants",side=1,line=5,at=4)
##
plot(WLAareamn~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(10.25,13.5),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAareamn-WLAarease,y1=WLAareamn+WLAarease,length=0,lwd=2)
	text(2,13.45,"local adaptation: *")
plot(100*WLAgreenmn~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= 100*c(0.437,0.467),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,100*(WLAgreenmn-WLAgreense),y1=100*(WLAgreenmn+WLAgreense),length=0,lwd=2)
	text(2,46.65,"local adaptation: n.s.")
plot(WLAod420mnp~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(-2.48,-2.12),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAod420mnp-WLAod420sep,y1=WLAod420mnp+WLAod420sep,length=0,lwd=2)
	text(2,-2.13,"local adaptation: n.s.")
plot(WLAod600mnp~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(-2.8,-2.53),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAod600mnp-WLAod600sep,y1=WLAod600mnp+WLAod600sep,length=0,lwd=2)
	text(2,-2.535,"local adaptation: n.s.")
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
##
plot(WLAareamnm~WLAxaxt,pch=NA,ylab="",xlab="",xlim=c(0.5,3.25),xaxt="n", ylim=  range(inocOnly$sqmm))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$sqmm~jitter(plotnumM),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAareamnm~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAareamnm-WLAareasem,y1=WLAareamnm+WLAareasem,length=0,lwd=2)
	text(2,27,"local adaptation: n.s.")
plot(100*WLAgreenmnm~WLAxaxt,pch=NA,ylab="",xlab="",xlim=c(0.5,3.25),xaxt="n",ylim= c(25,100*max(inocOnly$greenness)))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(100*inocOnly$greenness~jitter(plotnumM),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(100*WLAgreenmnm~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,100*(WLAgreenmnm-WLAgreensem),y1=100*(WLAgreenmnm+WLAgreensem),length=0,lwd=2)
	text(2,56,"local adaptation: n.s.")
plot(WLAod420mn~WLAxaxt,pch=NA,ylab="",xlab="",xlim=c(0.5,3.25),xaxt="n", ylim= range(inocOnly$lnOD420,na.rm=T))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$lnOD420~jitter(plotnumM),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAod420mn~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAod420mn-WLAod420se,y1=WLAod420mn+WLAod420se,length=0,lwd=2)
	text(2,-0.25,"local adaptation: n.s.")
plot(WLAod600mn~WLAxaxt,pch=NA,ylab="",xlab="",xlim=c(0.5,3.25),xaxt="n",ylim= range(inocOnly$lnOD600,na.rm=T))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	points(inocOnly$lnOD600~jitter(plotnumM),pch=16,cex=0.25,col=rgb(0,0,0,alpha=0.4))
	points(WLAod600mn~WLAxaxt,pch=21,cex=1.5,bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	arrows(WLAxaxt,WLAod600mn-WLAod600se,y1=WLAod600mn+WLAod600se,length=0,lwd=2)
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
	mtext("Treatments match microbes",side=1,line=5,at=4)
	text(2,-0.4,"local adaptation: n.s.")
#
plot(WLAareamnm~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(10.25,13.5),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$sqmm),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAareamnm-WLAareasem,y1=WLAareamnm+WLAareasem,length=0,lwd=2)
	text(2,13.45,"local adaptation: n.s.")
plot(100*WLAgreenmnm~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= 100*c(0.437,0.467),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=100*mean(inocOnly$greenness),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,100*(WLAgreenmnm-WLAgreensem),y1=100*(WLAgreenmnm+WLAgreensem),length=0,lwd=2)
	text(2,46.65,"local adaptation: n.s.")
plot(WLAod420mn~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(-2.48,-2.12),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$lnOD420,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAod420mn-WLAod420se,y1=WLAod420mn+WLAod420se,length=0,lwd=2)
	text(2,-2.13,"local adaptation: n.s.")
plot(WLAod600mn~WLAxaxt,pch=21,cex=2,ylab="",xlab="",xaxt="n",xlim=c(0.5,3.25),
		ylim= c(-2.8,-2.53),bg=rgb(WLAcolv,(1-WLAcolv)/2,WLAcolv,alpha=0.5))
	abline(h=mean(inocOnly$lnOD600,na.rm=T),lty=3,col=rgb(0,0,0,alpha=0.5))
	arrows(WLAxaxt,WLAod600mn-WLAod600se,y1=WLAod600mn+WLAod600se,length=0,lwd=2)
	text(2,-2.535,"local adaptation: n.s.")
	axis(side=1,at=c(1:3),labels=c("Neither", "One", "Both"),las=2)
dev.off()


#####
#Contingent tests for local adaptation
# this section considers experiment wide patterns:
	#for fitness in the plant: whether local vs nonlocal water matters in the absence of microbes
	#for fitness in the plant: whether local vs nonlocal microbes matter only when in plant-local water, or only when in nonlocal water
	#for fitness in the plant: whether local vs nonlocal water matters only when in no microbes, in microbes from the plant's local site, or in microbes from nonlocal sites
	#for total microbial cells: whether local vs nonlocal plants matter only when in water from the microbes' local site, or water from a nonlocal site
	#for total microbial cells: whether local vs nonlocal water matters only when with plants from the microbes' local site, or plants from a nonlocal site
	#for both fitness in the plant and total microbial cells: whether all treatments matching is different from when all treatments are mismatching

summary(MCMCglmm(sqmm~PWsame, random=~ Plant + Water, data=treatfinal[treatfinal$isinoc=="no",],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#n.s.
summary(MCMCglmm(sqmm~PWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(sqmm~PWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(sqmm~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(sqmm~PMsame, random=~ Plant + Microbe, data=inocOnly[inocOnly$PWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
#all n.s.
summary(MCMCglmm(sqmm~Allsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$Alldiff==1 | inocOnly$Allsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
#also n.s.
#
summary(MCMCglmm(greenness~PWsame, random=~ Plant + Water, data=treatfinal[treatfinal$isinoc=="no",],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#n.s.
summary(MCMCglmm(greenness~PWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(greenness~PWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(greenness~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(greenness~PMsame, random=~ Plant + Microbe, data=inocOnly[inocOnly$PWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) #marg + (may sometimes above the marginal threshold, no random seed set)
#all but one n.s.
summary(MCMCglmm(greenness~Allsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$Alldiff==1 | inocOnly$Allsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
#also n.s.
#
# uninoculated context does not make sense for OD response
summary(MCMCglmm(lnOD420~PWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))#marg - (may sometimes above the marginal threshold, no random seed set)
summary(MCMCglmm(lnOD420~PWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD420~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD420~PMsame, random=~ Plant + Microbe, data=inocOnly[inocOnly$PWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all but one n.s.
summary(MCMCglmm(lnOD420~Allsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$Alldiff==1 | inocOnly$Allsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
#also n.s.
#
# uninoculated context does not make sense for OD response
summary(MCMCglmm(lnOD600~PWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) #marg - (may sometimes above the marginal threshold, no random seed set)
summary(MCMCglmm(lnOD600~PWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD600~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD600~PMsame, random=~ Plant + Microbe, data=inocOnly[inocOnly$PWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all but one n.s.
summary(MCMCglmm(lnOD600~Allsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$Alldiff==1 | inocOnly$Allsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100))
#also n.s.

##stepping through match to the microbe
summary(MCMCglmm(sqmm~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$MWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(sqmm~PMsame, random=~ Plant + Microbe , data=inocOnly[inocOnly$MWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(sqmm~MWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(sqmm~MWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all n.s.
#all match and all not match already covered
#
summary(MCMCglmm(greenness~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$MWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) #marg + (may sometimes above the marginal threshold, no random seed set)
summary(MCMCglmm(greenness~PMsame, random=~ Plant + Microbe , data=inocOnly[inocOnly$MWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(greenness~MWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(greenness~MWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all but one n.s.
#all match and all not match already covered
#
summary(MCMCglmm(lnOD420~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$MWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD420~PMsame, random=~ Plant + Microbe , data=inocOnly[inocOnly$MWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(lnOD420~MWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(lnOD420~MWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all n.s.
#all match and all not match already covered
summary(MCMCglmm(lnOD600~PMsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$MWsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100))
summary(MCMCglmm(lnOD600~PMsame, random=~ Plant + Microbe , data=inocOnly[inocOnly$MWsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(lnOD600~MWsame, random=~ Plant + Microbe + Water, data=inocOnly[inocOnly$PMsame==0,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
summary(MCMCglmm(lnOD600~MWsame, random=~ Plant + Water, data=inocOnly[inocOnly$PMsame==1,],verbose=F,nitt=100000,burnin=10000,thin=100)) 
#all n.s.
#all match and all not match already covered


#####
###Correlations among response variables, etc.
###Reported rho values
###main text figure

traitmat <- inocOnly[,colnames(inocOnly)%in%c("sqmm","greenness","lnOD420","lnOD600")]
matrix( sapply(1:ncol(traitmat), function(z) sapply(1:ncol(traitmat), function(x) cor(traitmat[,z],traitmat[,x], use="complete.obs"))),ncol=4,nrow=4,byrow=T )
#OD600 and 420 nearly perfectly correlated.
rhodat <- data.frame( rhoOG=c(),rhoOS = c(), rhoGS = c(),
						microbe = c(), plant = c(), water=c())
#get correlations for the traits of 12 replicates of every treatment (for the 64 inoculated treatments = 64 correlations)
for(m in sort(unique(inocOnly$Microbe))){
	for(p in sort(unique(inocOnly$Plant))){
		for(w in sort(unique(inocOnly$Water))){
			datsub <- traitmat[inocOnly$Water==w & inocOnly$Plant == p & inocOnly$Microbe==m,]
			rhodat <- rbind(rhodat, 
					data.frame(rhoOG=cor(datsub$lnOD600,datsub$greenness,use="complete.obs"),
						rhoOS=cor(datsub$lnOD600,datsub$sqmm,use="complete.obs"),
						rhoGS=cor(datsub$greenness,datsub$sqmm,use="complete.obs"),
						microbe = m, plant = p, water=w))
		}
	}

}

siglabs <- array(NA, dim=c(4,3,3))

#check if correlations vary across the 64 treatments (n=64)
green600M <- MCMCglmm(rhoOG~microbe-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#rho higher at 95% CI in microbes DR and UM than microbes M, rho not different from 0 in microbes M
	#at 90% CI, rho in WDM higher than in M, but still not different than DR or UM
	siglabs[,1,1] <- c("a*","b","a*","ab*")	
green600P <- MCMCglmm(rhoOG~plant-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, but rho not different from 0 for plants WDM
	siglabs[,2,1] <- c("*","*","*","")	
green600W <- MCMCglmm(rhoOG~water-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, but rho not different from 0 in water UM, and only marginally pos in water WDM, sig pos in others
	siglabs[,3,1] <- c("*","*","","")	

sqmm600M <- MCMCglmm(rhoOS~microbe-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, but rho not different from 0 in microbes DR, but sig pos in all others
	siglabs[,1,2] <- c("","*","*","*")	
sqmm600P <- MCMCglmm(rhoOS~plant-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, but rho not different from 0 for plants M, but sig pos in all others
	siglabs[,2,2] <- c("*","","*","*")	
sqmm600W <- MCMCglmm(rhoOS~water-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences at 95% CI, but rho not different from 0 in water UM, but sig pos in all others
	siglabs[,3,2] <- c("*","*","","*")	

sqmmgreenM <- MCMCglmm(rhoGS~microbe-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, and all positive rho
	siglabs[,1,3] <- c("*","*","*","*")	
sqmmgreenP <- MCMCglmm(rhoGS~plant-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
	#no differences, and all positive rho
	siglabs[,2,3] <- c("*","*","*","*")	
sqmmgreenW <- MCMCglmm(rhoGS~water-1,data=rhodat, verbose=F, nitt=100000,burnin=10000,thin=100)
 #nothing sig diff from each other at 95% CI level, and all pos rho; but at 90%, rho in water WDM less positive than in water DR
	siglabs[,3,3] <- c("*","*","*","*")	

# calculate means and standard errors 
rhomns <- array(NA,dim=c(4,3,3)) #empty matrix for holding means
rhoses <- array(NA,dim=c(4,3,3)) #empty matrix for holding standard errors
rhonames <- c("rhoOG","rhoOS","rhoGS")
varnames <- c("microbe","plant","water")
for(i in 1:length(rhonames)){  #fill matrices
	for(v in 1:length(varnames)){
		rhomns[,v,i] <- tapply(rhodat[,rhonames[i]], rhodat[,varnames[v]],mean)
		rhoses[,v,i] <- tapply(rhodat[,rhonames[i]], rhodat[,varnames[v]],std.error)
	}
} 

#set up to plot all data observations
sitetoplot <- data.frame(site = c("DR","M","UM","WDM"), plotn = c(1,4,3,2))
rhoplot <- data.frame(microbe = sapply(rhodat$microbe,  function(z) sitetoplot$plotn[sitetoplot$site==z]),
	plant = sapply(rhodat$plant,  function(z) sitetoplot$plotn[sitetoplot$site==z]),
	water = sapply(rhodat$water,  function(z) sitetoplot$plotn[sitetoplot$site==z])
)
sourcecols0 <- c(rgb(1,0,0,alpha=0.4),rgb(1,0.5,0,alpha=0.4),rgb(0,0,1,alpha=0.4),rgb(0,0.5,1,alpha=0.4))
sourcecols_alpha <- c(rgb(1,0,0,alpha=0.4),rgb(0,0.5,1,alpha=0.4),rgb(0,0,1,alpha=0.4),rgb(1,0.5,0,alpha=0.4))

#plot
rholongnames <- c("ln(OD600) & Frond greenness","ln(OD600) & Frond area","Frond greenness & area")
varlongnames <- c("Microbe Source","Plant Source","Water Source")
pdf("ResponseCorrelationsByTrt.pdf",height=5,width=5)
par(mar=c(0,0,0,0))
par(oma=c(8,4,2,2))
par(mfrow=c(3,3))
for(v in 1:length(varnames)){
	for(i in 1:length(rhonames)){
plot(rhomns[,v,i]~c(1,4,3,2),pch=NA,ylim=c(-0.8,1.2),xlim=c(0.5, 4.5),xaxt="n",yaxt="n") #c(-0.2,0.75)
	arrows(c(1,4,3,2),rhomns[,v,i]-rhoses[,v,i],y1=rhomns[,v,i]+rhoses[,v,i],length=0,lwd=2) 
	points(rhomns[,v,i]~c(1,4,3,2),cex=1.75,pch=21,bg=sourcecols_alpha)
	points(rhodat[,rhonames[i]]~jitter(rhoplot[,varnames[v]]),pch=16,cex=0.75,col=rgb(0,0,0,alpha=0.25)) #bg=sourcecols0[rhoplot[,varnames[v]]],
	abline(h=0,lty=3)
	text(c(1,4,3,2),1.1, siglabs[,v,i] )
	if(v==1){mtext(rholongnames[i], side=3,cex=0.7,line=0)}
	if(i==1){axis(side=2, at=c(-0.5,0,0.5,1))}
	if(i==1 & v==2){mtext("Correlation",side=2,line=2)}
	if(i==3){mtext(varlongnames[v],side=4,line=0.25)}
	if(v==3){axis(side=1,labels=sitenames,at=1:4,las=2)}
	}
}  
dev.off()
