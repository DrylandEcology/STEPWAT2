library(RSQLite)
#/home/ryan/workspaceLuna/StepWat_SoilWat31/testing/Stepwat Inputs/Output/
sxwdebug <- "sxwdebug.sqlite3"

drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname = sxwdebug)
Tables <- dbListTables(con)

plotMAT <- function(Iteration) {
	sql <- paste("SELECT Year,MAT_C FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep="")
	res <- dbGetQuery(con, sql)
	plot(res, type="l")
}

plotMAP <- function(Iteration) {
	sql <- paste("SELECT Year,MAP_mm FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep="")
	res <- dbGetQuery(con, sql)
	plot(res, type="l")
}

plotAET <- function(Iteration) {
	sql <- paste("SELECT Year,AET_cm FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep="")
	res <- dbGetQuery(con, sql)
	plot(res, type="l")
}

plotAT <- function(Iteration) {
	sql <- paste("SELECT Year,AT_cm FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep="")
	res <- dbGetQuery(con, sql)
	plot(res, type="l")
}

plotAET_AT <- function(Iteration) {
	years <- dbGetQuery(con, paste("SELECT Year FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))$Year
	AET <- dbGetQuery(con, paste("SELECT Year, AET_cm FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))
	AT <- dbGetQuery(con, paste("SELECT Year,AT_cm FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))
	plot(AET, type="o", col="blue",ylim=c(min(AT$AT_cm),max(AET$AET_cm)))
	lines(AT, type="o",pch=22,lty=2,col="red")
}

plotRGroupTotals <- function(Iteration) {
	TotalRelsize <- dbGetQuery(con, paste("SELECT Year, TotalRelsize FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))
	TotalPR <- dbGetQuery(con, paste("SELECT Year,TotalPR FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))
	TotalTransp <- dbGetQuery(con, paste("SELECT Year,TotalTransp FROM sxwoutputvars WHERE Iteration=",Iteration,";",sep=""))
	minValue <- min(c(TotalRelsize$TotalRelsize,mean(TotalPR$TotalPR),TotalTransp$TotalTransp))
	maxValue <- max(c(TotalRelsize$TotalRelsize,mean(TotalPR$TotalPR),TotalTransp$TotalTransp))
	plot(TotalRelsize, type="o", col="blue",ylim=c(minValue,maxValue))
	lines(TotalPR, type="o",pch=22,lty=2,col="red")
	lines(TotalTransp, type="o",pch=22,lty=2,col="green")
}

plotInputFractions <- function(Iteration) {
	fracs<-dbGetQuery(con, paste("SELECT FracGrass,FracShrub,FracTree,FracForb,FracBareGround FROM sxwinputvars WHERE Iteration=",Iteration,";",sep=""))
	par(xpd=T, mar=par()$mar+c(0,0,0,6))
	barplot(t(fracs), main="Fractions", ylab="Fraction",col=c("green","red","orange","blue","black"))
	legend(x=dim(fracs)[1]+25,y=.5, names(fracs), cex=0.8, fill=c("green","red","orange","blue","black"));
	par(mar=c(5, 4, 4, 2) + 0.1)
	
}

plotRGroupBiomass <- function(RGroupID, Iteration) {
	minRsize <- dbGetQuery(con, paste("SELECT MIN(Biomass) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(Biomass) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	data <- dbGetQuery(con, paste("SELECT Biomass FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$Biomass
	plot(x=years,y=data,type="l",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="Biomass",xlab="Years")
	title(main=rgroups$NAME[RGroupID])
}

plotRGroupRealSize <- function(RGroupID, Iteration) {
	minRsize <- dbGetQuery(con, paste("SELECT MIN(Realsize) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(Realsize) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	data <- dbGetQuery(con, paste("SELECT Realsize FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$Realsize
	plot(x=years,y=data,type="l",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="Realsize",xlab="Years")
	title(main=rgroups$NAME[RGroupID])
}


plotRGroupsRealSize <- function(Iteration) {
	colors <- c("green","red","blue","yellow","orange","brown","purple")
	minRsize <- dbGetQuery(con, paste("SELECT MIN(Realsize) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(Realsize) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	plot(1,type="n",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="Realsize",xlab="Years" )
	for(i in rgroups$ID) {
		data <- dbGetQuery(con, paste("SELECT Realsize FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",i,";",sep=""))$Realsize
		lines(y=data,x=years, type="l",lty=2,lwd=2,col=colors[i])
	}
	
	legend(x="topleft", rgroups$NAME, cex=.8,col=colors,lty=rgroups$ID, lwd=2,bty="n")
}

plotRGroupPR <- function(RGroupID, Iteration) {
	minRsize <- dbGetQuery(con, paste("SELECT MIN(PR) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(PR) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	data <- dbGetQuery(con, paste("SELECT PR FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))$PR
	plot(x=years,y=data,type="l",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="PR",xlab="Years")
	title(main=rgroups$NAME[RGroupID])
}

plotRGroupsPR <- function(Iteration) {
	colors <- c("green","red","blue","yellow","orange","brown","purple")
	minRsize <- dbGetQuery(con, paste("SELECT MIN(PR) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(PR) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	plot(1,type="n",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="PR",xlab="Years" )
	for(i in rgroups$ID) {
		data <- dbGetQuery(con, paste("SELECT PR FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",i,";",sep=""))$PR
		lines(y=data,x=years, type="l",lty=2,lwd=2,col=colors[i])
	}
	
	legend(x="topleft", rgroups$NAME, cex=.8,col=colors,lty=rgroups$ID, lwd=2,bty="n")
}

plotRGroupsTransp <- function(Iteration) {
	colors <- c("green","red","blue","yellow","orange","brown","purple")
	minRsize <- dbGetQuery(con, paste("SELECT MIN(Transpiration) AS min FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$min
	maxRsize <- dbGetQuery(con, paste("SELECT MAX(Transpiration) AS max FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$max
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$YEAR
	rgroups<-dbReadTable(con,"rgroups")
	
	plot(1,type="n",ylim=c(minRsize,maxRsize), xlim=c(min(years),max(years)),ylab="PR",xlab="Years" )
	for(i in rgroups$ID) {
		data <- dbGetQuery(con, paste("SELECT Transpiration FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",i,";",sep=""))$Transpiration
		lines(y=data,x=years, type="l",lty=2,lwd=2,col=colors[i])
	}
	
	legend(x="topleft", rgroups$NAME, cex=.8,col=colors,lty=rgroups$ID, lwd=2,bty="n")
}

plotRGroupSummary <- function(RGroupID, Iteration) {
	colors <- c("green","red","blue","yellow","orange","brown","purple")
	years <- dbGetQuery(con, paste("SELECT DISTINCT Year FROM sxwoutputrgroup WHERE Iteration=",Iteration,";",sep=""))$YEAR
	RGroup <- dbGetQuery(con, paste("SELECT Realsize,PR,Transpiration FROM sxwoutputrgroup WHERE Iteration=",Iteration," AND RGroupID=",RGroupID,";",sep=""))
	rgroups<-dbReadTable(con,"rgroups")
	
	plot.ts(RGroup)
	title(rgroups$NAME[RGroupID])
}


#dev.new()
#plotMAT(1)
#dev.new()
#plotMAP(1)
#dev.new()
#plotAET_AT(1)
#dev.new()
#plotRGroupTotals(1)


