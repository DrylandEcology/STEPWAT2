source("StepWatDebug.R")

#This is the Year you want to get the data from
Year <- 1980
Iteration <- 1

#Some tables we will need
rgroups<-dbReadTable(con,"rgroups")
prod_conv <- dbReadTable(con,"sxwprod")#This is sxwprod.in
plotsize <- dbGetQuery(con, "SELECT PlotSize FROM info;")$PlotSize

#some variable for the calculations
totbmass <- 0
bmassg <- numeric(nrow(rgroups))

#This gets the biomass / plotsize and puts it in bmassg and sums them all up
for(g in rgroups$ID) {
	temp <- dbGetQuery(con, paste("SELECT Biomass FROM sxwoutputrgroup WHERE Year=",Year," AND Iteration=",Iteration," AND RGroupID=",g,";",sep=""))$Biomass
	temp <- temp / plotsize
	bmassg[g] <- temp
	totbmass <- totbmass + temp

}
#some variable for the calculations
props <- numeric(12)
props1 <- numeric(12)
props2 <- numeric(12)
props3 <- numeric(12)
props4 <- numeric(12)

cumprop <- numeric(12)
cumprop1 <- numeric(12)
cumprop2 <- numeric(12)
cumprop3 <- numeric(12)
cumprop4 <- numeric(12)

biomass1 <- 0
biomass2 <- 0
biomass3 <- 0
biomass4 <- 0
#The input prod files that are written to SOILWAT
prod_tree <- matrix(data=0, nrow=12, ncol=4, dimnames=list(c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec"),c("Litter","Biomass","PCTLive","LAI_conv")))
prod_shrub <- matrix(data=0, nrow=12, ncol=4, dimnames=list(c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec"),c("Litter","Biomass","PCTLive","LAI_conv")))
prod_grass <- matrix(data=0, nrow=12, ncol=4, dimnames=list(c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec"),c("Litter","Biomass","PCTLive","LAI_conv")))
prod_forb <- matrix(data=0, nrow=12, ncol=4, dimnames=list(c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec"),c("Litter","Biomass","PCTLive","LAI_conv")))

#m is month g is rgourpid
for(m in 1:12) {
	if(totbmass > 0) {
		for(g in rgroups$ID) {
			phen <- dbGetQuery(con,paste("SELECT GrowthPCT FROM sxwphen WHERE RGroupID=",g," AND Month=",m,";",sep=""))$GrowthPCT
			props[m] <- props[m] + phen * bmassg[g]
			if(rgroups[g,"VegProdType"] == 1) {
				props1[m] <- props1[m] + phen * bmassg[g]
			} else if(rgroups[g,"VegProdType"] == 2) {
				props2[m] <- props2[m] + phen * bmassg[g]
			} else if(rgroups[g,"VegProdType"] == 3) {
				props3[m] <- props3[m] + phen * bmassg[g]
			} else if(rgroups[g,"VegProdType"] == 4) {
				props4[m] <- props4[m] + phen * bmassg[g]
			}
		}
		props[m] <- props[m] / totbmass;
		props1[m] <- props1[m] / totbmass;
		props2[m] <- props2[m] / totbmass;
		props3[m] <- props3[m] / totbmass;
		props4[m] <- props4[m] / totbmass;
		
		prod_tree[m,3] <- props1[m]^(.2)
		prod_shrub[m,3] <- props2[m]^(.2)
		prod_grass[m,3] <- props3[m]^(.2)
		prod_forb[m,3] <- props4[m]^(.2)
		
		cumprop[m] <- cumprop[m] + props[m]
		cumprop1[m] <- cumprop1[m] + props1[m]
		cumprop2[m] <- cumprop2[m] + props2[m]
		cumprop3[m] <- cumprop3[m] + props3[m]
		cumprop4[m] <- cumprop4[m] + props4[m]
		
		prod_tree[m,2] <- prod_tree[m,3] * cumprop1[m] * totbmass
		prod_tree[m,1] <- prod_tree[m,2] * prod_conv[m,"LITTER"]
		prod_shrub[m,2] <- prod_shrub[m,3] * cumprop2[m] * totbmass
		prod_shrub[m,1] <- prod_shrub[m,2] * prod_conv[m,"LITTER"]
		prod_grass[m,2] <- prod_grass[m,3] * cumprop3[m] * totbmass
		prod_grass[m,1] <- prod_grass[m,2] * prod_conv[m,"LITTER"]
		prod_forb[m,2] <- prod_forb[m,3] * cumprop4[m] * totbmass
		prod_forb[m,1] <- prod_forb[m,2] * prod_conv[m,"LITTER"]
		
		biomass1 <- biomass1 + prod_tree[m,2]
		biomass2 <- biomass2 + prod_shrub[m,2]
		biomass3 <- biomass3 + prod_grass[m,2]
		biomass4 <- biomass4 + prod_forb[m,2]
	} else {
		prod_tree[m,] <- c(0,0,0,0)
		prod_shrub[m,] <- c(0,0,0,0)
		prod_grass[m,] <- c(0,0,0,0)
		prod_forb[m,] <- c(0,0,0,0)
	}
}
#These are the fractions that are written to soilwat
fractions <- c(0,0,0,0,0)
names(fractions) <- c("Tree","Shrub","Grass","Frob","BareGround")

#Calculation of the fractions
if(totbmass > 0) {
	biomass <- biomass1 + biomass2 + biomass3 + biomass4
	fractions[1] <- biomass1 / biomass
	fractions[2] <- biomass2 / biomass
	fractions[3] <- biomass3 / biomass
	fractions[4] <- biomass4 / biomass
} else {
	fractions[4] <- 1
}

print("Calculated Values")
print("Tree")
print(prod_tree[,1:3])
print("Shrub")
print(prod_shrub[,1:3])
print("Grass")
print(prod_grass[,1:3])
print("Forb")
print(prod_forb[,1:3])
print(fractions)

#This pulls these values from the database. Notice I used the next year. The values might be off because of the positioning of the printDebug function in the main loop before some mort stuff happens.
print("Actual Values")
print(dbGetQuery(con,paste("SELECT Litter, Biomass, PLive FROM sxwinputprod WHERE YEAR=",Year+1," AND Iteration=",Iteration," AND VegProdType=1 ORDER BY Month;",sep="")))
print(dbGetQuery(con,paste("SELECT Litter, Biomass, PLive FROM sxwinputprod WHERE YEAR=",Year+1," AND Iteration=",Iteration," AND VegProdType=2 ORDER BY Month;",sep="")))
print(dbGetQuery(con,paste("SELECT Litter, Biomass, PLive FROM sxwinputprod WHERE YEAR=",Year+1," AND Iteration=",Iteration," AND VegProdType=3 ORDER BY Month;",sep="")))
print(dbGetQuery(con,paste("SELECT Litter, Biomass, PLive FROM sxwinputprod WHERE YEAR=",Year+1," AND Iteration=",Iteration," AND VegProdType=4 ORDER BY Month;",sep="")))
print(dbGetQuery(con,paste("SELECT FracTree,FracShrub,FracGrass,FracForb,FracBareGround FROM sxwinputvars WHERE YEAR=",Year+1," AND Iteration=",Iteration,";",sep="")))


