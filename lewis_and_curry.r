# Edited version of code from "The implications for climate sensitivity of AR5 forcing and heat uptake estimates"
# Authors: Nicholas Lewis · Judith A. Curry
# Supporting R-code to calculate results. Copyright Nicholas Lewis, September 2014

# Edited to account for temperature-dependent feedback by: Jonah Bloch-Johnson
# Supporting data and orginal code can be found at: https://niclewis.files.wordpress.com/2014/09/lewiscurry-ar5-eb-climdyn-2342-datacode.zip



# R can be downloaded from http://www.r-project.org/ This code was run on the 64 bit Windows version on a machine with 8 GB of memory


# Computation of ECS and TCR using 1859-82 to 1995-2011/2012 periods, AR5 heat uptake estimates
################################################################################################


######################################
# Directory setup and copying of files
######################################

# Set the working directory to whatever directory (folder) is desired, using setwd('directory path')

# Create subdirectory and path structure
dir.create('data')
dir.create('code')
dir.create('plots')
dir.create('GMST')

path.temps= 'GMST'
path.heat= 'data' 
path.home=getwd()

# Manually copy the supplied 'heat_ascii_v6.txt' file to the 'data' subdirectory
# Manually copy the supplied 'AR5_TabAII1.2ForcFin.tab' file to the 'data' subdirectory
# Manually copy the supplied 'AMO_MEI.ext_AR5.F.volc_1851-2012.txt' file to the 'data' subdirectory


##########################
# Input the functions used
##########################

eBudg5= function(n=1e6, pds, base.pds, AR5_Forc=AR5_Forc, AR5_Forc.range=AR5_Forc.range, f.aer2=0, vol.scale=1, BCsnow.scale=1, dQ.base, dQsd.base, intSD.pd=0.08/sqrt(2), intSD.dOHC= 0.045, s.ghg.sc=1, s.aer.sc=1, s.aer.ind=0, s.nonGABC.sc=1, s.nonGABC.ind=0.25, s.nonGABC.fxd=0, s.sol.sc=1, s.sol.ind=0.5, s.sol.fxd=0, s.vol.sc=1, s.vol.ind=1, s.vol.fxd=0, LUC=1, yrs.more=0, dGHG.pa=0.04, s.dT.1750=NA, f2x=3.71, p.points=c(0.05, 0.17, 0.50,0.83,0.95), ECS.lim=c(0,10), TCR.lim=c(0,5), ECS.div=0.01, TCR.div=0.01, R_ECS.lim=c(-1,2.5), R_TCR.lim=c(-0.5,3), R_ECS.div=0.01, R_TCR.div=0.01,
path.temps=path.temps, path.heat=path.heat, path.home=path.home) {

# History and description
#########################
# Code to compute energy budget ECS and TCR estimates and ancilliary information using AR5 forcing and heat uptake estimates and HadCRUT4v2 global surface temperature data

# v5a 14Sep14: minor error affecting reporting of Tot prob by plotBoxCIs3 fixed 
# v5 5Apr14: computes reciprocals of ECS and TCR and their PDFs & makes 2nd aerosol forcing optional
# v4a 2Mar14: enables substituting 1750 zero forcings for any base period, with extra dT uncertainty
# v4 19Jan14: reads new final AR5 heat content data file

# requires plotBoxCIs3 to be loaded or input as text

# n is the number of samples to use.

# pds and base.pds are 2 column matrices with each row consisting of the start and end years of respectively each end and base period to use.
# AR5_Forc is the matrix of forcing data to use, intended to come from Table AII.1.2 in AR5 WG1; rownames should be the years, colnames as in that table: Total, Tot.Anthro, CO2, GHG.Other, O3.Trop, O3.Strat, Aerosol.Tot, LUC, H2O.Strat, BC.Snow, Contrails, Solar, Volcano, Aero.Dir, TotAnth_.78, Total_.78. (The last 3 columns are not in Table AII.1.2: they are not used.)
# AR5_Forc.range is a 3 row matrix with the first column for GHG and the other columns as AR5.forc cosl 5 to 13, with the first two rows the 5% and 95% points of the uncertainty ranges for the corresponding 2011 AR5 best forcing estimates and row 3 scaling factors from RF to ERF uncertainty. Alternatively AR5_Forc.range may have 4 rows, with the first row being median 2011 forcing values to use in place of those in AR5_Forc and other rows all shifted down one.
#if f.aer2 is positive it is treated as a scaling factor to apply to the entire AR5 2011 total aerosol forcing distribution, not just to its central value. If f.aer2 is negative, it is treated as a value to substitute for the -0.9 estimate (median) of the AR5 2011 distribution, whilst maintaining the 5% and 95% points unchanged (as they already allow for the true value being very uncertain). If f.aer2 is zero no alternative aerosol forcing level is used.
# vol.scale and BCsnow.scale respectively scale forcing from volcanic and Black carbon on snow and ice forcings to allow for their ERFs having efficacies that are substantially <1 or >1. If BCsnowscale is set to 3, the AR3 efficacy 5-95% range of 2-4x is used. Default is efficacy of unity.
# dQ.base and dQsd.base are vectors of estimated system heat uptake and its SE (standard deviation) in each of the specified base periods.
# intSD.pd is an estimate of global mean surface temperature internal variability in both the base and end periods, assumed independent. If the estimate used is for the nternal variability in inter-period mean GMST, it should be divided by sqrt(2) to adjust for the estimate being applied twice
# intSD.dOHC is an estimate of internal variability in total system heat uptake duiring both the base and end periods. The same comments apply as for intSD.pd
# s.ghg.sc, s.aer.sc, s.nonGAB.sc, s.sol.sc and s.vol.sc scale the default AR5 uncertainty ranges (without changing their shapes). 
# s.aer.ind, s.nonGAB.ind, s.sol.ind and s.vol.ind set what proportion of the 2011 variance in each of those forcings relates to the level of the cause (e.g., aerosol emissions, ozone emissions, TSI, stratospheric AOD) as opposed to what forcing it results in, but is proportional to the forcing.
# At s.aer.ind= 0.25, the median aer value is ~0.025 more -ve than the specified AR5 median, after adjusting by median*(sqrt(0.75)+sqrt(0.25)-1), even though the 5% and 95% points are very close. Making a positive adjustment of 0.025 to the AR5 median fixes this.
# s.aer.fxd, s.nonGAB.fxd, s.sol.fxd and s.vol.fxd set what proportion of the 2011 variance in each of those forcings is not only independent but also fixed rather than proportionate.
# the rest of the uncertainty is treated as having proportionately the same effect in all years.
# the AR5 Table 8.6 uncertainties are used, with 0.17 thereof added in quadrature where only an RF uncertainty range is given, per Fig.8.16. i.e. Trop O3, Strat O3, Strat H2O ad Land use change 
# LUC: if zero rather than one, land use change 2011 forcing is set to zero and the same again uncertainty added in quadrature, to reflect statement in AR5 that the sign of overall LUC forcing is as likely +ve as -ve. If LUC is < 0, land use change 2011 forcing is set to zero but its uncertainty is left unchanged rather than increased.
# yrs.more and dGHG.pa: if yrs.more > 0 then the end period is extended by that number of years and GHG forcing is treated as rising by dGHG.pa during each extension year, with all other forcings staying fixed. Uncertainty in end-period dQ is scaled down pro rata to its increased length. ECS and TCR median values are kept unchanged by using the increase in the change in total forcing to scale down uncertainty in non-GHG forcing and in heat uptake pro-rata rather than scaling up dQ and dT. Uncertainty in dT is not scaled down since it is relatively small. 
# p.points are 5 one-sided CI intervals at which CDFs are computed for pp.dTFQ (the 3rd dimension of which relates to dT, dF.1, dF.2 and dQ sample values) and for the resulting ECS and TCR estimates.
# s.dT.1750: if a vector of same length as no. of base periods, for those base periods for which s.dT.1750 has a non-zero, non-NA, component forcings are taken as having their 1750 zero start values & temperature as being the same as in the actual base period but with added uncertainty of s.dT.1750
# f2x should usually be left at 3.71 as the AR5 CO2 forcing estimates are based on F_2xCO2=3.71 W/m2
# ECS.lim, ECS.div, etc. specify the PDF limits and resolutions for ECS and TCR and their reciprocals

# Initial setting up
####################
n.pds= nrow(pds)
yr.fin= max(pds)
n.base.pds= nrow(base.pds)
n.pds.comb= n.pds * n.base.pds
n.dT= rnorm(n,0,1); n.dQ= rnorm(n,0,1); n.GHG= rnorm(n,0,1); n.nonGABC= rnorm(n,0,1); n.Aer= rnorm(n,0,1); n.sol= rnorm(n,0,1); n.vol= rnorm(n,0,1)

# check if some base periods are specified to be treated as 1750 for forcing puproses, with extra uncertainty in the base period temperature
if(length(s.dT.1750)==n.base.pds) { s.dT.1750[is.na(s.dT.1750)]= 0 
} else { if(identical(s.dT.1750,NA)) { s.dT.1750= rep(0, n.base.pds) } else { print('Error: s.dT.1750 must either be the default of NA or a vector of same length as no. of base periods'); break }
}

# Generate the uncertainty distributions
########################################
# take the median 2011 forcing values from AR5.Forc unless given in AR5_Forc.range matrix
if( !is.matrix(AR5_Forc.range) || nrow(AR5_Forc.range) > 4 || nrow(AR5_Forc.range) < 3) {
	print("Error: AR5_Forc.range is not a 3 or 4 row matrix"); break() 
} else {
	if( nrow(AR5_Forc.range)==3 ) {
		AR5_Forc.range= rbind( c( sum(AR5_Forc['2011',c('CO2','GHG.Other')]), AR5_Forc['2011',colnames(AR5_Forc.range)[-1]] ), AR5_Forc.range )
		colnames(AR5_Forc.range)[1]= 'GHG'
	} 
}

# First fit shifted lognormal distributions to the 3 skewed forcing uncertainties and sample them
q= AR5_Forc.range[1:3,'BC.Snow']
s= (q[3]-q[2])/(q[1]-q[2])
logsd= abs(log(s-1) / 1.644854)
a= q[1] - (q[3]-q[1])/(s-2) 	
logmu= log((q[3]-q[1])/abs(s-2))
type= sign(s-2)
shift.log= c(logmu, logsd, a, type)
f.BCsnow.2011n= type * rlnorm(n, logmu, logsd) + a 
if(BCsnow.scale>1) {f.BCsnow.2011n= f.BCsnow.2011n * rnorm(n, BCsnow.scale, BCsnow.scale/(3*1.645))}
	
q= AR5_Forc.range[1:3,'Contrails']
s= (q[3]-q[2])/(q[1]-q[2])
logsd= abs(log(s-1) / 1.644854)
a= q[1] - (q[3]-q[1])/(s-2) 	
logmu= log((q[3]-q[1])/abs(s-2))
type= sign(s-2)
shift.log= cbind( shift.log, c(logmu, logsd, a, type) )
f.Contrails.2011n= type * rlnorm(n, logmu, logsd) + a 

q= AR5_Forc.range[1:3,'Aerosol.Tot']
# if part of the aerosol forcing uncertainty is independent, adjust median by +0.025 at s.aer.ind= 0.25, scaled by sqrt(min(s.aer.ind, 1-s.aer.ind)/0.25).  
q[1]= q[1] + 0.025 * sqrt(min(s.aer.ind, 1-s.aer.ind)/0.25)
# create adjustment to restore median for change caused by adding two components summing to > 1. For convenience this is all added to the non-independent component, after scaling
f.aer.adj= ( 1 - sqrt(s.aer.ind) - sqrt(1- s.aer.ind) ) * q[1]

s= (q[3]-q[2])/(q[1]-q[2])
logsd= abs(log(s-1) / 1.644854)
a= q[1] - (q[3]-q[1])/(s-2) 	
logmu= log((q[3]-q[1])/abs(s-2))
type= sign(s-2)
shift.log= cbind( shift.log, c(logmu, logsd, a, type) )
rownames(shift.log)= c('logmu', 'logsd', 'shift', 'sign')
colnames(shift.log)= c('BC.Snow','Contrails','Aerosol.Tot') 

f.aer1.2011n= ( ((type*rlnorm(n,logmu,logsd)+a) / AR5_Forc['2011','Aerosol.Tot'] - 1) * s.aer.sc + 1) * sqrt(1 - s.aer.ind) + f.aer.adj / AR5_Forc['2011','Aerosol.Tot']
# if f.aer2 > 0 scale entire normalised Faer sampling distribution by it, else shift 2011 Faer distribution by difference between -0.90 AR5.Forc median and f.aer2 value and renormalise to AR5.Forc Faer timeseries
if(f.aer2 > 0) { f.aer2.2011n= f.aer1.2011n * f.aer2 } else { if(f.aer2<0) { f.aer2.2011n= f.aer1.2011n + f.aer2/ AR5_Forc['2011','Aerosol.Tot'] -1 } }
f.aer1.2011.com.med= median(f.aer1.2011n)

# if part of the aerosol forcing uncertainty is independent, create two independent set of samples
f.aer1.2011ind.base= f.aer1.2011ind.end= f.aer2.2011ind.base= f.aer2.2011ind.end= 0
if(s.aer.ind > 0) {
	f.aer1.2011ind.base= ( ((type*rlnorm(n,logmu,logsd)+a) / AR5_Forc['2011','Aerosol.Tot'] - 1) * s.aer.sc + 1)* sqrt(s.aer.ind)
	if(f.aer2 > 0) { f.aer2.2011ind.base= f.aer1.2011ind.base* f.aer2 } else { if(f.aer2<0) f.aer2.2011ind.base= f.aer1.2011ind.base }
	f.aer1.2011ind.end= ( ((type*rlnorm(n,logmu,logsd)+a) / AR5_Forc['2011','Aerosol.Tot'] - 1) * s.aer.sc + 1) * sqrt(s.aer.ind)
	if(f.aer2 > 0) { f.aer2.2011ind.end= f.aer1.2011ind.end* f.aer2 } else { if(f.aer2<0) f.aer2.2011ind.end= f.aer1.2011ind.end }
}
f.aer1.2011.tot.med= median(f.aer1.2011n + f.aer1.2011ind.end)

# scaled SE in 2011 GHG forcing & f2x (3.71 per AR5 F_CO2) and in 2011 nonGABC anthropogenic forcings
s.ghg= (AR5_Forc.range[3,'GHG'] - AR5_Forc.range[2,'GHG'])/3.29 * s.ghg.sc  
s.f2x= s.ghg/sum(AR5_Forc['2011',c('CO2','GHG.Other')]) # fractional SE in f2x & total GHG forcing
f2x.n= f2x * (1 + n.GHG * s.f2x)
nonGABC.ranges= AR5_Forc.range[3,c('O3.Trop','O3.Strat','LUC','H2O.Strat')] - AR5_Forc.range[2,c('O3.Trop','O3.Strat','LUC','H2O.Strat')]
# double LUC uncertainty if specified to be symmetric about 0, as => unknown part offsetting albedo
if(LUC==0) { nonGABC.ranges['LUC'] = sqrt(2) * nonGABC.ranges['LUC'] }
s.nonGABC= sqrt( sum(nonGABC.ranges^2 * AR5_Forc.range[4,c('O3.Trop','O3.Strat','LUC','H2O.Strat')]^2 ) ) / 3.29 * s.nonGABC.sc

# create separate common, independent & fixed nonGABC anthropogenic, solar and volcanic uncertainties: common scales with forcing and is the same in both periods; independent scales but differs; fixed differs but stays same

AR5_nonGABC.2011= sum( AR5_Forc['2011', c('O3.Trop','O3.Strat','H2O.Strat')] ) + AR5_Forc['2011','LUC'] * max(0, LUC) 
sfxd.nonGABC= s.nonGABC * sqrt(s.nonGABC.fxd)
sind.nonGABC= s.nonGABC * sqrt(1 - s.nonGABC.fxd) * sqrt(s.nonGABC.ind) / AR5_nonGABC.2011
scom.nonGABC= s.nonGABC * sqrt(1 - s.nonGABC.fxd) * sqrt(1 - s.nonGABC.ind) / AR5_nonGABC.2011

s.sol= (AR5_Forc.range[3,'Solar'] - AR5_Forc.range[2,'Solar'])/3.29 * s.sol.sc
sfxd.sol= s.sol * sqrt(s.sol.fxd)
sind.sol= s.sol * sqrt(1 - s.sol.fxd) * sqrt(s.sol.ind) / AR5_Forc['2011','Solar']
scom.sol= s.sol * sqrt(1 - s.sol.fxd) * sqrt(1 - s.sol.ind) / AR5_Forc['2011','Solar']

s.vol= (AR5_Forc.range[3,'Volcano'] - AR5_Forc.range[2,'Volcano'])/3.29 * s.vol.sc
sfxd.vol= s.vol * sqrt(s.vol.fxd)
sind.vol= s.vol * sqrt(1 - s.vol.fxd) * sqrt(s.vol.ind) / AR5_Forc['2011','Volcano']
scom.vol= s.vol * sqrt(1 - s.vol.fxd) * sqrt(1 - s.vol.ind ) / AR5_Forc['2011','Volcano']

# Surface temperature
#####################
# HadCRUT4 file structure: Col 1 is the date; Col 2 is the temperature anomaly time series for the ensemble member (deg C w.r.t. 1961-1990); Col 3 is the 1 sigma uncertainty in the ensemble member's  time series from the combination of measurement and sampling uncertainty and coverage uncertainty. Uncertainties from systematic biases are described by the set of 100 ensemble member time series.

# allocate storage for 100 HadCRUT4 annual realisations in ensemble, 164 years (1850-2013)
temps.mn= array( dim=c(164, 101) )

# input and parse the downloaded HadCRUT4 data
for (i in 1:100) {
	temp= scan(paste(path.temps,'/HadCRUT.4.3.0.0.annual_ns_avg.', i, '.txt', sep=""), quiet=TRUE)
	if(i==1) { temps.sd= matrix(temp, ncol=6, byrow=TRUE)[,c(1,3)] } # SDs same for all in ensemble
	temps.mn[,c(1,i+1)]= matrix(temp, ncol=6, byrow=TRUE)[1:164,1:2]
}
temp.mn.base= temp.sd.base= vector() 
for (i in 1:n.base.pds) {
	yrRows= intersect(which(!temps.mn[,1] < base.pds[i,1]),which(!temps.mn[,1] > base.pds[i,2]))
	temp.mn.base[i]= mean(temps.mn[yrRows,-1])			 
	temp.sd.base[i]= sqrt(sum(diag(temps.sd[yrRows,-1]^2) + cov(t(temps.mn[yrRows,-1]))))/length(yrRows + s.dT.1750[i]^2 ) 
}
temp.mn.end= temp.sd.end= vector()
for (i in 1:n.pds) {
	yrRows= intersect(which(!temps.mn[,1] < pds[i,1]),which(!temps.mn[,1] > pds[i,2]))
	temp.mn.end[i]= mean(temps.mn[yrRows,-1])						
	temp.sd.end[i]= sqrt(sum(diag(temps.sd[yrRows,-1]^2) + cov(t(temps.mn[yrRows,-1]))))/length(yrRows)							
}

# create initial dTsd & dT vectors (all period combinations, vary fastest by base period) then sample
dTmn= dTsd= vector()
for (i in 1:n.pds) {
	dTmn= c(dTmn, temp.mn.end[i] - temp.mn.base)									
	dTsd= c(dTsd, sqrt( temp.sd.base^2 + temp.sd.end[i]^2 + intSD.pd^2 + intSD.pd^2 ) )	
}
dT= t( t(n.dT %o% dTsd) + dTmn )

# Heat uptake (W/m2; raw data is in ZJ)
#############
AR5hc= matrix(scan(paste(path.heat, '/heat_ascii_v6.txt', sep=""), skip=1, quiet=TRUE), ncol=13, byrow=TRUE)
colnames(AR5hc)= scan(paste(path.heat, '/heat_ascii_v6.txt', sep=""), what='character', nlines=1)
dQ.all= AR5hc[,c('year', 'total_h', 'total_u')]
dQ.all= rbind( dQ.all, c(2012.5, 289.71, 16.30) )	# ssq= 217850464
dQ.end= dQsd.end= vector()
for (i in 1:n.pds) {
	yrRows= intersect(which(!dQ.all[,1] < pds[i,1] + 0.4),which(!dQ.all[,1] > pds[i,2] + 0.6))
	dQ.end[i]= ( dQ.all[rev(yrRows)[1],2] - dQ.all[yrRows[1],2] ) / length(yrRows[-1])
	dQsd.end[i]= sqrt( dQ.all[rev(yrRows)[1],3]^2 + dQ.all[yrRows[1],3]^2 ) / ( length(yrRows[-1]) + yrs.more )	# ading extra years to the end period spreads uncertainty in dQ over a longer period
}
ZJ.yrToW.m2 = 1e21/ (365*24*3600 * 510072e9)
dQ.end= dQ.end * ZJ.yrToW.m2
dQsd.end= dQsd.end * ZJ.yrToW.m2
dQ= dQsd= vector()
for (i in 1:n.pds) {
	dQ= c( dQ, dQ.end[i] - dQ.base )
	dQsd= c( dQsd, sqrt(dQsd.end[i]^2 + dQsd.base^2 + 2 * intSD.dOHC^2) )
}
dQ= t( t(n.dQ %o% dQsd) + dQ )

# Forcing (W/m2)
###############
f.base1= f.base2= array(dim=c(n, n.base.pds))
means.base= array(dim=c(7, n.base.pds)) 
medians.base= array(dim=c(10, n.base.pds)) 
rownames(means.base)= c('f.ghg','f.aer','f.bcs', 'f.con', 'f.nonGABC','f.sol','f.vol')
rownames(medians.base)= c('f.ghg','f.aeroAR5.1','f.aeroAR5.2','f.bcs', 'f.con', 'f.nonGABC','f.sol','f.vol', 'F.1','F.2')
if(f.aer2==0) { medians.base= as.matrix(medians.base[-c(3,10),]) }
for (i in 1:n.base.pds) {
	yrRows.base= intersect(which(!1750:2012 < base.pds[i,1]),which(!1750:2012 > base.pds[i,2]))
	if(s.dT.1750[i] > 0) { yrRows.base= c(1,1) }
	f.ghg.base= rowSums(AR5_Forc[yrRows.base, c('CO2','GHG.Other')])
	f.aer.base= AR5_Forc[yrRows.base, 'Aerosol.Tot']
	f.bcs.base= AR5_Forc[yrRows.base, 'BC.Snow']
	f.con.base= AR5_Forc[yrRows.base, 'Contrails']
	f.nonGABC.base= rowSums( AR5_Forc[yrRows.base, c('O3.Trop','O3.Strat','H2O.Strat')] ) + AR5_Forc[yrRows.base, 'LUC'] * max(0, LUC)
	f.sol.base= AR5_Forc[yrRows.base, 'Solar']
	f.vol.base= AR5_Forc[yrRows.base, 'Volcano']
	means.base[,i]= c(mean(f.ghg.base), mean(f.aer.base), mean(f.bcs.base), mean(f.con.base), mean(f.nonGABC.base), mean(f.sol.base), mean(f.vol.base))

	f.ghg.base= mean(f.ghg.base) * (1 + n.GHG * s.f2x)
	f.aeroAR5.base1= mean(f.aer.base) * (f.aer1.2011n + f.aer1.2011ind.base)
	if(!f.aer2==0) { f.aeroAR5.base2= mean(f.aer.base) * (f.aer2.2011n + f.aer2.2011ind.base) }
	f.bcs.base= mean(f.bcs.base) * f.BCsnow.2011n / AR5_Forc['2011','BC.Snow']
	f.con.base= mean(f.con.base) * f.Contrails.2011n / AR5_Forc['2011','Contrails']
	f.nonGABC.base= mean(f.nonGABC.base) *  ( (1 + n.nonGABC * scom.nonGABC) + rnorm(n,0,1) * sind.nonGABC) + rnorm(n,0,1)*sfxd.nonGABC	
	f.sol.base= mean(f.sol.base) * ((1+n.sol*scom.sol)+rnorm(n,0,1)*sind.sol) + rnorm(n,0,1)*sfxd.sol
	f.vol.base= (mean(f.vol.base) * ( (1 + n.vol * scom.vol) + rnorm(n,0,1) * sind.vol) + rnorm(n,0,1)*sfxd.vol ) * vol.scale
	f.base1[,i]= f.ghg.base + f.nonGABC.base + f.aeroAR5.base1 + f.bcs.base + f.con.base + f.sol.base + f.vol.base
	if(!f.aer2==0) { f.base2[,i]= f.ghg.base + f.nonGABC.base + f.aeroAR5.base2 + f.bcs.base + f.con.base + f.sol.base + f.vol.base }
	if(!f.aer2==0) {
		medians.base[,i]= c(median(f.ghg.base), median(f.aeroAR5.base1), median(f.aeroAR5.base2), median(f.bcs.base), median(f.con.base), median(f.nonGABC.base), median(f.sol.base), median(f.vol.base), median(f.base1[,i]), median(f.base2[,i])) 
		} else { 
		medians.base[,i]= c(median(f.ghg.base), median(f.aeroAR5.base1), median(f.bcs.base), median(f.con.base), median(f.nonGABC.base), median(f.sol.base), median(f.vol.base), median(f.base1[,i]))
	}

}

# create the mean end period forcings, sample uncertainty and then deduct the base period samples
#########################################################
dF.1= array( dim=c(n,n.pds.comb) )
if(!f.aer2==0) { dF.2= array( dim=c(n,n.pds.comb) ) }
means.end= array(dim=c(7, n.pds)) 
medians.end= array(dim=c(10, n.pds)) 
if(f.aer2==0) { medians.end= array(dim=c(8, n.pds)) }
rownames(means.end)= rownames(means.base)
rownames(medians.end)= rownames(medians.base)
z.1= z.2= vector()
for (j in 1:n.pds) {
	yrRows.end= intersect(which(!1750:2012 < pds[j,1]),which(!1750:2012 > pds[j,2]))
	f.ghg.end= rowSums(AR5_Forc[yrRows.end, c('CO2','GHG.Other')])
	f.aer.end= AR5_Forc[yrRows.end, 'Aerosol.Tot']
	f.bcs.end= AR5_Forc[yrRows.end, 'BC.Snow']
	f.con.end= AR5_Forc[yrRows.end, 'Contrails']
	# create nonGABC forcing, adjusting LUC weighting if required
	f.nonGABC.end= rowSums( AR5_Forc[yrRows.end, c('O3.Trop','O3.Strat','H2O.Strat')] ) + AR5_Forc[yrRows.end, 'LUC'] * max(0, LUC)
	f.sol.end= AR5_Forc[yrRows.end, 'Solar']
	f.vol.end= AR5_Forc[yrRows.end, 'Volcano']
	means.end[,j]= c(mean(f.ghg.end), mean(f.aer.end), mean(f.bcs.end), mean(f.con.end), mean(f.nonGABC.end), mean(f.sol.end), mean(f.vol.end))

	f.ghg.end= mean(f.ghg.end) * (1 + n.GHG * s.f2x)
	f.aeroAR5.end1= mean(f.aer.end) * (f.aer1.2011n + f.aer1.2011ind.end)
	if(!f.aer2==0) { f.aeroAR5.end2= mean(f.aer.end) * (f.aer2.2011n + f.aer2.2011ind.end) }
	f.bcs.end= mean(f.bcs.end) * f.BCsnow.2011n / AR5_Forc['2011','BC.Snow']
	f.con.end= mean(f.con.end) * f.Contrails.2011n / AR5_Forc['2011','Contrails']
	f.nonGABC.end= mean(f.nonGABC.end) *  ( (1 + n.nonGABC * scom.nonGABC) + rnorm(n,0,1) * sind.nonGABC) + rnorm(n,0,1)*sfxd.nonGABC	
	f.sol.end= mean(f.sol.end) * ((1+n.sol*scom.sol)+rnorm(n,0,1)*sind.sol) + rnorm(n,0,1)*sfxd.sol
	f.vol.end= (mean(f.vol.end) * ( (1 + n.vol * scom.vol) + rnorm(n,0,1) * sind.vol) + rnorm(n,0,1)*sfxd.vol ) * vol.scale

	# create dF.1 and dF.2, the sampled changes in total forcing in the base and alternative case
	cols.now= ((j-1)*n.base.pds+1):(j*n.base.pds)
	dF.1[,cols.now]= (f.ghg.end + f.nonGABC.end + f.aeroAR5.end1 + f.bcs.end + f.con.end + f.sol.end + f.vol.end) - f.base1
	if(!f.aer2==0) { dF.2[,cols.now]= (f.ghg.end + f.nonGABC.end + f.aeroAR5.end2 + f.bcs.end + f.con.end + f.sol.end + f.vol.end) - f.base2 }
	if(!f.aer2==0) {
		medians.end[,j]= c(median(f.ghg.end), median(f.aeroAR5.end1), median(f.aeroAR5.end2), median(f.bcs.end), median(f.con.end), median(f.nonGABC.end), median(f.sol.end), median(f.vol.end), median(f.ghg.end + f.nonGABC.end + f.aeroAR5.end1 + f.bcs.end + f.con.end + f.sol.end + f.vol.end), median(f.ghg.end + f.nonGABC.end + f.aeroAR5.end2 + f.bcs.end + f.con.end + f.sol.end + f.vol.end))
	} else { 
		medians.end[,j]= c(median(f.ghg.end), median(f.aeroAR5.end1), median(f.bcs.end), median(f.con.end), median(f.nonGABC.end), median(f.sol.end), median(f.vol.end), median(f.ghg.end + f.nonGABC.end + f.aeroAR5.end1 + f.bcs.end + f.con.end + f.sol.end + f.vol.end))
	}

	# simulate the effect of projected data for future years if required
	if(yrs.more>0) {
		# treat extra years added to the end period as leaving all non-GHG forcings and their associated uncertainties unchanged, but now applied to a period with higher average forcing. Actually scale down those uncertainties and leave dF, dQ and dT unchanged so that ECS and TCR stay the same.
		dGHG= dGHG.pa * (yrs.more/2 + 0.5) * yrs.more / (length(yrRows.end) + yrs.more)
		z.1[cols.now]= 1 / ( 1 + dGHG/(medians.end['F.1',j] - medians.base[9,]) )
		if(!f.aer2==0) { z.2[cols.now]= 1 / ( 1 + dGHG/(medians.end[10,j] - medians.base[10,]) ) }
  		sd.GHG= mean(f.ghg.end) *  n.GHG * s.f2x # sd in GHG: add to median to offset downscaling
		# scale down dF then restore missing medians (no uncertainty) and uncertainty in end pd GHG
		dF.1[,cols.now]= t( t(dF.1[,cols.now]) * z.1[cols.now]  +  (1-z.1[cols.now]) %o% (medians.end[9,j] + sd.GHG) - (1-z.1[cols.now]) * medians.base['F.1',] )
		if(!f.aer2==0) { dF.2[,cols.now]= t( t(dF.2[,cols.now]) * z.2[cols.now]  + (1-z.2[cols.now]) %o% (medians.end[10,j] + sd.GHG) - (1-z.2[cols.now]) * medians.base[10,] ) }
		# scale down uncertainty in dQ using the (larger) % change in dF.1
		dQ[,cols.now]= t( t(dQ[,cols.now]) * z.1[cols.now]	 )
	}
}

# allocate summary results storage for all base - end period combinations 
pp.dTFQ= array( dim=c(length(p.points), n.pds.comb, ifelse(f.aer2==0,3,4) ) )
ECSpdf.1= array( dim=c((ECS.lim[2]- ECS.lim[1])/ECS.div+2, n.pds.comb) )
if(!f.aer2==0) { ECSpdf.2= array( dim=c((ECS.lim[2]- ECS.lim[1])/ECS.div+2, n.pds.comb) ) }
TCRpdf.1= array( dim=c((TCR.lim[2]-TCR.lim[1])/TCR.div+2, n.pds.comb) )
if(!f.aer2==0) { TCRpdf.2= array( dim=c((TCR.lim[2]-TCR.lim[1])/TCR.div+2, n.pds.comb) ) }
R_ECSpdf.1= array( dim=c((R_ECS.lim[2]-R_ECS.lim[1])/R_ECS.div+2, n.pds.comb) )
if(!f.aer2==0) { R_ECSpdf.2= array( dim=c((R_ECS.lim[2]-R_ECS.lim[1])/R_ECS.div+2, n.pds.comb) ) }
R_TCRpdf.1= array( dim=c((R_TCR.lim[2]-R_TCR.lim[1])/R_TCR.div+2, n.pds.comb) )
if(!f.aer2==0) { R_TCRpdf.2= array( dim=c((R_TCR.lim[2]-R_TCR.lim[1])/R_TCR.div+2, n.pds.comb) ) }

# create names for each combination of periods
periods= vector()
for (i in 0:(n.pds-1)) {
	for(j in 1:n.base.pds) {
		periods[i*n.base.pds+j]= paste(base.pds[j,1],'-',base.pds[j,2],' to ', pds[i+1,1],'-',pds[i+1,2], sep='')
	}
}

# populate data summary statistics arrays
for (i in 1:n.pds.comb) {
	if(!f.aer2==0) { dat= cbind(dT[,i], dF.1[,i], dF.2[,i], dQ[,i]) } else { dat= cbind(dT[,i], dF.1[,i], dQ[,i]) }
	pp.dTFQ[,i,]= apply(dat, 2, quantile, probs=p.points, names=FALSE, type=8) 
}

# compute sampled temperature changes for doubled CO2
dT.2xCO2= f2x.n * dT

##############################
# Edit for feedback temperature dependence begins:


# create ECS and TCR estimates for 1859-1882/1995-2011 combination.
##############################
F_4x = 2 * f2x.n
ΔF = dF.1[,1]
ΔQ = dQ[,1]
ΔT = dT[,1]
F_1 = f.base1[,1]
Q_1 = dQ.base[1]

print(dim(f.base1))
print(length(dQ.base))

# linear
λ = -(ΔF - ΔQ) / ΔT
λ_linear = λ
ΔT_4x = -F_4x / λ
ΔT_4x[is.na(ΔT_4x)] = Inf
write.csv(ΔT_4x[(λ_linear <= 0) & !is.na(λ_linear)], file = "data/lewis_and_curry_linear.csv")

# a_min
a = -.035
λ = -((((ΔF - ΔQ) + a * (ΔT ** 2)) / ΔT) ** 2 + 4 * a * (F_1 - Q_1)) ** .5
λ_linear = -(ΔF - ΔQ) / ΔT
dT_det = (λ ** 2 - 4 * a * F_4x) ** .5
ΔT_4x = (-λ - dT_det) / (2 * a)
ΔT_4x[is.na(ΔT_4x)] = Inf
write.csv(ΔT_4x[(λ_linear <= 0) & !is.na(λ_linear)], file = "data/lewis_and_curry_a_min.csv")

# a_max
a = .058
λ = -((((ΔF - ΔQ) + a * (ΔT ** 2)) / ΔT) ** 2 + 4 * a * (F_1 - Q_1)) ** .5
λ_linear = -(ΔF - ΔQ) / ΔT
dT_det = (λ ** 2 - 4 * a * F_4x) ** .5
ΔT_4x = (-λ - dT_det) / (2 * a)
ΔT_4x[is.na(ΔT_4x)] = Inf
write.csv(ΔT_4x[(λ_linear <= 0) & !is.na(λ_linear)], file = "data/lewis_and_curry_a_max.csv")

# a_dist
a = runif(2000000,-.06,.06)
λ = -((((ΔF - ΔQ) + a * (ΔT ** 2)) / ΔT) ** 2 + 4 * a * (F_1 - Q_1)) ** .5
λ_linear = -(ΔF - ΔQ) / ΔT
dT_det = (λ ** 2 - 4 * a * F_4x) ** .5
ΔT_4x = (-λ - dT_det) / (2 * a)
ΔT_4x[is.na(ΔT_4x)] = Inf
write.csv(ΔT_4x[(λ_linear <= 0) & !is.na(λ_linear)], file = "data/lewis_and_curry_a_dist.csv")

# a_pos_dist
a = runif(2000000,0,.06)

λ = -((((ΔF - ΔQ) + a * (ΔT ** 2)) / ΔT) ** 2 + 4 * a * (F_1 - Q_1)) ** .5
λ_linear = -(ΔF - ΔQ) / ΔT
dT_det = (λ ** 2 - 4 * a * F_4x) ** .5
ΔT_4x = (-λ - dT_det) / (2 * a)
ΔT_4x[is.na(ΔT_4x)] = Inf
write.csv(ΔT_4x[(λ_linear <= 0) & !is.na(λ_linear)], file = "data/lewis_and_curry_a_pos_dist.csv")


# Create output list
####################
out=list()
out$AR5_Forc.range= AR5_Forc.range
out$n= n 
out$pds= pds
out$base.pds= base.pds
out$f.aer2= f.aer2
out$vol.scale= vol.scale
out$BCsnow.scale= BCsnow.scale
out$dQ.base= dQ.base
out$dQsd.base= dQsd.base
out$LUC= LUC
out$intSD.pd= intSD.pd
out$intSD.dOHC= intSD.dOHC
out$s.ghg.sc= s.ghg.sc
out$s.aer.sc= s.aer.sc
out$s.sol.sc= s.sol.sc 
out$s.nonGABC.sc= s.nonGABC.sc
out$s.nonGABC= s.nonGABC 
out$s.aer.ind= s.aer.ind
out$s.nonGABC.ind= s.nonGABC.ind 
out$s.nonGABC.fxd= s.nonGABC.fxd 
out$s.sol.ind= s.sol.ind 
out$s.sol.fxd= s.sol.fxd 
out$s.vol.sc= s.vol.sc 
out$s.vol.ind= s.vol.ind 
out$s.vol.fxd= s.vol.fxd 
out$shift.log= shift.log
out$yrs.more= yrs.more 
out$dGHG.pa= dGHG.pa
out$s.dT.1750= s.dT.1750
out$s.ghg= s.ghg
out$f2x= f2x
out$s.f2x= s.f2x
out$f.aer.adj= f.aer.adj
out$s.sol= s.sol
out$sfxd.sol= sfxd.sol
out$sind.sol= sind.sol
out$scom.sol= scom.sol
out$s.vol= s.vol
out$sfxd.vol= sfxd.vol
out$sind.vol= sind.vol
out$scom.vol= scom.vol
out$f.aer1.2011n= f.aer1.2011n
out$f.aer1.2011ind.base= f.aer1.2011ind.base
out$f.aer1.2011ind.end= f.aer1.2011ind.end
out$f.aer1.2011.com.med= f.aer1.2011.com.med
out$f.aer1.2011.tot.med= f.aer1.2011.tot.med
out$means.base= means.base
out$means.end= means.end
out$medians.base= medians.base
out$medians.end= medians.end
out$f2x.n= f2x.n
out$dT= dT
out$dTmn= dTmn
out$dTsd= dTsd
out$dF.1= dF.1
if(!f.aer2==0) { out$dF.2= dF.2 }
out$dQ= dQ
out$dQsd= dQsd
out$pp.dTFQ= pp.dTFQ
out$ECS.lim= ECS.lim
out$TCR.lim= TCR.lim
out$ECS.div= ECS.div
out$TCR.div= TCR.div
out$R_ECS.lim= R_ECS.lim
out$R_TCR.lim= R_TCR.lim
out$R_ECS.div= R_ECS.div
out$R_TCR.div= R_TCR.div
out$yrs.more= yrs.more
out$dGHG.pa= dGHG.pa
if(yrs.more>0) out$z.1= z.1
if(yrs.more>0 && !f.aer2==0) { out$z.2= z.2 }

out
}

save(eBudg5, file='code/eBudg5.Rd') 


plotBoxCIs3= function(pdfsToPlot, profLikes=NA, divs=NA, lower, upper, boxsize=0.75, col, yOffset=0.05, profLikes.yOffset=0, spacing=0.75, lwd=3, medlwd=lwd, medcol=col, boxfill=NA, boxlty=1, whisklty=1, plot=TRUE, points=c(0.05,0.25,0.5,0.75,0.95), plotPts=NA, centredPDF=TRUE, pkLike=TRUE, dof=NA, horizontal=TRUE, cdf.pts=FALSE) {

# History and description
#########################
# 1Sep14: added boxlty argument
# 9Nov13. Now accepts single pdfsToPlot and profLikes as vectors as well as 1 column matrices
#v3 30Aug13. Now allows vertical bars (horizontal=FALSE) and inputting of CDF point values (cdf.pts=TRUE) rather than pdfs. One column of pdfsToPlot data per box to be plotted
#v2 1Aug13. Now allows a different number of % points ('points') than 5. If >5 points used and plot=TRUE, specify as plotPts which 5 of the points to use for plotting box and whiskers; the 3rd point should be 50%. Total probability in CDFs now included as final column of output (stats), replaced by peak of profile likelihood where using that method; if pkLike=TRUE the peak of the profile likelihood is plotted instead of the central CDF point.
# v1 18Jul13. Improved function to replace plotBoxCIs. Interpolates properly between grid values

# function to tabulate CDF percentage points from PDFs and profile likelihoods and produce boxplots 


# default of centredPDF=TRUE works on basis that the PDF value at each point is for that point, so that the CDF at each point is the sum of all the lower PDF values + 1/2 the current PDF value, multiplied by the grid spacing. If centredPDF=FALSE then CDF at each value is sum of that and all previous PDF values, implying that the PDF values are an average for the bin ending at that parameter value. Also now plots CIs using signed root log profile likelihood ratios if profLikes supplied. If pkLike=TRUE the vertical bar in the box shows the peak of the profile likelihood estimated from a quadratic fit (given in final column of the returned stats); this is more accurate than the 50th percentile of the CDF estimated using the SRLR method.

#Note re 'TotProb'. The CDF at each PDF point is, when centredPDF= TRUE, computed on the basis that only 50% of the probability for that division is realised by the CDF point located at the PDF point, the PDF point value applying symmetrically about that point. TotProb shows how much probability has been included in the CDF up to the upper limit. This will give a misleading answer if the final point is a sweeper up that includes all values above a limit one division below upper.

# if dof is set to a value, a t-distribution with that many degrees of freedom is used for the signed root profile likelihood ratio test instead of a normal distribution. This is non-standard.
# lower and upper are values of the parameter at start and end of pdfsToPlot columns

	n.points= length(points)
	range= upper - lower 
	if(is.vector(pdfsToPlot)) { pdfsToPlot= matrix(pdfsToPlot, ncol=1) }
	cases= ncol(pdfsToPlot)
	if(!identical(profLikes,NA) && is.vector(profLikes)) { profLikes= matrix(profLikes,ncol=1) }
	likeCases= ifelse(identical(profLikes,NA), 0, ncol(profLikes) )
	if(identical(divs,NA)) { divs= range / (nrow(pdfsToPlot)-1) }
	box=boxplot.stats(runif(100,0,1)) # sets up structure required by bxp
	stats=matrix(NA,cases,n.points+1) 
	if( identical(points, c(0.05,0.25,0.5,0.75,0.95)) ) { colnames(stats)= c('5%', '25%',' 50%', '75%', '95%','TotProb') } else { colnames(stats)= c(as.character(points),'TotProb') }
	if(!length(plotPts)==5) { plotPts=1:n.points }

	if(cdf.pts==FALSE) {
		if(!n.points==5 && plot && !length(plotPts)==5) { print('Error: where if plot==TRUE and points>5, must specify (by position) which 5 of those percentage points are to be box-plotted'); break }
	}
	# loop through posterior PDF cases finding where all the specified CDF points are and plotting
	for (j in 1:cases) {
		if(cdf.pts==FALSE) {
			z= ( cumsum(pdfsToPlot[,j]) - 0.5*pdfsToPlot[,j] * centredPDF ) * divs
			for(k in 1:n.points) {
				i= rev(which(z < points[k]))[1]	# last point at which cumprob < this prob. point
				stats[j,k]= ( i + (points[k] - z[i]) / (z[i+1] - z[i]) ) * divs + lower - divs
			}
			stats[j,n.points+1]= max(z)
		} else {
			stats[j,]= c(pdfsToPlot[,j], 1)
		}

		# create and plot the box for the current posterior PDF case
		box$stats=matrix(stats[j,plotPts], ncol=1)
		if( plot==TRUE) { bxp(box, width = NULL, varwidth = FALSE, notch = FALSE, outline = FALSE, names="", plot = TRUE, border = col[j], col = NULL, boxfill=boxfill[j], pars = list(xaxs='i', boxlty=boxlty, boxwex=boxsize, lwd=lwd, staplewex=0.5, outwex=0.5, whisklty=whisklty[j], medlwd=medlwd, medcol=medcol[j]), horizontal=horizontal, add= TRUE, at=boxsize*yOffset+j*boxsize*spacing, axes=FALSE) }
	}

	# loop through the profile likelihood cases, if any, finding all the specified CI points
	if(likeCases>0) { for(j in 1:likeCases) {
		percPts= vector()
		srlr= log(profLikes[,j])
		srlr[is.infinite(srlr)]= -1e12
		srlr.like= sign( 1:length(srlr) - which.max(srlr) ) * sqrt( 2 * (max(srlr) - srlr) ) 
		if(identical(dof,NA)) {srlr.cdf= pnorm(srlr.like)} else { srlr.cdf= pt(srlr.like,dof)}
		for (k in 1:n.points) {
			test= which(srlr.cdf>points[k])
			if(length(test)>0) {
				percPts[k]= test[1]	
				if(percPts[k]==1) {
					percPts[k]= 0
				} else {
					percPts[k]= percPts[k] + (points[k]-srlr.cdf[percPts[k]])/(srlr.cdf[percPts[k]] - srlr.cdf[percPts[k]-1])
				}
			} else {
				percPts[k]= length(srlr) + 1
			}
		}

		# find peak of the likelihood function, assumed quadratic near there and put in final column
		x= which.max(srlr)
		if(x>1 & x< nrow(profLikes) ) {
			X_mat= matrix(c(1,1,1,x-1,x,x+1,(x-1)^2,x^2,(x+1)^2), ncol=3)
			coeffs= solve(X_mat) %*% c(srlr[x-1], srlr[x], srlr[x+1])
			# sub-divide by 20 cells above and below peak and find which is peak, assuming quadratic
			x.cands= seq(x-1, x+1, length= 41)
			y.cands= matrix(c(rep(1,41),x.cands,x.cands^2), ncol=3) %*% coeffs
			percPts[n.points+1]= x + (which.max(y.cands)-21)/20 
		} else {
			percPts[n.points+1]= NA
		}
				
		percPts= (percPts-1) * divs + lower 	# convert into parameter positions into values
		stats=rbind( stats, percPts )			# add the current profLikes case percPts to stats

		# create and plot the box for the current profile likelihood case
		box$stats=matrix(percPts[plotPts], ncol=1)
		if(pkLike) { box$stats[3,1]= percPts[n.points+1] }	# if pkLike=TRUE replace median with peak
		if( plot==TRUE) { bxp(box, width = NULL, varwidth = FALSE, notch = FALSE, outline = FALSE, names="", plot = TRUE, border = col[j+cases], col = NULL, boxfill=boxfill[j], pars = list(xaxs='i', boxlty=1, boxwex=boxsize, lwd=lwd, staplewex=0.5, outwex=0.5, whisklty=whisklty[j+cases], medlwd=medlwd, medcol=medcol[j]), horizontal=horizontal, add=TRUE, at=boxsize*yOffset+(j+cases)*boxsize*spacing+profLikes.yOffset, axes=FALSE) }

} }
stats
}

save(plotBoxCIs3, file='code/plotBoxCIs3.Rd') 


##################
# prepare the data
##################

# download the HadCRUT4 data
############################

# download.file(url='http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.3.0.0.annual_ns_avg_realisations.zip', destfile='GMST/HadCRUT4_ensemble.zip')
# unzip('GMST/HadCRUT4_ensemble.zip', exdir='GMST')

# parse the AR5 forcing data from Table AII.1.2
###############################################

AR5_Forc= as.matrix(read.table(paste(path.home, '/data/AR5_TabAII1.2ForcFin.tsv', sep=""), header=TRUE))
rownames(AR5_Forc)= as.character(AR5_Forc[,1])
AR5_Forc= AR5_Forc[,-1]

# create matrix of forcing uncertainty ranges & scaling factor for ERF to RF uncertainty: AR5 8.SM.1
#the RF to ERF uncertainty addition is 17% of the RF, not 17% of the uncertainty. It is added to the 5-50% and 50-95% ranges
# ozone ERF is given as 0.141 to 0.559 = 0.418 for combined O3.Trop & O3.Strat. Adding their RF ranges in quadrature would give sqrt(0.4^2+0.2^2)= 0.447. So they need to be scaled down, by 0.9347
# from their combined 0.4 5-95% RF range, add 0.17*0.35*2=0.119,^2 => 0.4173 ~agrees to Table 8.SM.5
# for LUC with a 0.2 5-95% RF range, add 0.17*0.15*2=0.051; sqrt(0.2^2 + 0.051^2)= 0.2064 => x 1.032
# for H2O.Strat with a 0.1 5-95% RF range in AR5, add 0.17*0.07*2=0.0238 => 0.10279 = x 1.028
# for BC.Snow the 0.06 RF range in Table 8.SM.5 has already been raised to ERF 0.07 in Table 8.6 
# for solar with a 0.10 5-95% RF range in AR5, add 0.17*0.05*2=0.017 => 0.10143 = x 1.014
# Volcano: 0.07 5-95% RF range in AR5 for 2008-11 RF of -0.11: raise range to 0.16 as discussed in the paper. Then add 0.17*0.11*2=0.0374, which added in quadrature => 0.164 = x 1.025 

AR5_Forc.range= rbind( c(2.26,0.20,-0.15,-1.9,-0.25,0.023,0.02,0.02,0-0.02,-0.205), c(3.40,0.60,0.05,-0.1,-0.05,0.123,0.09,0.15,0.1-0.02,-0.045), c(1,0.9347,0.9347,1,1.032,1.028,1,1,1.014,1.025))
colnames(AR5_Forc.range)= c('GHG', 'O3.Trop', 'O3.Strat', 'Aerosol.Tot', 'LUC', 'H2O.Strat', 'BC.Snow', 'Contrails', 'Solar', 'Volcano')

# specify base and end periods to use
#####################################

base.pds= rbind( c(1859,1882), c(1850,1900), c(1930, 1950) ) 
end.pds= rbind( c(1995, 2011), c(1987, 2011), c(1971, 2011) )

############################
# run the main results cases
############################

eb= eBudg5(n=2e6, pds=end.pds, base.pds=base.pds, AR5_Forc=AR5_Forc, AR5_Forc.range=AR5_Forc.range, f.aer2=0, vol.scale=1, BCsnow.scale=3, dQ.base=c(0.15,0.1,0.2), dQsd.base= c(0.15,0.1,0.2)*0.5, intSD.pd= 0.08/sqrt(2), intSD.dOHC= 0.045, s.ghg.sc=1, s.aer.sc=1, s.aer.ind=0.25, s.nonGABC.sc=1, s.nonGABC.ind=0.5, s.nonGABC.fxd=0, s.sol.sc=sqrt(2), s.sol.ind=0.5, s.sol.fxd=0.5, s.vol.sc=1, s.vol.ind=0.5, s.vol.fxd=1-0.4375^2, LUC=1, p.points=c(0.01, 0.05, 0.17, 0.50,0.83,0.95, 0.99), path.temps=path.temps, path.heat=path.heat, path.home=path.home) 
