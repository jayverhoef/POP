  library(POP)
 	library(rgeos)
  library(rgdal)
	library(sp)
	library(classInt)
 	library(viridis)
  
#	library(Matrix)

	data(POPhexagons)
 	data(AKpoly)
	data(GOA_NavyArea)
 
  
	data(nb_list)
	data(SigiList)
	data(detsig)
	data(Dmat)
	data(W)
	data(MCMC_eff_zmat)
	data(MCMC_eff_beta0)
	data(MCMC_CU_zmat)
	data(MCMC_CU_beta0)
	data(MCMC_EJ_zmat)
	data(MCMC_EJ_beta0)
	data(MCMC_CU_smoo)
	data(MCMC_EJ_smoo)

# get the centroids of the polygons for faster plotting
	trueCentroids = gCentroid(POPhexagons,byid=TRUE)
# extract the data from POPhexagons
	ptDF = POPhexagons@data
#	pts = SpatialPointsDataFrame(ptDF[, c('Long','Lat')], data = ptDF, 
#		proj4string = CRS("+init=epsg:4326")) #latlong projection
	
setwd(paste0('/home/jay/data/2019_papers/POP',
  '/POP_package/POP/inst/doc/figures'))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-study_area
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

	pdf(file = 'Fig-study_area.pdf', width = 8, height = 11)
		layout(matrix(1:2, nrow = 2, byrow = TRUE))
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, cex = .6)
		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
		plot(GOA_NavyArea, border = 'black', add = TRUE, lty = 2, lwd = 2)
		text(300000,1696705,'Alaska', cex = 2.7, col = 'white')
		text(-1200000,800000,'Bering Sea', cex = 1.5, col = 'black')
		text(640000,800000,'Gulf of Alaska', cex = 1.5, col = 'black')
		lines(c(15000, 15000, 535000, 535000,15000), 
			c(990000, 1338000, 1338000, 990000, 990000), col = 'white', lwd = 3)
		text(-2250000,2600000, 'A', cex = 4)
				
		par(mar = c(6,8,0,8))
		plot(AKpoly[1:100,], col = 'grey50', xlim = c(15000,535000), 
			ylim = c(990000,1338000),)
		plot(POPhexagons,  border = 'black', add = TRUE, lwd = 1.5)
		plot(GOA_NavyArea, border = 'black', add = TRUE, lty = 2, lwd = 5)
		mtext('B', cex = 4, adj = -.19)
	dev.off()

	# crop figure
  system('pdfcrop Fig-study_area.pdf')


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-raw_data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


	pdf(file = paste0('Fig-raw_data.pdf'), width = 8, height = 14)
			
		layout(matrix(1:3, nrow = 3, byrow = TRUE))
		par(mar = c(0,0,0,0))
    #create some break points for shipdays
		fixedBreaks = c(1:2, 10, 30, 70, 
			max(ptDF$shipdays, na.rm = TRUE) + 1e-10) - 1e-12
    # find class intervals for shipdays based on breakpoints
		f10 = classIntervals(ptDF$shipdays, style = 'fixed', 
			fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Ship-Days', cex = 2.8, col = 'white')
		text(650000,1600000,'Raw Effort', cex = 2.8, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.8, printFormat = "4.0")
		text(-2250000,2600000,'A', cex = 4)

		fixedBreaks = c(1:2, 6, 12, 25, 
			max(ptDF$CU, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(ptDF$CU , style = 'fixed', fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Northern Fur Seals', cex = 2.8, col = 'white')
		text(650000,1600000,'Raw Counts', cex = 2.8, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.8, printFormat = "4.0")
		text(-2250000,2600000,'B', cex = 4)

		fixedBreaks = c(1:2, 10, 20, 50, 
			max(ptDF$EJ, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(ptDF$EJ , style = 'fixed', fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Steller Sea Lions', cex = 2.8, col = 'white')
		text(650000,1600000,'Raw Counts', cex = 2.8, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.8, printFormat = "4.0")
		text(-2250000,2600000,'C', cex = 4)
	dev.off()

	# crop figure
  system("pdfcrop Fig-raw_data.pdf")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-Fig-Effort
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

	#define mode function
	mode = function(x)
	{
		densout = density(x)
		densout$x[densout$y == max(densout$y)]	
	}
  
  shipdays = ptDF$shipdays
	# find which cells have no effort, nor any of their neighbors
	ind_no_nei_effort = sapply(nb_list,function(x){all(is.na(shipdays[x]))}) & 
		is.na(shipdays)
	# set these to zero effort
	shipdays[ind_no_nei_effort] = 0
	# create indicator of those with certain effort or zero effort
	ind_samp = which(!is.na(shipdays))
	# the rest are missing - unknown effort
	ind_miss = which(is.na(shipdays))

	#get posterior distribution of effort per cell, and mode per cell
	MCMC_effdays = matrix(NA, nrow = dim(MCMC_eff_zmat)[1], 
		ncol = dim(MCMC_eff_zmat)[2])
	for(i in 1:dim(MCMC_eff_zmat)[2])
		MCMC_effdays[,i] = exp(MCMC_eff_zmat[,i] + MCMC_eff_beta0[i])
	MCMC_eff_mean = apply(MCMC_effdays,1,mode)
	# create an indicator of effort greater than or equal to 1
	ind_eff = !is.na(ptDF$shipdays)  |  MCMC_eff_mean >= 1

	pdf(file = paste0('Fig-effort.pdf'), width = 10, height = 14)
			
		layout(matrix(1:2, nrow = 2, byrow = TRUE))
		
		fixedBreaks = c(0:2, 10, 30, 70, 
			max(shipdays, na.rm = TRUE) + 1e-10) - 1e-12
		ffb = classIntervals(shipdays, style = 'fixed', fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		pal[1] = 'maroon3'
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Ship-Days', cex = 1.7, col = 'white')
		text(650000,1600000,'Zeros Added', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.2, printFormat = "4.0")
		text(-2250000,2600000,'A', cex = 4)

		fixedBreaks = c(0:2, 10, 30, 70,
			max(MCMC_eff_mean, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(MCMC_eff_mean, style = 'fixed', 
			fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) -1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Ship-Days', cex = 1.7, col = 'white')
		text(650000,1600000,'Model Fit', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.0")
		text(-2250000,2600000,'B', cex = 4)
	dev.off()

	# crop figure
  system("pdfcrop Fig-effort.pdf")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-FurSealFit
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
 
 
	CU = ptDF$CU
	# if sufficient effort, but species not observed, set grid cell to zero
  CU[ind_eff & is.na(CU)] = 0
  # create indicator of those with certain effort or zero effort
  ind_samp = which(!is.na(CU))
  # the rest are missing - unknown effort
  ind_miss = which(is.na(CU))

	#posterior summary -- mode for each grid cell
  temp_CU = matrix(NA, nrow = dim(MCMC_CU_zmat)[1], ncol = dim(MCMC_CU_zmat)[2])
  for(i in 1:dim(MCMC_CU_zmat)[2])
	  temp_CU[,i] = exp(MCMC_CU_beta0[i] + MCMC_CU_zmat[,i])
  CU_post = apply(temp_CU,1,mode)
  CU_post = pmax(CU_post,0)

	pdf(file = paste0('Fig-FurSealFit.pdf'), width = 10, height = 14)
			
		layout(matrix(1:2, nrow = 2, byrow = TRUE))
		
		fixedBreaks = fixedBreaks = c(0:2, 6, 12, 25, 
			max(CU, na.rm = TRUE) + 1e-10) - 1e-12
		ffb = classIntervals(CU, style = 'fixed', fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		pal[1] = 'maroon3'
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Northern Fur Seals', cex = 1.7, col = 'white')
		text(650000,1600000,'Zeros Added', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.2, printFormat = "4.0")
		text(-2250000,2600000,'A', cex = 4)

    fixedBreaks = c(0, .1, .2, .4, .8, 3, max(CU_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(CU_post, style = 'fixed', fixedBreaks = fixedBreaks)
	#  f10 = classIntervals(CU_post, n = 6, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Northern Fur Seals', cex = 1.7, col = 'white')
		text(650000,1600000,'Model Fit', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'B', cex = 4)
	dev.off()

	# crop figure
  system("pdfcrop Fig-FurSealFit.pdf")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-SeaLionFit
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
 
 EJ = ptDF$EJ
	# if sufficient effort, but species not observed, set grid cell to zero
  EJ[ind_eff & is.na(EJ)] = 0
  # create indicator of those with certain effort or zero effort
  ind_samp = which(!is.na(EJ))
  # the rest are missing - unknown effort
  ind_miss = which(is.na(EJ))

	#posterior summary -- mode for each grid cell
  temp_EJ = matrix(NA, nrow = dim(MCMC_EJ_zmat)[1], ncol = dim(MCMC_EJ_zmat)[2])
  for(i in 1:dim(MCMC_EJ_zmat)[2])
	  temp_EJ[,i] = exp(MCMC_EJ_beta0[i] + MCMC_EJ_zmat[,i])
  EJ_post = apply(temp_EJ,1,mode)
  EJ_post[EJ_post < 0] = 0

	pdf(file = paste0('Fig-SeaLionFit.pdf'), width = 10, height = 14)
			
		layout(matrix(1:2, nrow = 2, byrow = TRUE))
		
		fixedBreaks = fixedBreaks = c(0:2, 10, 20, 50,  
			max(EJ, na.rm = TRUE) + 1e-10) - 1e-12
		ffb = classIntervals(EJ, style = 'fixed', fixedBreaks = fixedBreaks)
		pal = viridis(length(fixedBreaks) - 1)
		pal[1] = 'maroon3'
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = 'grey80', add = TRUE, pch = 19, cex = .3)
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Steller Sea Lions', cex = 1.7, col = 'white')
		text(650000,1600000,'Zeros Added', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.2, printFormat = "4.0")
		text(-2250000,2600000,'A', cex = 4)

		fixedBreaks = c(0, .1, .2, .5, 1, 5, 
			max(EJ_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(EJ_post, style = 'fixed', fixedBreaks = fixedBreaks)
	#  f10 = classIntervals(EJ_post, n = 9, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .3)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Steller Sea Lions', cex = 1.7, col = 'white')
		text(650000,1600000,'Model Fit', cex = 1.7, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'B', cex = 4)
	dev.off()

	# crop figure
  system("pdfcrop Fig-SeaLionFit.pdf")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-Uncertainty
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

	
  sd_zmat_eff = apply(MCMC_eff_zmat,1,sd)
  sd_eff = apply(MCMC_effdays,1,sd)
  sd_zmat_CU = apply(MCMC_CU_zmat,1,sd)
  sd_CU = apply(temp_CU,1,sd)
  sd_zmat_EJ = apply(MCMC_EJ_zmat,1,sd)
  sd_EJ = apply(temp_EJ,1,sd)

	pdf(file = paste0('Fig-Uncertainty.pdf'), width = 10, height = 14)
		
		layout(matrix(1:6, nrow = 3, byrow = TRUE))

		fixedBreaks = c(0, .3, .42, .54, .65, .82, 
			max(sd_zmat_eff) + 1e-10) - 1e-12
		ffb = classIntervals(sd_zmat_eff, style = 'fixed', 
			fixedBreaks = fixedBreaks)
	#  ffb = classIntervals(sd_zmat_eff, n = 6, style = 'fisher')
		pal = viridis(length(ffb$brks) - 1)
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Ship-Days', cex = 1.5, col = 'white')
		text(650000,1600000,'R-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'A', cex = 4)


		fixedBreaks = c(0, .7, 1.6, 2.8, 4.4, 7.5,
			max(sd_eff) + 1e-10) - 1e-12
		f10 = classIntervals(sd_eff, style = 'fixed', fixedBreaks = fixedBreaks)
	#  f10 = classIntervals(sd_eff, n = 6, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Ship-Days', cex = 1.5, col = 'white')
		text(650000,1600000,'Mean-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'B', cex = 4)

		fixedBreaks = c(0, .7, .9, 1.1, 1.3, 1.7, 
			max(sd_zmat_CU) + 1e-10) - 1e-12
		ffb = classIntervals(sd_zmat_CU, style = 'fixed', fixedBreaks = fixedBreaks)
	#  ffb = classIntervals(sd_zmat_CU, n = 6, style = 'fisher')
		pal = viridis(length(ffb$brks) - 1)
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Northern Fur Seals', cex = 1.5, col = 'white')
		text(650000,1600000,'Z-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'C', cex = 4)

		fixedBreaks = c(0, .2, .5, 1.2, 2.2, 3.4, 
			max(sd_CU) + 1e-10) - 1e-12
		ffb = classIntervals(sd_CU, style = 'fixed', fixedBreaks = fixedBreaks)
	#  ffb = classIntervals(sd_CU, n = 6, style = 'fisher')
		pal = viridis(length(ffb$brks) - 1)
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Northern Fur Seals', cex = 1.5, col = 'white')
		text(650000,1600000,'Mean-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'D', cex = 4)

		fixedBreaks = c(0, .9, 1.1, 1.3, 1.5, 2,
			max(sd_zmat_EJ) + 1e-10) - 1e-12
		ffb = classIntervals(sd_zmat_EJ, style = 'fixed', fixedBreaks = fixedBreaks)
	#  ffb = classIntervals(sd_zmat_EJ, n = 6, style = 'fisher')
		pal = viridis(length(ffb$brks) - 1)
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Steller Sea Lions', cex = 1.5, col = 'white')
		text(650000,1600000,'Z-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'E', cex = 4)

		fixedBreaks = c(0, .3, .8, 1.5, 2.5, 4.6, 
			max(sd_EJ) + 1e-10) - 1e-12
		ffb = classIntervals(sd_EJ, style = 'fixed', fixedBreaks = fixedBreaks)
	#  ffb = classIntervals(sd_EJ, n = 6, style = 'fisher')
		pal = viridis(length(ffb$brks) - 1)
		ffbColours = findColours(ffb, pal)
		par(mar = c(0,0,0,0))
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = ffbColours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], add = TRUE, border = 'black')
		text(650000,1850000,'Steller Sea Lions', cex = 1.5, col = 'white')
		text(650000,1600000,'Mean-Uncertainty', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, ffb$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'F', cex = 4)

	dev.off()

	# crop figure
  system("pdfcrop Fig-Uncertainty.pdf")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-GOAabu
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

  overlay_indx = unlist(over(GOA_NavyArea, trueCentroids, returnList = TRUE))

  Nhat_CU = 620660
	CU_norm = matrix(NA, nrow = dim(MCMC_CU_zmat)[1], ncol = dim(MCMC_CU_zmat)[2])
  for(i in 1:dim(MCMC_CU_zmat)[2])
	  CU_norm[,i] = temp_CU[,i]/sum(temp_CU[,i])
	CU_dens = CU_norm*Nhat_CU
  CU_dens_post = apply(CU_dens,1,mode)
  CU_dens_post = pmax(CU_dens_post,0)

	CU_GOA_abu = NULL
  for(i in 1:dim(MCMC_CU_zmat)[2])
	  CU_GOA_abu = c(CU_GOA_abu, sum(CU_dens[overlay_indx,i]))

  Nhat_E = 20756 + 7838
  Nhat_W =  54267
  Nhat_EJ = Nhat_E + Nhat_W

	EJ_norm = matrix(NA, nrow = dim(MCMC_EJ_zmat)[1], ncol = dim(MCMC_EJ_zmat)[2])
  for(i in 1:dim(MCMC_EJ_zmat)[2])
	  EJ_norm[,i] = temp_EJ[,i]/sum(temp_EJ[,i])
	EJ_dens = EJ_norm*Nhat_EJ
  EJ_dens_post = apply(EJ_dens,1,mode)
  EJ_dens_post[EJ_dens_post < 0] = 0

	EJ_GOA_abu = NULL
  for(i in 1:dim(MCMC_EJ_zmat)[2])
	  EJ_GOA_abu = c(EJ_GOA_abu, sum(EJ_dens[overlay_indx,i]))
 
	pdf(file = paste0('Fig-GOAabu.pdf'), width = 9, height = 9)

		layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(3,1,3))

		fixedBreaks = c(0, 15, 30, 50, 100, 
			max(CU_dens_post[overlay_indx], na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(CU_dens_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
	#  f10 = classIntervals(CU_dens_post[overlay_indx], n = 6, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		par(mar = c(0,5,7,0))
		plot(AKpoly[1:100,], col = 'grey20', xlim = c(-115000,655000), 
				ylim = c(10000,1600000))
		plot(trueCentroids, col = 'grey60', add = TRUE, 
			lwd = .5, pch = 19, cex = .6)
		plot(trueCentroids[overlay_indx,], col = f10Colours[overlay_indx], 
			add = TRUE, lwd = .5, pch = 19, cex = .6)
#		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
# 	text(-390000, 2528054,'A', cex = 6)
		mtext('A', adj = -.1, cex = 3)
		par(mar = c(0,0,0,0))
		plot(c(0,1), c(0,1), type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n')
		addBreakColorLegend(.01, .2, .3, .8, f10$brks, 
				colors = pal, cex = 1.5, printFormat = "4.0")

		par(mar = c(5,7,7,1))
		hist(CU_GOA_abu, col = 'blue',
			main = '', cex.lab = 2, cex.axis = 1.5,
			xlab = 'Total Northern Fur Seals')
		mtext('B', adj = -.15, cex = 3)


		fixedBreaks = c(0, 2, 6, 17, 37,
			max(EJ_dens_post[overlay_indx], na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(EJ_dens_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
	#  classIntervals(EJ_dens_post[overlay_indx], n = 8, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		par(mar = c(0,5,7,0))
		plot(AKpoly[1:100,], col = 'grey20', xlim = c(-115000,655000), 
				ylim = c(10000,1600000))
		plot(trueCentroids, col = 'grey60', add = TRUE, 
			lwd = .5, pch = 19, cex = .6)
		plot(trueCentroids[overlay_indx,], col = f10Colours[overlay_indx], 
			add = TRUE, lwd = .5, pch = 19, cex = .6)
#		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
# 	text(-390000, 2528054,'A', cex = 6)
		mtext('C', adj = -.1, cex = 3)
		par(mar = c(0,0,0,0))
		plot(c(0,1), c(0,1), type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n')
		addBreakColorLegend(.01, .2, .3, .8, f10$brks, 
				colors = pal, cex = 1.5, printFormat = "4.0")

		par(mar = c(5,7,7,1))
		hist(EJ_GOA_abu, col = 'blue',
			main = '', cex.lab = 2, cex.axis = 1.5,
			xlab = 'Total Steller Sea Lions')
		mtext('D', adj = -.15, cex = 3)

	dev.off()
	
	# crop figure
  system("pdfcrop Fig-GOAabu.pdf")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-HotSpots
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
  
	l
	smoo_CU_post = apply(MCMC_CU_smoo,1,mean)
	smoosd_CU_post = apply(MCMC_CU_smoo,1,sd)
	smoo_EJ_post = apply(MCMC_EJ_smoo,1,mean)
	smoosd_EJ_post = apply(MCMC_EJ_smoo,1,sd)

 
	pdf(file = 'Fig-SmooHotSpots.pdf', width = 11, height = 11)

		layout(matrix(1:4, nrow = 2, byrow = TRUE))
			
		par(mar = c(0,0,0,0))
		fixedBreaks = c(0, .13, .30, .57, .98, 1.92, 
			max(smoo_CU_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(smoo_CU_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
	# f10 = classIntervals(smoo_CU_post, n = 6, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, 
			add = TRUE, lwd = .5, pch = 19, cex = .2)
#		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
		text(650000,1850000,'Northern Fur Seals', cex = 1.5, col = 'white')
		text(650000,1600000,'Smoothed Model Fit', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'A', cex = 4)

		high_CU_post = smoo_CU_post
		high_CU_post[smoo_CU_post < quantile(smoo_CU_post,.90)] = 0
	
		fixedBreaks = c(0, abs(qnorm(.025/length(smoosd_CU_post))), 6, 8, 
			max(high_CU_post/smoosd_CU_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(high_CU_post/smoosd_CU_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
		pal = viridis(length(f10$brks) -1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
#		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
		text(650000,1850000,'Northern Fur Seals', cex = 1.5, col = 'white')
		text(650000,1600000,'Certain Hotspots', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
			colors = pal, cex = 1.1, printFormat = "4.2")
		text(-2250000,2600000,'B', cex = 4)

		fixedBreaks = c(0, .17, .40, .71, 1.10, 1.62, 
			max(smoo_EJ_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(smoo_EJ_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
	# f10 = classIntervals(smoo_EJ_post, n = 6, style = 'fisher')
		pal = viridis(length(f10$brks) - 1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, 
			add = TRUE, lwd = .5, pch = 19, cex = .2)
#		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
		text(650000,1850000,'Steller Sea Lions', cex = 1.5, col = 'white')
		text(650000,1600000,'Smoothed Model Fit', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
				colors = pal, cex = 1.1, printFormat = "4.1")
		text(-2250000,2600000,'C', cex = 4)

		high_EJ_post = smoo_EJ_post
		high_EJ_post[smoo_EJ_post < quantile(smoo_EJ_post,.90)] = 0
	
   fixedBreaks = c(0, abs(qnorm(.025/length(smoosd_EJ_post))), 5.5, 7, 
		max(high_EJ_post/smoosd_EJ_post, na.rm = TRUE) + 1e-10) - 1e-12
		f10 = classIntervals(high_EJ_post/smoosd_EJ_post, style = 'fixed', 
			fixedBreaks = fixedBreaks)
		pal = viridis(length(f10$brks) -1)
		f10Colours = findColours(f10, pal)
		plot(AKpoly[1:100,], col = 'grey20')
		plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, 
			pch = 19, cex = .2)
		plot(AKpoly[1:100,], col = 'grey20', add = TRUE)
		text(650000,1850000,'Steller Sea Lions', cex = 1.5, col = 'white')
		text(650000,1600000,'Certain Hotspots', cex = 1.5, col = 'white')
		addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
			colors = pal, cex = 1.1, printFormat = "4.2")
		text(-2250000,2600000,'D', cex = 4)

	dev.off()

	# crop figure
  system("pdfcrop Fig-SmooHotSpots.pdf")

