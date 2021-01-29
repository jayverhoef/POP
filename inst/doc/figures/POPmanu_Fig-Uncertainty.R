# Rscript /home/xverhoef/data/2019_papers/POP/POP_package/POP/inst/doc/POPmanu_Fig-Uncertainty.R
  library(POP)
	library(viridis)
	library(classInt)
	library(Matrix)
	library(sp)
	library(rgeos)

	data(grid)
	data(nb_list)
	data(AKpoly)
	data(GOA_NavyArea)
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

	trueCentroids = gCentroid(grid,byid=TRUE)
	
setwd(paste0('/home/xverhoef/data/2019_papers/POP/POP_package/POP',
		'/inst/doc/'))
setwd(paste0('/media/jay/data/desktop_data/2019_papers/POP/POP_package/POP',
		'/inst/doc/figures/'))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Fig-Uncertainty
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

	mode = function(x)
	{
		densout = density(x)
		densout$x[densout$y == max(densout$y)]	
	}

	#get posterior distribution of effort per cell, and mode per cell
	MCMC_effdays = matrix(NA, nrow = dim(MCMC_eff_zmat)[1], 
		ncol = dim(MCMC_eff_zmat)[2])
	for(i in 1:dim(MCMC_eff_zmat)[2])
		MCMC_effdays[,i] = exp(MCMC_eff_zmat[,i] + MCMC_eff_beta0[i])
	MCMC_eff_mean = apply(MCMC_effdays,1,mode)

  temp_CU = matrix(NA, nrow = dim(MCMC_CU_zmat)[1], ncol = dim(MCMC_CU_zmat)[2])
  for(i in 1:dim(MCMC_CU_zmat)[2])
	  temp_CU[,i] = exp(MCMC_CU_beta0[i] + MCMC_CU_zmat[,i])
  CU_post = apply(temp_CU,1,mode)

  temp_EJ = matrix(NA, nrow = dim(MCMC_EJ_zmat)[1], ncol = dim(MCMC_EJ_zmat)[2])
  for(i in 1:dim(MCMC_EJ_zmat)[2])
	  temp_EJ[,i] = exp(MCMC_EJ_beta0[i] + MCMC_EJ_zmat[,i])
  EJ_post = apply(temp_EJ,1,mode)

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
  system("pdfcrop /media/jay/data/desktop_data/2019_papers/POP/POP_package/POP/inst/doc/figures/Fig-Uncertainty.pdf")


