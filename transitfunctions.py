#####
#File: transitfunctions.py
#Purpose: Contains functions that relate to the transit curves of transits.
#:-) Started 3-14-15
#
#####

#Below imports necessary packages
import numpy as np
import math as calc
import random as rand
import matplotlib.pyplot as graph
import matplotlib.patches as mpatch
import matplotlib.cm as cm
import astropy.constants as const
#Below imports necessary functions from elsewhere
from starutils.populations import BGStarPopulation_TRILEGAL
from starutils.populations import MultipleStarPopulation
from starutils.populations import DARTMOUTH
bandlist = ['g', 'r', 'i', 'z', 'J', 'H', 'K', 'Kepler']

#Below imports necessary constants
#Below are some constant values
daysecs = 86400.0 #Number of seconds in a day
yearsecs = daysecs*365.256 #Number of seconds in a year
yeardays = 365.256 #Number of days in a year
masssun = const.M_sun.cgs.value #Mass of Sun in grams
massearth = const.M_earth.cgs.value #Mass of Earth in grams
Gconst = const.G.cgs.value #Gravitational constant, cgs
rsun = const.R_sun.cgs.value #Radius of Sun, cgs
rearth = const.R_earth.cgs.value #Radius of Earth, cgs
eccearth = 0.01671123 #Eccentricty of Earth's orbit; ###FIND EXACT VALUE!
au = const.au.cgs.value #Astronomical unit, cgs
pi = np.pi #Value of pi


##FUNCTION: gensim
#Purpose: This function is meant to run a star simulation to generate a population of stars and statistics
#NOTE: Sections of code heavily (very heavily) adapted from Tim Morton's vespa.populations code; thanks!
#Dependencies: starutils.populations.BGStarPopulation_TRILEGAL
def gensim(trilegal_filename, ra=None, dec=None, n=2e4, ichrone=DARTMOUTH, MAfn=None, mags=None, maxrad=10, f_binary=0.4, **kwargs):
	##
	#BELOW SECTION: Borrowed code from Tim Morton
	bgpop = BGStarPopulation_TRILEGAL(trilegal_filename, ra=ra, dec=dec, mags=mags, maxrad=maxrad, n=n, **kwargs)

	# Make sure that
    # properties of stars are within allowable range for isochrone.
    # This is a bit hacky, admitted.
	mass = bgpop.stars['m_ini'].values
	age = bgpop.stars['logAge'].values
	feh = bgpop.stars['[M/H]'].values

	pct = 0.05 #pct distance from "edges" of ichrone interpolation
	mass[mass < ichrone.minmass*(1+pct)] = ichrone.minmass*(1+pct)
	mass[mass > ichrone.maxmass*(1-pct)] = ichrone.maxmass*(1-pct)
	age[age < ichrone.minage*(1+pct)] = ichrone.minage*(1+pct)
	age[age > ichrone.maxage*(1-pct)] = ichrone.maxage*(1-pct)
	feh[feh < ichrone.minfeh+0.05] = ichrone.minfeh+0.05
	feh[feh > ichrone.maxfeh-0.05] = ichrone.maxfeh-0.05

	distance = bgpop.stars['distance'].values

	#Generate binary population to draw eclipses from
	pop = MultipleStarPopulation(mA=mass, age=age, feh=feh, f_triple=0, f_binary=1, distance=distance, n=n, ichrone=ichrone)
	
	#End of borrowed code
	##
	
	#Below returns generated population
	return pop
##


#####
#BELOW SECTION: Opens up the star distribution and records necessary values
#Below generates a sample population
filename = 'samplepop.h5'
pop = gensim(filename)
#####



##FUNCTION: genstarpop
#Purpose: This function is meant to generate a population of stars from passed in population of stars; will randomize the values by shuffling the given arrays
def genstarpop(population=pop, periodp=3, ecca=2, eccb=2, numruns=int(1e4), band='Kepler', plot=True, savename='exrunstar', **kwargs):
	#BELOW SECTION: Puts together the magnitudes
	if band not in bandlist:
		raise ValueError("Oh no!  Looks like you passed an invalid band.  Please pass one of the these bands:\n" + str(bandlist) + "\n(Note that the default is 'Kepler'.)")
	elif band in bandlist:
		magkeyword = str(band) + '_mag'
		magraw = population.stars[magkeyword]
	
	#Below are star parameters, unshuffled
	rAraw = population.stars.radius_A
	rBraw = population.stars.radius_B
	massAraw = population.stars.mass_A
	massBraw = population.stars.mass_B
	uAraw = population.stars.Teff_A
	uBraw = population.stars.logg_A
	
	#BELOW SECTION: Generates shuffled arrays of population data; allows duplicates
	#Below generates an array of random indices; allows replacement; similar to bootstrapping technique
	shuffledinds = np.random.randint(len(rAraw), size=(numruns*5))
	
	#Below shuffles the orbital parameters with the shuffled indices
	rAshuffle = np.array(rAraw[shuffledinds])
	rBshuffle = np.array(rBraw[shuffledinds])
	massAshuffle = np.array(massAraw[shuffledinds])
	massBshuffle = np.array(massBraw[shuffledinds])
	uAshuffle = np.array(uAraw[shuffledinds])
	uBshuffle = np.array(uBraw[shuffledinds])
	magshuffle = np.array(magraw[shuffledinds])

	
	#BELOW SECTION: Runs model on the shuffled values
	model = runmodel(rAraw=rAshuffle, rBraw=rBshuffle, massAraw=massAshuffle, massBraw=massBshuffle, uAraw=uAshuffle, uBraw=uBshuffle, magraw=magshuffle, periodp=periodp, ecca=ecca, eccb=eccb, numruns=numruns, band=band, plot=plot, savename=savename, poptype='Star')
	return model
	
##


##FUNCTION: genplanetpop
#Purpose: This function is meant to generate a population of planets based upon read-in primary star values
def genplanetpop(population=pop, rArange=.10, massArange=.10, periodp=3, ecca=2, eccb=2, numruns=int(1e4), plot=True, radiusratio=1.0, radiusbinsize=0.3, savename='exrunplanet', band='Kepler', **kwargs):	
	#BELOW SECTION: Puts together the magnitudes
	if band not in bandlist:
		raise ValueError("Oh no!  Looks like you passed an invalid band.  Please pass one of the these bands:\n" + str(bandlist) + "\n(Note that the default is 'Kepler'.)")
	elif band in bandlist:
		magkeyword = str(band) + '_mag'
		magraw = population.stars[magkeyword]

	#Below are star parameters, unshuffled
	uAraw = population.stars.Teff_A
	uBraw = population.stars.logg_A
	
	#BELOW SECTION: Generates shuffled arrays of population data; allows duplicates
	#Below generates an array of random indices; allows replacement; similar to bootstrapping technique
	shuffledinds = np.random.randint(len(uAraw), size=(numruns*5))
	
	#Below shuffles the orbital parameters using the shuffled indices, with replacement
	uAshuffle = np.array(uAraw[shuffledinds])
	uBshuffle = np.array(uBraw[shuffledinds])
	magshuffle = np.array(magraw[shuffledinds])
	
	#BELOW SECTION: Generates random primary radii and masses
	#Follows a Gaussian distribution about Sun's mass and Sun's radius using passed in uncertainties above
	rAshuffle = np.random.normal(loc=rsun, scale=(rArange/4.0)*rsun, size=(numruns*5))/(rsun*1.0) #In units of solar radii
	massAshuffle = np.random.normal(loc=masssun, scale=(massArange/4.0)*masssun, size=(numruns*5))/(masssun*1.0) #In units of solar masses
	
	#BELOW SECTION: Generates random planetary radii and masses, based upon ratio of Earth's radius to Sun's radius
	#NOTE: Mostly follows the procedures for generating planetary populations as given in T. Morton's 'PlanetPopulation' code
	#Below generates radius bin from which to uniformly draw radii
	baseratio = (rearth*radiusratio/(rsun*1.0))*np.mean(np.ma.masked_array(rAshuffle, np.isnan(rAshuffle))) #Excludes nans
	radiusbinmax = baseratio*(1 + radiusbinsize) #Max radius in bin
	radiusbinmin = baseratio*(1 - radiusbinsize) #Min radius in bin
	
	#Below pulls radii randomly from given bin
	rBshuffle = (np.random.random(len(rAshuffle)))*(radiusbinmax - radiusbinmin) + radiusbinmin
	#Below next generates masses based upon radii; same formula as used by T. Morton's planetary mass generation
	massBshuffle = ((rBshuffle*rsun/(1.0*rearth))**2.06) * (massearth/(1.0*masssun))

	
	#BELOW SECTION: Runs model on the shuffled values
	model = runmodel(rAraw=rAshuffle, rBraw=rBshuffle, massAraw=massAshuffle, massBraw=massBshuffle, uAraw=uAshuffle, uBraw=uBshuffle, magraw=magshuffle, periodp=periodp, ecca=ecca, eccb=eccb, numruns=numruns, band='Kepler', poptype='Planet', radiusratio=radiusratio, plot=plot, savename=savename)
	return model
##



##FUNCTION: varypop
#Purpose: This function is meant to plot populations with varied parameters, such as varied eccentricity inputs or varied period distributions.
#Requires that the parameters periodp, ecca, and eccb be passed in as arrays of the same size.  Will treat the parameter arrays as a giant matrix and will overplot each set of entries for each singular plot.
#For example, the first entries of each parameter array would be used to generate the first population; would plot for that population.
def varypop(population=pop, poptype=['Star', 'Star', 'Star'], periodp=[3,3,3], ecca=[2,2,2], eccb=[2,2,2], numruns=int(1e4), radiusratio=[1,1,1], radiusbinsize=0.3, savenameroot='exrunmulti', band='Kepler', **kwargs):
	#BELOW SECTION: Some error control
	#Below throws an error if passed-in values are not in list or np.array form
	for part in [periodp, ecca, eccb]:
		if not isinstance(part, list) and not isinstance(part, np.ndarray):
			raise ValueError("Wait!  Please pass in all parameters periodp, ecca, and eccb in either list or numpy array form.  The parameter poptype should of course be a list.")
	#Below throws an error if mismatching parameter array lengths were passed in
	if len(periodp) != len(ecca) or len(ecca) != len(eccb) or len(eccb) != len(periodp) or len(poptype) != len(periodp):
		raise ValueError("Whoa!  It seems you passed in parameter arrays with mismatching lengths.  Please make sure your parameter arrays all have the same length.")
	numhere = len(periodp) #Number of parameter sets
	cols = cm.brg(np.linspace(0, 1, numhere)) #Colors to use
	graph.close() #Just in case
	
	
	#BELOW SECTION: Generates plots for each parameter combination
	savename = savenameroot #For recording graphs #For saving graphs
	popdict = {}
	for a in range(0, numhere):
		#BELOW SECTION: Generates results for populations for the current parameters
		pophere = None #Declares variable
		if poptype[a] == 'Star':
			pophere = genstarpop(population=pop, periodp=periodp[a], ecca=ecca[a], eccb=eccb[a], numruns=numruns, savename=(str(savenameroot)+'star'), band=band, plot=False)
		elif poptype[a] == 'Planet':
			pophere = genplanetpop(population=pop, periodp=periodp[a], ecca=ecca[a], eccb=eccb[a], numruns=numruns, savename=(str(savenameroot)+'planet'), radiusratio=float(radiusratio[a]), band=band, plot=False)
		else:
			raise ValueError("Looks like you passed an invalid value for poptype.")
		
		#Below adds these populations to the dictionary
		popdict[str(a)] = pophere
		
	
	#BELOW SECTION: Carries out the graphs to compare each set of parameters
	#Hist: Periods
	for b in range(0, numhere):
		popformhere = None #Population form name
		if poptype[b] == 'Star':
			popformhere = poptype[b] + ': '
		elif poptype[b] == 'Planet':
			popformhere = poptype[b] + ' (' + str(radiusratio[b]) + '$R_E$): '
		labelhere = (popformhere+'p='+str(popdict[str(b)]['periodp']))
		graph.hist(popdict[str(b)]['period']/yearsecs, histtype='step', color=cols[b], label=labelhere)
		
	graph.title('Period Distribution (in Years)')
	graph.xlabel('Period (in Years)')
	graph.ylabel('Counts')
	graph.legend(loc=8, prop={'size':14})
	graph.savefig(savename+'PerHist.png')
	graph.close()
	#
	
	
	#Hist: Eccentricities
	for c in range(0, numhere):
		popformhere = None #Population form name
		if poptype[c] == 'Star':
			popformhere = poptype[c] + ': '
		elif poptype[c] == 'Planet':
			popformhere = poptype[c] + ' (' + str(radiusratio[c]) + '$R_E$): '
		labelhere = (popformhere+'a='+str(popdict[str(c)]['ecca'])+', b='+str(popdict[str(c)]['eccb']))
		graph.hist(popdict[str(c)]['ecc'], histtype='step',  color=cols[c],label=labelhere)
	
	graph.title('Eccentricity Distribution')
	graph.xlabel('Eccentricity')
	graph.ylabel('Counts')
	graph.legend(loc='best', prop={'size':14})
	graph.savefig(savename+'EccHist.png')
	graph.close()
	#
	
	
	#Hist: T14 (in days); Outliers cut out
	for c in range(0, numhere):
		popformhere = None #Population form name
		if poptype[c] == 'Star':
			popformhere = poptype[c] + ': '
		elif poptype[c] == 'Planet':
			popformhere = poptype[c] + ' (' + str(radiusratio[c]) + '$R_E$): '
		xhere = popdict[str(c)]['T14']/daysecs
		labelhere = (popformhere+'a='+str(popdict[str(c)]['ecca'])+', b='+str(popdict[str(c)]['eccb'])+', p='+str(popdict[str(c)]['periodp']))
		graph.hist(xhere[xhere < 5], histtype='step',  color=cols[c],label=labelhere)
	
	graph.title('T14 (in Days) Distribution')
	graph.xlabel('T14 (in Days)')
	graph.ylabel('Counts')
	graph.legend(loc='best', prop={'size':14})
	graph.savefig(savename+'T14Hist.png')
	graph.close()
	#
	
	
	#Hist: T14/tau; Outliers cut out
	for c in range(0, numhere):
		popformhere = None #Population form name
		if poptype[c] == 'Star':
			popformhere = poptype[c] + ': '
		elif poptype[c] == 'Planet':
			popformhere = poptype[c] + ' (' + str(radiusratio[c]) + '$R_E$): '
		xhere = popdict[str(c)]['T14']/popdict[str(c)]['tau']
		labelhere = (popformhere+'a='+str(popdict[str(c)]['ecca'])+', b='+str(popdict[str(c)]['eccb'])+', p='+str(popdict[str(c)]['periodp']))
		graph.hist(xhere[xhere < 10], histtype='step',  color=cols[c],label=labelhere)
	
	graph.title('T14/Tau Distribution')
	graph.xlabel('T14/Tau')
	graph.ylabel('Counts')
	graph.legend(loc='best', prop={'size':14})
	graph.savefig(savename+'T14-tauRatioHist.png')
	graph.close()
	#
	

	
	
	#BELOW SECTION: Plots based in form upon background research papers (specifically, on the 'Exoplanet Transit Validations...' Paper	
	#Scatter: T14 (in days) vs. log10(Exact Depth), Zoomed
	patches = []
	for d in range(0, numhere):
		popformhere = None #Population form name
		if poptype[d] == 'Star':
			popformhere = poptype[d] + ': '
		elif poptype[d] == 'Planet':
			popformhere = poptype[d] + ' (' + str(radiusratio[d]) + '$R_E$): '
		labelhere = (popformhere+'a='+str(popdict[str(d)]['ecca'])+', b='+str(popdict[str(d)]['eccb'])+', p='+str(popdict[str(d)]['periodp']))
		
		#Below takes care of labeling
		patches.append(mpatch.Patch(color=cols[d], label=labelhere))
		x = popdict[str(d)]['T14']/daysecs
		y = np.log10(popdict[str(d)]['depthexact'])
		lins = np.linspace(0,1,numhere)
		graph.scatter(x, y, facecolors='black', edgecolors=cols[d], s=2,alpha=(1-((d+1)/(len(lins)*1.0+1)))*.2)
		
	graph.title('T14 (in Days) vs. Log10(Depth, Exact)')
	graph.xlabel('T14 in Days')
	graph.ylabel('Log10(Depth, Exact)')
	graph.xlim([-1, 5])
	graph.legend(loc='best', handles=patches, prop={'size':10})
	graph.savefig(savename+'T14vsLog10DepthExactZoom.png')
	graph.close()
	#
	
	

	#Scatter: T14 (in days) vs. T14/tau, Zoomed
	patches = []
	for e in range(0, numhere):
		popformhere = None #Population form name
		if poptype[e] == 'Star':
			popformhere = poptype[e] + ': '
		elif poptype[e] == 'Planet':
			popformhere = poptype[e] + ' (' + str(radiusratio[e]) + '$R_E$): '
		labelhere = (popformhere+'a='+str(popdict[str(e)]['ecca'])+', b='+str(popdict[str(e)]['eccb'])+', p='+str(popdict[str(e)]['periodp']))
	
		#Below takes care of labeling
		patches.append(mpatch.Patch(color=cols[e], label=labelhere))
		x = popdict[str(e)]['T14']/daysecs
		y = (popdict[str(e)]['T14'])/(popdict[str(e)]['tau'])
		lins = np.linspace(0,1,numhere)
		graph.scatter(x, y, facecolors='black', edgecolors=cols[e], s=2,alpha=(1-((e+1)/(len(lins)*1.0+1)))*.2)
	
	graph.title('T14 (in Days) vs. T14/tau Ratio')
	graph.xlabel('T14 in Days')
	graph.ylabel('T14/tau Ratio')
	graph.xlim([-.5, 5])
	graph.ylim([-1, 150])
	graph.legend(loc='upper left', handles=patches, prop={'size':10})
	graph.savefig(savename+'T14vsT14-tauRatio.png')
	graph.close()
	
	
	
	#Scatter: Log10(Depth, Exact) vs. T14/tau, Zoomed
	patches = []
	for f in range(0, numhere):
		popformhere = None #Population form name
		if poptype[f] == 'Star':
			popformhere = poptype[f] + ': '
		elif poptype[f] == 'Planet':
			popformhere = poptype[f] + ' (' + str(radiusratio[f]) + '$R_E$): '
		labelhere = (popformhere+'a='+str(popdict[str(f)]['ecca'])+', b='+str(popdict[str(f)]['eccb'])+', p='+str(popdict[str(f)]['periodp']))
		
		#Below takes care of labeling
		patches.append(mpatch.Patch(color=cols[f], label=labelhere))
		x = np.log10(popdict[str(f)]['depthexact'])
		y = (popdict[str(f)]['T14'])/(popdict[str(f)]['tau'])
		lins = np.linspace(0,1,numhere)
		graph.scatter(x, y, facecolors='black', edgecolors=cols[f], s=2,alpha=(1-((e+1)/(len(lins)*1.0+1)))*.2)
		
	graph.title('Log10(Depth, Exact) vs. T14/tau Ratio')
	graph.xlabel('Log10(Depth, Exact)')
	graph.ylabel('T14/tau Ratio')
	graph.xlim([-8, 1])
	graph.ylim([-10, 100])
	graph.legend(loc='best', handles=patches, prop={'size':10})
	graph.savefig(savename+'Log10DepthExactvsT14-tauRatioZoom.png')
	graph.close()
	
	
	#Returns finished dictionary of populations
	return popdict
##



##FUNCTION: runmodel
#Purpose: This function is meant to run a simulation to generate a population of statistics based upon passed in radii, masses, and (if needed) darkening parameters, and then perform transit calculations.  It accepts one value to generate period distributions (power law distribution to power periodp) and two values to generate eccentricity distributions (eccentricity distribution with alpha ecca and beta eccb).
#Numruns indicates the number of stars.
#Dependencies:
def runmodel(rAraw, rBraw, massAraw, massBraw, uAraw, uBraw, magraw, periodp=3, ecca=2, eccb=2, numruns=int(1e4), radiusratio=None, poptype=None, plot=True, savename='exrun', band='Kepler', **kwargs):
	#Note: totalruns specifies more than enough data to be used
	#This allows room to screen out 'nan' values later in calctransit
	totalruns = calc.ceil(numruns*1.01)
	
	#BELOW SECTION: Generates random parameter distributions
	eccs = makeecc(totalruns, ecca, eccb) #Beta distribution
	omegas = makeangle(totalruns)
	impsraw = makeimp(rAraw, rBraw)
	
	#Below calculates omegas; throws out any invalid omega values
	angles = np.zeros(totalruns)
	rA = np.zeros(totalruns)
	rB = np.zeros(totalruns)
	massA = np.zeros(totalruns)
	massB = np.zeros(totalruns)
	uA = np.zeros(totalruns)
	uB = np.zeros(totalruns)
	mag = np.zeros(totalruns)
	periods = np.zeros(totalruns)
	imps = np.zeros(totalruns)
	semitesting = np.zeros(totalruns) #For TESTING purposes
	
	#BELOW SECTION: Ensures only valid periods are produced
	#Also ensures only valid radii are used
	track = 0 #Track number of valid calculations
	tryhere = 0 #Track number of tried periods
	totalplace = 0 #Track place in original total arrays
	while track < totalruns:
		#Below selects larger radius from rA and rB
		rhere = 0
		rnothere = 0
		#Below skips current place in total arrays if zero-valued radii
		if rAraw[totalplace] == 0 and rBraw[totalplace] == 0:
			totalplace = totalplace + 1
			continue
		#Below skips current place if nans in radii
		if calc.isnan(rAraw[totalplace]) or calc.isnan(rBraw[totalplace]):
			totalplace = totalplace + 1
			continue
		
		#Below assigns larger and smaller radii at this spot
		if rAraw[totalplace] > rBraw[totalplace]:
			rhere = rAraw[totalplace]
			rnothere = rBraw[totalplace]
		elif rAraw[totalplace] <= rBraw[totalplace]:
			rhere = rBraw[totalplace]
			rnothere = rAraw[totalplace]
		
		
		#Below generates period and semi here
		periodhere = makeperiod(1, periodp, rangestart=yeardays*4, rangeend=yeardays*20, dayunits=False) #Period between about 4 and 20 years
		semihere = calcsemimajor(period=periodhere, massstar=massAraw[totalplace]+massBraw[totalplace])
		
		#Below ensures period as within the roche limit
		while withinroche(semihere, mass1=massAraw[totalplace], mass2=massBraw[totalplace], r1=rhere, r2=rnothere):
			#Below generates new periods and semis until valid distance
			periodhere = makeperiod(1, periodp, rangestart=yeardays*4, rangeend=yeardays*20, dayunits=False) #Period between about 4 and 20 years
			semihere = calcsemimajor(period=periodhere, massstar=massAraw[totalplace]+massBraw[totalplace])
			
			#Below increments count of tries
			tryhere = tryhere + 1
			#Skips ahead if too many tries
			if tryhere >= 1000:
				print('No periods within roche range, apparently.')
				#print('tryhere reached ', tryhere)
				#print('totalplace reached ', totalplace)
				#print('Current track as ', track)
				print('So, skipping this data point at totalplace ', totalplace)
				print('')
				totalplace = totalplace + 1
				tryhere = 0
				break
				#raise ValueError("Oh no!  Used too many tries for generating a valid period at this point.")
				


		#Below computes semi and angle with current values
		semi = semihere
		#anglehere = calcminangle(period=periodhere, mass1=massAraw[totalplace], mass2=massBraw[totalplace], r1=rhere, r2=rnothere)
		anglehere = calcangle(period=periodhere, r1=rhere, semi=semihere, omega=omegas[track]*pi/180.0, imp=impsraw[totalplace], ecc=eccs[track])

		#Below keeps data if valid angle parameter was produced; throws out otherwise
		if not calc.isnan(anglehere):
			angles[track] = anglehere
			semitesting[track] = semi
			imps[track] = impsraw[totalplace]
			rA[track] = rAraw[totalplace]
			rB[track] = rBraw[totalplace]
			massA[track] = massAraw[totalplace]
			massB[track] = massBraw[totalplace]
			uA[track] = uAraw[totalplace]
			uB[track] = uBraw[totalplace]
			periods[track] = periodhere
			mag[track] = magraw[totalplace]
			
			#Below increments counts
			track = track + 1
			totalplace = totalplace + 1
			tryhere = 0

		else: #If too many tries used, skips to next data
			tryhere = tryhere + 1
			if tryhere >= 1000:
				print('tryhere reached ', tryhere, ' at totalplace of ', totalplace)
				#print('totalplace reached ', totalplace)
				#print('Current track as ', track)
				print('So, skipping this data point.')
				print('')
				totalplace = totalplace + 1
				#raise ValueError("Oh no!  Used too many tries for generating a valid period at this point.")

		
				
	#BELOW SECTION: Calculates transit values
	transitvals = calctransit(mass1=massA, massp=massB, r1=rA, r2=rB, period=periods, ecc=eccs, angle=angles, imp=imps, omega=omegas, u1=uA, u2=uB, mag=mag, numruns=numruns, ecca=ecca, eccb=eccb, band=band, periodp=periodp, poptype=poptype, radiusratio=radiusratio)
	
	#BELOW SECTION: Graphs transit values if passed-in parameter 'plot' set as 'True'
	if plot == True:
		plotmodel(transitvals, savename=savename)
	
	#Below returns the calculated values
	return transitvals
##




###BASE CALCULATIONS
##FUNCTION: calctransit
#Purpose: This function is meant to accept several inputs in cgs units (except for period, which should be given in days), including primary and orbiting masses (mass1 and massp), primary and orbiting radii (r1 and r2), a period (period), eccentricity (ecc), omega (omega), angle (angle), impact parameter (imp), and limb darkening (u1 and u2).
#It calculates analytically several aspects of the transit curve, including full width (T14), mid width (T23), transit probability (prob), exact depth (depthexact), and approximate depth (depthappr).
#Dependencies: calcsemimajor
def calctransit(mass1=None, massp=None, r1=None, r2=None, period=None, ecc=None, angle=None, imp=None, omega=None, u1=None, u2=None, mag=None, numruns=100, ecca=None, eccb=None, periodp=None, band='Kepler', poptype=None, radiusratio=None, **kwargs):
	#Below makes sure necessary parameters have been passed
	if mass1 is None or massp is None or r1 is None or r2 is None or period is None or ecc is None or angle is None or imp is None or omega is None:
		raise ValueError('Please pass in all of the following parameters with correct values using cgs units: mass1 (primary mass), massp (orbiting mass), r1 and r2 (primary and orbiting radius), period, ecc (eccentricity), angle (given as i usually), omega, and imp (impact parameter b).')
	if isinstance(mass1, float) or isinstance(mass1, int):
		if numruns != 1:
			raise ValueError("It appears that the number of runs you specified is greater than the number of values (i.e, the number of mass1s) that you passed into the function.  Please specify a valid number of runs under the parameter 'numruns'.")

	elif numruns > len(np.asarray(mass1)):
		raise ValueError("It appears that the number of runs you specified is greater than the number of values (i.e, the number of mass1s) that you passed into the function.  Please specify a valid number of runs under the parameter 'numruns'.")

	#BELOW SECTION: Calculates some constants for calculations
	#Below converts everything to numpy arrays
	if not isinstance(mass1, float) and not isinstance(mass1, int):
		mass1 = np.asarray(mass1)
		massp = np.asarray(massp)
		r1 = np.asarray(r1)
		r2 = np.asarray(r2)
		period = np.asarray(period)
		ecc = np.asarray(ecc)
		angle = np.asarray(angle)
		imp = np.asarray(imp)
		
		if omega is not None:
			omega = np.asarray(omega)
		if u1 is not None and u2 is not None:
			u1 = np.asarray(u1)
			u2 = np.asarray(u2)
			
	else: #If singular values passed in
		mass1 = np.asarray([mass1])
		massp = np.asarray([massp])
		r1 = np.asarray([r1])
		r2 = np.asarray([r2])
		period = np.asarray([period])
		ecc = np.asarray([ecc])
		angle = np.asarray([angle])
		imp = np.asarray([imp])
		
		if omega is not None:
			omega = np.asarray([omega])
		if u1 is not None and u2 is not None:
			u1 = np.asarray([u1])
			u2 = np.asarray([u2])


	#Below sets larger and smaller radii
	rs = np.zeros(len(mass1))
	rp = np.zeros(len(mass1))
	for d in range(0, len(mass1)):
		if (r1[d] > r2[d]):
			rs[d] = r1[d]
			rp[d] = r2[d]
		elif (r2[d] >= r1[d]):
			rs[d] = r2[d]
			rp[d] = r1[d]


	#Below formulates some basic values
	semi = calcsemimajor(period, (mass1+massp)) #Semimajor axis
		
		
	#BELOW SECTION: Calculates transit times, probability, and whatnot
	#Note: Formulas from Winn paper on Transits and Occultations
	#Below calculates transit times

	#For T14: in seconds
	sqrt14 = np.sqrt((1 + (rp/rs))**2.0 - imp**2.0)
	const14 = np.arcsin((rs*rsun/semi)*(sqrt14/np.sin(angle*pi/180.0)))
	approx14 = np.sqrt(1 - ecc**2.0) / (1 + ecc*np.sin(omega*pi/180.0))
	T14 = (period/pi)*const14*approx14
	
	#For T23: in days
	sqrt23 = np.sqrt((1 - (rp/rs))**2.0 - imp**2.0)
	const23 = np.arcsin((rs*rsun/semi)*(sqrt23/np.sin(angle*pi/180.0)))
	approx23 = np.sqrt(1 - ecc**2.0) / (1 + ecc*np.sin(omega*pi/180.0))
	T23 = (period/pi)*const23*approx23
	
	#Below calculates transit probability
	prob = "You didn't specify omega.  Therefore the transit probability was not generated."
	if omega is not None:
		constprr = (rs + rp)*rsun/semi
		constpre = (1 + (ecc*np.sin(omega*pi/180.0))) / (1 - ecc**2.0)
		prob = constprr * constpre
	
	#Below calculates depths
	#Approximate depth
	depthappr = (rp/rs)**2.0
	
	#Exact depth; requires u1 and u2 be specified
	#Formula given in Appendix B of Csizmadia et al paper, entitled "Transit Parameters, Stellar Spots and Limb Darkening"
	depthexact = "You didn't specify u1 and u2.  Therefore depthexact was not created."
	if u1 is not None and u2 is not None:
		constex = (1 - np.sqrt(1 - imp**2.0)) #At max light lost
		numeratorex = 1.0 - (u1*constex) - (u2*constex)**2.0
		denominatorex = 1.0 - (u1/3.0) - (u2/6.0)
		depthexact = ((rp/rs)**2.0)*numeratorex/denominatorex
	
	
	#BELOW SECTION: Trims out any 'nan' values or invalid probs
	#Below sets up arrays to hold trimmed values
	mass1done = np.zeros(numruns)
	masspdone = np.zeros(numruns)
	r1done = np.zeros(numruns)
	r2done = np.zeros(numruns)
	perioddone = np.zeros(numruns)
	semidone = np.zeros(numruns)
	eccdone = np.zeros(numruns)
	angledone = np.zeros(numruns)
	impdone = np.zeros(numruns)
	omegadone = np.zeros(numruns)
	u1done = np.zeros(numruns)
	u2done = np.zeros(numruns)
	T14done = np.zeros(numruns)
	T23done = np.zeros(numruns)
	taudone = np.zeros(numruns)
	probdone = np.zeros(numruns)
	depthapprdone = np.zeros(numruns)
	depthexactdone = np.zeros(numruns)
	magdone = np.zeros(numruns)
	
	#Below fills above arrays with non-'nan' values
	track = 0 #Track place in trimmed arrays
	placehere = 0 #Track place in original arrays
	while track < numruns:
		#Below skips current original values if 'nans' involved
		if calc.isnan(T14[placehere]) or prob[placehere] > 1.0:
			if prob[placehere] > 1.0:
				print('A probability greater than 1 was calculated...')
				print('Found at index ', placehere, '. It was discarded.')
				print('')
			placehere = placehere + 1
			
			#Below throws an error if no more data to screen through
			if placehere == len(mass1):
				print('Oh no, too many nans for given amount of data!')
				print('At a track of ', track)
				raise ValueError('Oh no!  No more data to screen through.  Seems there were too many nans or too few data to start with.')
			continue
		
		#Below adds non-'nan' values to trimmed arrays
		mass1done[track] = mass1[placehere]
		masspdone[track] = massp[placehere]
		r1done[track] = r1[placehere]
		r2done[track] = r2[placehere]
		perioddone[track] = period[placehere]
		semidone[track] = semi[placehere]
		eccdone[track] = ecc[placehere]
		angledone[track] = angle[placehere]
		impdone[track] = imp[placehere]
		omegadone[track] = omega[placehere]
		u1done[track] = u1[placehere]
		u2done[track] = u2[placehere]
		T14done[track] = T14[placehere]
		T23done[track] = T23[placehere]
		probdone[track] = prob[placehere]
		depthapprdone[track] = depthappr[placehere]
		depthexactdone[track] = depthexact[placehere]
		if mag is not None:
			magdone[track] = mag[placehere]
		
		#Below increments counts and tracking variables
		track = track + 1
		placehere = placehere + 1
		
	
	######################
	print('finished! with track as ', track, ' and placehere as ', placehere)
	
	######################
	
	#BELOW SECTION: Nan formatting for T23; also constructs tau array
	#Below replaces the nans within T23 with values of zero
	nanplaces = np.isnan(T23done)
	T23done[nanplaces] = 0.0
	
	#Below fills in tau array; tau = (T14-T23)/2
	taudone = (T14done - T23done)/2.0
	#Below replaces all taus below 30 minutes with the minimum cutoff of 30 minutes
	cutoff = 60.0*30.0 #Number of seconds in thirty minutes
	taudone[taudone < cutoff] = cutoff
	
	
	#BELOW SECTION: Returns the calculated values
	done = {'T14':T14done, 'T23':T23done, 'tau':taudone, 'prob':probdone, 'depthapprox':depthapprdone, 'depthexact':depthexactdone, 'mass1':mass1done, 'mass2':masspdone, 'r1':r1done, 'r2':r2done, 'period':perioddone, 'semi':semidone, 'ecc':eccdone, 'angle':angledone, 'imp':impdone, 'omega':omegadone, 'u1':u1done, 'u2':u2done, 'ecca':ecca, 'eccb':eccb, 'periodp':periodp, 'numruns':numruns, 'mag':mag, 'band':band, 'poptype':poptype, 'radiusratio':radiusratio}

	#Below returns generated values
	return done
##


##FUNCTION: calcsemimajor
#Purpose: This function calculates the semimajor axis (in centimeters) using a given star mass (in grams) and period (in days)
def calcsemimajor(period, massstar):
	#Semimajor axis formula, by Kepler's Law:
	#T^2/a^3 = (4*pi^2 / (G*Mstar))
	
	#Below calculates the semimajor axis
	const = (Gconst*massstar*masssun)/(4.0*pi**2.0)
	semimajoraxis = (const*period**2.0)**(1.0/3.0)

	#Below returns calculated value
	return semimajoraxis #In centimeters
##


##FUNCTION: calcomega; WRONG?!? OUTDATED
#Purpose: This function calculates omega based upon the given semimajor axis (semi), impact parameter (imp), eccentricity (ecc), primary radius (rs), and angle (angle).
def calcomega(semi, imp, ecc, rs, angle, shift=False):
	#Omega formula derived from Winn paper on Transits and Occultations
	frac1 = (semi*np.cos(angle*pi/180.0)) * (1.0 - ecc**2.0) / (rs*imp*rsun)
	frac2 = (frac1 - 1.0) / ecc
	omega = np.arcsin(frac2)*180.0/pi
	
	#BELOW SECTION: Phase shifts omega into 0 to 360 degree frame, if 'shift' passed in as 'True'
	if shift == True:
		#Below phase shifts omega into the 0 to 180 degree frame, if needed
		import random as rand
		if -180 < omega < 0:
			omega = omega + 180
		#Below next randomly assigns value to 0-t0-180 frame or 180-to-360 frame
		randnum = rand.random()
		if randnum >= 0.5:
			omega = omega + 180
		
	#Below returns calculated omega
	return omega
##


##FUNCTION: calcangle
#Purpose: This function calculates the angle i for a given orbit using already determined parameters
def calcangle(period, semi, r1, omega, imp, ecc):
	#Below calculates angle from impact parameter formula given in Winn paper
	val1 = (imp*rsun*r1)/semi
	val2 = (1 + ecc*np.sin(omega))/(1.0 - ecc**2.0)
	anglehere = np.arccos(val1*val2)*180.0/(1.0*pi)
	
	#Below returns calculated angle
	return anglehere
##


##FUNCTION: calcminangle
#Purpose: This function calculates the minimum allowed angle i for a given orbit
def calcminangle(period, mass1, mass2, r1, r2):
	#Below calculates minimum allowed angle i
	semi = calcsemimajor(period, (mass1+mass2))
	val = ((r1+r2)*rsun) / semi
	minangle = np.arccos(val)*180.0/pi
	
	#Below returns calculated angle
	return minangle
##


##FUNCTION: calcimp
#Purpose: This function calculates the impact parameter (usually given as b) analytically using a given angle (usually denoted by i), eccentricity (ecc), semimajor axis (semi), primary radius (rs), and omega (omega).
def calcimp(ecc, angle, semi, rs, omega):
	#Impact parameter formula derived from Winn paper on Transits and Occultations
	frac1 = (semi*np.cos(angle*pi/180.0)) / (rs * rsun)
	frac2 = (1 - ecc**2.0)
	bot1 = (1.0 + (ecc*np.sin(omega*pi/180.0)))

	#Below returns the calculated impact parameter
	imp = frac1*frac2/bot1
	return imp
##


##FUNCTION: calcimpactsemi
#Purpose: Below as a method to calculate semimajor axis from given impact parameter values.
def calcimpactsemi(ecc, angle, rs, omega, imp):
	#Formula derived from Winn paper on Transits and Occultations
	part1 = (imp*rs*rsun) / np.cos(angle*pi/180.0)
	part2 = (1 - ecc**2.0) / (1 + ecc*np.sin(omega*pi/180.0))
	
	#Returns calculated impact-parameter-based semimajor axis
	return (part1/part2)
	

##FUNCTION: calcperiod
#Purpose: Below as a method to calculate period based on Kepler's Third Law.
def calcperiod(semi, massstar):
	#Calculates period from Kepler's 3rd Law
	part1 = 4*pi**2.0/(Gconst*massstar*masssun)
	part2 = part1*semi**3.0
	
	#Returns calculated period
	return (abs(part2))**(1/2.0)
	
	
##FUNCTION: makeimp
#Purpose: Below as a method to return a given number (n) of impact parameters, sampled uniformly between 0 and (rs + rp)/rs, where rs = star radius, rp = orbiting radius.
def makeimp(r1, r2):
	#Below imports basic packages
	import random as rand
	
	#Below sets the primary and orbiting radii based on which one is larger
	if isinstance(r1, int):
		r1 = np.array([r1])
		r2 = np.array([r2])
	rs = np.zeros(len(r1))
	rp = np.zeros(len(r1))

	for f in range(0, len(r1)):
		if (r1[f] > r2[f]):
			rs[f] = r1[f]
			rp[f] = r2[f]
		elif (r2[f] >= r1[f]):
			rs[f] = r2[f]
			rp[f] = r1[f]
	
	#Below samples random numbers between 0 and endrange
	endrange = (rs + rp)/rs
	randdist = np.zeros(len(r1))
	for b in range(0, len(r1)):
		randdist[b] = rand.random() * endrange[b]
	
	#Below returns uniform unscaled impact parameter results
	return randdist
##


##FUNCTION: makeangle
#Purpose: Below as a method to return a given number (n) of angles, sampled uniformly between 0 and 360 degrees.
def makeangle(n):
	#Below imports basic packages
	import random as rand
	n = int(n)
	
	#Below samples random numbers between 0 and endrange
	endrange = 360
	randdist = np.zeros(n)
	for c in range(0, n):
		randdist[c] = rand.random() * endrange
	
	#Below generates uniform angle results
	return randdist
##

	
##FUNCTION: makeperiod
#Purpose: Below as method to return a given number (n) of periods, sampled from a Power Law Distribution, to the power of p; rangestart and rangeend should be passed through in days
#Notes: If dayunits passed as False, will work in units of seconds
def makeperiod(n, p, rangestart=None, rangeend=None, dayunits=True):
	#Below imports basic packages
	import random as rand
	n = int(n)
	
	#BELOW SECTION: Samples random numbers in a certain range
	#If no range given, then samples between zero and endrange generated below
	#Throws an error if only one range endpoint given
	if rangestart is None and rangeend is not None or rangestart is not None and rangeend is None:
		raise ValueError("Oh no!  You only specified one range endpoint.  Please pass in both parameters 'rangestart' and 'rangeend' simultaneously, or do not give a range at all.")
		
	if rangestart is None and rangeend is None:
		endrange = 10 #Gives the cutoff for random sampling, so random numbers drawn from sample [0, endrange)
		randnum = np.zeros(n) #To hold random numbers
		for a in range(0, n):
			randnum[a] = rand.random() * endrange
	
		#Below generates the power law results
		if dayunits is True: #In units of days
			perdone = (randnum**p)
			return perdone
		
		elif dayunits is False: #In units of seconds
			perdone = (randnum**p)*daysecs
			return perdone
		
	#If range is given, then samples from between endrange
	elif rangestart is not None and rangeend is not None:
		randinrange = np.zeros(n) #To hold random numbers
		for b in range(0, n):
			randinrange[b] = rand.uniform((rangestart*daysecs)**(1.0/p), (rangeend*daysecs)**(1.0/p))
			
		#Below generates the power law results
		perdone = (randinrange**p)
		return perdone
##


##FUNCTION: makeecc
#Purpose: Below as method to return a given number of eccentricities, sampled from a Beta Distribution
def makeecc(n, a, b):
	#Below imports basic packages
	from scipy.stats import beta
	n = int(n)
	
	#Below returns an array of values sampled from beta distribution with given parameters (a, b) and of size n
	eccdone = beta.rvs(a, b, size = n)
	return eccdone
##
	
	
	

###UTILITIES
##UTILITY: rochelobe; borrowed from T. Morton's work; thanks!
def rochelobe(q):
    return 0.49*q**(2./3)/(0.6*q**(2./3) + np.log(1+q**(1./3)))
## #End of borrowed code


##UTILITY: Checks if values within roche limit
def withinroche(semi, mass1, mass2, r1, r2):
	ratio = mass1/(mass2*1.0)
	bool = ((r1 + r2)*rsun) > (rochelobe(ratio)*semi)

	return bool
##

	

	
##########################
###TESTS
##TEST: testcalctransit; OUTDATED TEST
#Purpose: This test is meant to test the calctransit function
def testcalctransit(filename, numruns=100, **kwargs):
	#Below calls an error if given file as nonexistent
	import os.path
	if os.path.isfile(filename) is False:
		raise ValueError("Please pass in an existent file name.")

	#Below generates sample star population using vespa.populations BEBPopulation to read in the file
	from vespa.populations import BEBPopulation as beb
	origpop = beb(trilegal_filename=filename)
	
	#Below calculates a few values
	semi = calcsemimajor(period=origpop.stars.P, massstar=origpop.stars.mass_1+origpop.stars.mass_2) #Semimajor axis; add primary and orbiting masses?
	imps = calcimp(angle=origpop.stars.inc, ecc=origpop.stars.ecc, omega=origpop.stars.w, rs=origpop.stars.radius_1, semi=semi) #Impact parameter
	
	#BELOW SECTION: Feeds in calculated periods and other values into calctransit function, to test to make sure that the transit calculations are identical
	testvals = calctransit(mass1=origpop.stars.mass_1, massp=origpop.stars.mass_2, r1=origpop.stars.radius_1, r2=origpop.stars.radius_2, period=origpop.stars.P, ecc=origpop.stars.ecc, angle=origpop.stars.inc, imp=imps, omega=origpop.stars.w, u1=origpop.stars.u1_1, u2=origpop.stars.u2_1, numruns=numruns)
	
	#Below prints the differences between the first handful of transit calculations
	origT14 = np.array(origpop.stars.T14_pri)
	origT23 = np.array(origpop.stars.T23_pri)
	testT14 = testvals['T14']
	testT23 = testvals['T23']
	
	#Below tests the output
	testfile = open('testfile.txt', 'w')
	testfile.write('Testing.\n')
	testfile.write('First Original T14 Values as:\n')
	g = 0
	index = 0
	while index < 100:
		if np.sort(origT14)[g] == 0.0:
			g = g + 1
			continue
		testfile.write(str(np.sort(origT14)[g]) + '\n')
		g = g + 1
		index = index + 1
	testfile.write('\n\n')
	testfile.write('First Testing T14 Values as:\n')
	for g in range(0, 100):
		testfile.write(str(np.sort(testT14)[g]) + '\n')
	testfile.write('\n\n\n\n')
	testfile.write('First Original T23 Values as:\n')
	g = 0
	index = 0
	while index < 100:
		if np.sort(origT23)[g] == 0.0:
			g = g + 1
			continue
		testfile.write(str(np.sort(origT23)[g]) + '\n')
		g = g + 1
		index = index + 1
	testfile.write('\n\n')
	testfile.write('First Testing T23 Values as:\n')
	for g in range(0, 100):
		testfile.write(str(np.sort(testT23)[g]) + '\n')
	testfile.write('\n\n')	
##
