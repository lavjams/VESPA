#####
#File: transitfunctions.py
#Purpose: Contains functions that relate to the transit curves of transits.
#:-) Started 3-14-15
#
#####

#Below imports necessary packages
import numpy as np
import math as calc
import matplotlib.pyplot as graph
import astropy.constants as const

#Below imports necessary constants
#Below are some constant values
daysecs = 86400.0 #Number of seconds in a day
masssun = const.M_sun.cgs.value #Mass of Sun in grams
massearth = const.M_earth.cgs.value #Mass of Earth in grams
Gconst = const.G.cgs.value #Gravitational constant, cgs
rsun = const.R_sun.cgs.value #Radius of Sun, cgs
eccearth = 0.01671123 #Eccentricty of Earth's orbit; ###FIND EXACT VALUE!
au = const.au.cgs.value #Astronomical unit, cgs
pi = np.pi #Value of pi


##FUNCTION: gensim ###NOT TESTED
#Purpose: This function is meant to run a star simulation to generate a population of stars and statistics
#NOTE: Sections of code heavily (very heavily) adapted from Tim Morton's vespa.populations code; thanks!
#Dependencies: starutils.populations.BGStarPopulation_TRILEGAL
def gensim(trilegal_filename, ra=None, dec=None, n=2e4, ichrone=None, MAfn=None, mags=None, maxrad=10, f_binary=0.4, **kwargs):
	#Below imports necessary functions from elsewhere
	from starutils.populations import BGStarPopulation_TRILEGAL
	from starutils.populations import MultipleStarPopulation
	if ichrone is None:
		from starutils.populations import DARTMOUTH
		ichrone = DARTMOUTH
		
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


##FUNCTION: runmodel
#Purpose: This function is meant to run a star simulation to generate a population of stars and statistics, and then perform transit calculations.  It accepts one value to generate period distributions (power law distribution to power periodp) and two values to generate eccentricity distributions (eccentricity distribution with alpha ecca and beta eccb).
#Numruns indicates the number of stars.
#Dependencies:
def runmodel(filename='samplepop.h5', periodp=2, ecca=2, eccb=2, numruns=100, plot=False, **kwargs):
	#Below calls an error if given file as nonexistent
	import os.path
	if os.path.isfile(filename) is False:
		raise ValueError("Please pass in an existent file name.")
	
	#BELOW SECTION: Opens up the star distribution and records necessary values
#	totalruns = numruns*100
	pop = gensim(filename)#, n=totalruns)
	#Below are star parameters
	rAraw = pop.stars.radius_A#[0:totalruns]
	rBraw = pop.stars.radius_B#[0:totalruns]
	massAraw = pop.stars.mass_A#[0:totalruns]
	massBraw = pop.stars.mass_B#[0:totalruns]
	uAraw = pop.stars.Teff_A#[0:totalruns]
	uBraw = pop.stars.Teff_B#[0:totalruns]
	totalruns = len(rAraw)

	#BELOW SECTION: Shuffles generated star parameters
	#Below generates random indices
#	indices = np.arange(numruns) #Array of indices
#	np.random.shuffle(indices) #Shuffles indices
	
	#Below creates blank arrays to hold shuffled star parameters
#	rA = np.zeros(numruns)
#	rB = np.zeros(numruns)
#	massA = np.zeros(numruns)
#	massB = np.zeros(numruns)
#	uA = np.zeros(numruns)
#	uB = np.zeros(numruns)
	
	#Below fills the blank arrays with shuffled values
#	for a in range(0, numruns):
#		rA[a] = rAraw[indices[a]]
#		rB[a] = rBraw[indices[a]]
#		massA[a] = massAraw[indices[a]]
#		massB[a] = massBraw[indices[a]]
#		uA[a] = uAraw[indices[a]]
#		uB[a] = uBraw[indices[a]]
	
	#BELOW SECTION: Generates random parameter distributions
	eccs = makeecc(numruns, ecca, eccb) #Beta distribution

	#Below generates angles and imps
	angles = makeangle(numruns)
	impsraw = makeimp(rAraw, rBraw)
	
	#Below calculates omegas; throws out any invalid omega values
	omegas = np.zeros(numruns)
	rA = np.zeros(numruns)
	rB = np.zeros(numruns)
	massA = np.zeros(numruns)
	massB = np.zeros(numruns)
	uA = np.zeros(numruns)
	uB = np.zeros(numruns)
	periods = np.zeros(numruns)
	imps = np.zeros(numruns)
	semitesting = np.zeros(numruns) #For TESTING purposes
	
	track = 0 #Track number of valid calculations
	tryhere = 0 #Track number of tried periods
	totalplace = 0 #Track place in original total arrays
	while track < numruns:
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
		
		
		
		#######TRY to determine range of allowed period values
		#Below determines range of allowed period values
		highestsemi = calcimpactsemi(imp=impsraw[totalplace], angle=angles[track], rs=rhere, ecc=eccs[track], omega=360.0)
		lowestsemi = calcimpactsemi(imp=impsraw[totalplace], angle=angles[track], rs=rhere, ecc=eccs[track], omega=270.0)

		highestperiod = calcperiod(semi=highestsemi, massstar=(massAraw[totalplace]+massBraw[totalplace]))
		lowestperiod = calcperiod(semi=lowestsemi, massstar=(massAraw[totalplace]+massBraw[totalplace]))
	
		#Below generates period here
		periodhere = makeperiod(1, periodp, rangestart=lowestperiod, rangeend=highestperiod)#*daysecs
		#Power law distribution, converted to seconds
		#Below computes omega with current values
		semi = calcsemimajor(period=periodhere, massstar=massAraw[totalplace]+massBraw[totalplace])
		if highestsemi < 0: #Temp fix for sign issues of semi
			semi = semi*(-1)
		
		omegahere = calcomega(semi=semi, imp=impsraw[totalplace], angle=angles[track], rs=rhere, ecc=eccs[track])
		
		
		###################
#		print('')
#		print('currently, highest period as ', highestperiod)
#		print('lowest period as ', lowestperiod)
#		print('highest semi as ', highestsemi)
#		print('lowest semi as ', lowestsemi)
#		print('period here is ', periodhere)
#		print('semi as ', semi)
#		print('omegahere as ', omegahere)
		###################

		################################## TEMP BORROW
		"""
		tooclose = withinroche(semi/au*(1-eccs[track]),massA[track],rA[track],massB[track],rB[track])
		ntooclose = np.array(tooclose).sum()
		tries = 0
		maxtries = 5
		while ntooclose > 0:
			lastntooclose=ntooclose
			periodhere = makeperiod(1, periodp)*daysecs

			semimajors = calcsemimajor(periodhere, (massA[track]+massB[track]))
			
			tooclose = withinroche(semimajors/au*(1-eccs[track]),massA[track],rA[track],massB[track],rB[track])
			
			ntooclose = np.array(tooclose).sum()
			if ntooclose==lastntooclose:   #prevent infinite loop
				tries += 1
				if tries > maxtries:
					raise ValueError("Too many withinroche!")
		"""			                       

		########################## END OF TEMP BORROW
		
		
		
		if not calc.isnan(omegahere):
			
		#Below keeps data if valid omega parameter was produced; throws out otherwise
		
#		if calc.isnan(omegahere): #If not valid omega produced
#			placehere = placehere + 1
#			if placehere >= totalruns and track < (numruns - 1):
#				print('track as ', track) ###########
#				raise ValueError("Oh no!  Ran out of totalrun values to choose from.  Better make totalrun values larger, perhaps.")
#		else: #If valid omega produced
			omegas[track] = omegahere
			semitesting[track] = semi
			imps[track] = impsraw[totalplace]
			rA[track] = rAraw[totalplace]
			rB[track] = rBraw[totalplace]
			massA[track] = massAraw[totalplace]
			massB[track] = massBraw[totalplace]
			uA[track] = uAraw[totalplace]
			uB[track] = uBraw[totalplace]
			periods[track] = periodhere/daysecs
			
			print("Success!  At track ", track, ' and tryhere ', tryhere)
			print("Also with totalplace of ", totalplace)
			print("With period here of ", periodhere)
			track = track + 1
			totalplace = totalplace + 1
			tryhere = 0
		else:
			tryhere = tryhere + 1
			if tryhere >= 10000:
				print("Oh no!  Tryhere was reached.")
				print("Current values:")
				print('periodhere as ', periodhere)
				print('semi here as ', semi)
				print('highestperiod as ', highestperiod)
				print('lowestperiod as ', lowestperiod)
				print('highestsemi as ', highestsemi)
				print('lowestsemi as ', lowestsemi)
				print('ra here as ', rhere)
				print('rb here as ', rnothere)
				print('massa here as ', massAraw[totalplace])
				print('massb here as ', massBraw[totalplace])
				
				print('')
				print('tryhere reached ', tryhere)
				print('totalplace reached ', totalplace)
				print('Current track as ', track)
				raise ValueError("Oh no!  Used too many tries for generating a valid period at this point.")

			#Below increments counts and tracking variables
#			track = track + 1
#			placehere = placehere + 1
#			if placehere >= totalruns and track < (numruns - 1):
#				print('track as ', track) #############
#				raise ValueError("Oh no!  Ran out of totalrun values to choose from.  Better make totalrun values larger, perhaps.")

			
	
	print('new imps as ', imps) ############
	print('new omegas as ', omegas) ###############
	print('mass1*masssun as ', massA*masssun) #############
	print('r1*rsun as ', rA*rsun) ###############
	print('r2*rsun as ', rB*rsun) #################
	print('semi as ', semitesting) ##############
	print('periods as ', periods/daysecs) ###########
		
	#BELOW SECTION: Calculates transit values
	transitvals = calctransit(mass1=massA, massp=massB, r1=rA, r2=rB, period=periods, ecc=eccs, angle=angles, imp=imps, omega=omegas, u1=uA, u2=uB, numruns=numruns)
	
	#BELOW SECTION: Graphs transit values if passed-in parameter 'plot' set as 'True'
	if plot == True:
		graph.scatter(transitvals['T14'], transitvals['depthapprox'], s=4)
		graph.show()
		
		graph.scatter(transitvals['T23'], transitvals['depthapprox'], s=4)
		graph.show()
	
	#Below returns the calculated values
	return transitvals
##




###BASE CALCULATIONS
##FUNCTION: calctransit
#Purpose: This function is meant to accept several inputs in cgs units (except for period, which should be given in days), including primary and orbiting masses (mass1 and massp), primary and orbiting radii (r1 and r2), a period (period), eccentricity (ecc), omega (omega), angle (angle), impact parameter (imp), and limb darkening (u1 and u2).
#It calculates analytically several aspects of the transit curve, including full width (T14), mid width (T23), transit probability (prob), exact depth (depthexact), and approximate depth (depthappr).
#Dependencies: calcsemimajor
def calctransit(mass1=None, massp=None, r1=None, r2=None, period=None, ecc=None, angle=None, imp=None, omega=None, u1=None, u2=None, numruns=100, **kwargs):
	#Below makes sure necessary parameters have been passed
	if mass1 is None or massp is None or r1 is None or r2 is None or period is None or ecc is None or angle is None or imp is None:
		raise ValueError('Please pass in all of the following parameters with correct values using cgs units (except for period, which should be given in days): mass1 (primary mass), massp (orbiting mass), r1 and r2 (primary and orbiting radius), period (the period), ecc (eccentricity), angle (given as i usually), and imp (impact parameter b).')
	if numruns > len(mass1):
		raise ValueError("It appears that the number of runs you specified is greater than the number of values (i.e, the number of mass1s) that you passed into the function.  Please specify a valid number of runs under the parameter 'numruns'.")

	#BELOW SECTION: Calculates some constants for calculations
	#Below converts everything to numpy arrays
	mass1 = np.array(mass1)
	massp = np.array(massp)
	r1 = np.array(r1)
	r2 = np.array(r2)
	period = np.array(period)
	ecc = np.array(ecc)
	angle = np.array(angle)
	imp = np.array(imp)
	
	if omega is not None:
		omega = np.array(omega)
	if u1 is not None and u2 is not None:
		u1 = np.array(u1)
		u2 = np.array(u2)
	
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
	
	#BELOW SECTION: Trims values to fit specified number of runs
	if (len(mass1) > numruns):
		mass1 = mass1[0:numruns] #Primary mass
	if (len(massp) > numruns):
		massp = massp[0:numruns] #Orbiting mass
	if (len(rs) > numruns):
		rs = rs[0:numruns] #Primary radius
	if (len(rp) > numruns):
		rp = rp[0:numruns] #Orbiting radius
	if (len(period) > numruns):
		period = period[0:numruns] #Period
	if (len(ecc) > numruns):
		ecc = ecc[0:numruns] #Eccentricity
	if (len(angle) > numruns):
		angle = angle[0:numruns] #Inclination
	if (len(imp) > numruns):
		imp = imp[0:numruns] #Impact parameter
	if (len(omega) > numruns):
		omega = omega[0:numruns] #Omega
	
	if u1 is not None and u2 is not None: #Darkening parameters
		if (len(u1) > numruns):
			u1 = u1[0:numruns]
		if (len(u2) > numruns):
			u2 = u2[0:numruns]

	#Below formulates some basic values
	semi = calcsemimajor(period, (mass1+massp)) #Semimajor axis
		
		
	#BELOW SECTION: Calculates transit times, probability, and whatnot
	#Note: Formulas from Winn paper on Transits and Occultations
	#Below calculates transit times

	#For T14:
	sqrt14 = np.sqrt((1 + (rp/rs))**2.0 - imp**2.0)
	const14 = np.arcsin((rs*rsun/semi)*(sqrt14/np.sin(angle*pi/180.0)))
	approx14 = np.sqrt(1 - ecc**2.0) / (1 + ecc*np.sin(omega*pi/180.0))
	T14 = (period/pi)*const14*approx14
	
	#For T23:
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
	
	#BELOW SECTION: Returns the calculated values
	done = {'T14':T14, 'T23':T23, 'prob':prob, 'depthapprox':depthappr, 'depthexact':depthexact}

	#Below returns generated values
	return done
##


##FUNCTION: calcsemimajor
#Purpose: This function calculates the semimajor axis (in centimeters) using a given star mass (in grams) and period (in days)
def calcsemimajor(period, massstar):
	#Below converts the period from days to seconds
	periodsecs = period#*daysecs
	
	#Semimajor axis formula, by Kepler's Law:
	#T^2/a^3 = (4*pi^2 / (G*Mstar))
	
	#Below calculates the semimajor axis
	const = (Gconst*massstar*masssun)/(4.0*pi**2.0)
	semimajoraxis = (const*periodsecs**2.0)**(1.0/3.0)

	#Below returns calculated value
	return semimajoraxis #In centimeters
##


##FUNCTION: calcomega
#Purpose: This function calculates omega based upon the given semimajor axis (semi), impact parameter (imp), eccentricity (ecc), primary radius (rs), and angle (angle).
def calcomega(semi, imp, ecc, rs, angle):
	#Omega formula derived from Winn paper on Transits and Occultations
	frac1 = (semi*np.cos(angle*pi/180.0)) * (1.0 - ecc**2.0) / (rs*imp*rsun)
#	print('Calculating omega!') ##############
#	print('frac1 part1 as ', (semi*np.cos(angle*pi/180.0))) ########
#	print('frac1 as ', frac1) ############
#	print('') #############
	frac2 = (frac1 - 1.0) / ecc
#	print('rs*rsun as ', rs*rsun) ###############
#	print('')
#	print('semi as ', semi) #################
#	print('')
#	print('frac2 as ', frac2) ####################
	omega = np.arcsin(frac2)*180.0/pi
	
	#################
#	print('')
#	print('frac2 as ', frac2)
#	print('')
#	print('omega as ', omega) #####
	#################
	#Below returns calculated omega
	return omega
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
	
	return (part1/part2)
	

##FUNCTION: calcperiod
#Purpose: Below as a method to calculate period based on Kepler's Third Law.
def calcperiod(semi, massstar):
	part1 = 4*pi**2.0/(Gconst*massstar*masssun)
	part2 = part1*semi**3.0
		
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
#Purpose: Below as method to return a given number (n) of periods, sampled from a Power Law Distribution, to the power of p
def makeperiod(n, p, rangestart=None, rangeend=None):
	#Below imports basic packages
	import random as rand
	n = int(n)
	
	#BELOW SECTION: Samples random numbers in a certain range
	#If no range given, then samples between zero and endrange generated below
	if rangestart is None or rangeend is None:
		endrange = 10 #Gives the cutoff for random sampling, so random numbers drawn from sample [0, endrange)
		randnum = np.zeros(n) #To hold random numbers
		for a in range(0, n):
			randnum[a] = rand.random() * endrange
	
		#Below generates the power law results
		perdone = randnum**p
		return perdone
		
	#If range is given, then samples from between endrange
	elif rangestart is not None and rangeend is not None:
		randinrange = np.zeros(n) #To hold random numbers
		for b in range(0, n):
			randinrange[b] = rand.uniform(rangestart**(1.0/p), rangeend**(1.0/p))
			
		#Below generates the power law results
		perdone = randinrange**p
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
	
	
	

###TESTS
##TEST: testcalctransit
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

########################
def rochelobe(q):
    """returns r1/a; q = M1/M2"""
    return 0.49*q**(2./3)/(0.6*q**(2./3) + np.log(1+q**(1./3)))

def withinroche(semimajors,M1,R1,M2,R2):
    q = M1/M2
    return ((R1+R2)*rsun) > (rochelobe(q)*semimajors)

	