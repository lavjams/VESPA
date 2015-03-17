#####
#File: transitfunctions.py
#Purpose: Contains functions that relate to the transit curves of transits.
#:-) Started 3-14-15
#
#####

#Below imports necessary packages
import numpy as np
import astropy.constants as const

#Below imports necessary constants
#Below are some constant values
daysecs = 86400 #Number of seconds in a day
masssun = const.M_sun.cgs.value #Mass of Sun in grams
Gconst = const.G.cgs.value #Gravitational constant, cgs
rsun = const.R_sun.cgs.value #Radius of sun, cgs
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
def runmodel(filename, periodp=2, ecca=0.5, eccb=0.5, numruns=1e3, **kwargs):
	#Below calls an error if given file as nonexistent
	import os.path
	if os.path.isfile(filename) is False:
		raise ValueError("Please pass in an existent file name.")
	
	#Below opens up the star distribution
	pop = gensim(filename, n=numruns)
	
	#Below generates random parameter distributions
	periods = makeperiod(numruns, periodp) #Power law distribution
	eccs = makeecc(numruns, ecca, eccb) #Beta distribution
	angles = makeangle(numruns) #Uniform distribution
	imps = makeimp(numruns, pop.stars.radius_A, pop.stars.radius_B) #Uniform distribution
	semi = calcsemimajor(period=periods, massstar=pop.stars.mass_A)
	omegas = calcomega(semi=semi, imp=imps, angle=angles, rs=pop.stars.radius_A, ecc=eccs) #Uniform distribution
	
	#Below calculates transit values
	transitvals = calctransit(mass1=pop.stars.mass_A, massp=pop.stars.mass_B, r1=pop.stars.radius_A, r2=pop.stars.radius_B, period=periods, ecc=eccs, angle=angles, imp=imps, omega=omegas, u1=pop.stars.Teff_A, u2=pop.stars.logg_A, numruns=numruns)
	
	#Below returns the calculated values
	return transitvals
##


##FUNCTION: calctransit
#Purpose: This function is meant to accept several inputs in cgs units (except for period, which should be given in days), including primary and orbiting masses (mass1 and massp), primary and orbiting radii (r1 and r2), a period (period), eccentricity (ecc), omega (omega), angle (angle), impact parameter (imp), and limb darkening (u1 and u2).
#It calculates analytically several aspects of the transit curve, including full width (T14), mid width (T23), transit probability (prob), exact depth (depthexact), and approximate depth (depthappr).
#Dependencies: calcsemimajor
def calctransit(mass1=None, massp=None, r1=None, r2=None, period=None, ecc=None, angle=None, imp=None, omega=None, u1=None, u2=None, numruns=1e3, **kwargs):
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
	semi = calcsemimajor(period, mass1) #Semimajor axis
		
		
	#BELOW SECTION: Calculates transit times, probability, and whatnot
	#Note: Formulas from Winn paper on Transits and Occultations
	#Below calculates transit times

	#For T14:
	sqrt14 = np.sqrt((1 + (rp/rs))**2.0 - imp**2.0)
	const14 = np.arcsin((rs*rsun/semi)*(sqrt14/np.sin(angle*pi/180.0)))
	T14 = (period/pi)*const14
	
	#For T140: This is the version used in T. Morton's code; not sure why yet...
	T140 = (period/pi)*const14 *\
		np.sqrt(1-ecc**2.0)/(1+ecc*np.sin(omega*pi/180.0))

	#For T23:
	sqrt23 = np.sqrt((1 - (rp/rs))**2.0 - imp**2.0)
	const23 = np.arcsin((rs*rsun/semi)*(sqrt23/np.sin(angle*pi/180.0)))
	T23 = (period/pi)*const23
	
	#Below calculates transit probability
	prob = "You didn't specify omega.  Therefore the transit probability was not generated."
	if omega is not None:
		constprr = (rs + rp)/semi
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
	done = {'T14':T14, 'T23':T23, 'T140':T140, 'prob':prob, 'depthappr':depthappr, 'depthexact':depthexact}

	#Below returns generated values
	return done
##



###TESTS
##TEST: testcalctransit
#Purpose: This test is meant to test the calctransit function
def testcalctransit(filename, numruns=1e3, **kwargs):
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
	testT140 = testvals['T140']
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
	testfile.write('\n\n')
	testfile.write('First Testing T140 Values as:\n')
	for g in range(0, 100):
		testfile.write(str(np.sort(testT140)[g]) + '\n')
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


###BASE CALCULATIONS
##FUNCTION: calcsemimajor
#Purpose: This function calculates the semimajor axis (in centimeters) using a given star mass (in grams) and period (in days)
def calcsemimajor(period, massstar):
	#Below converts the period from days to seconds
	periodsecs = period*daysecs
	
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
	frac1 = (rs*imp) / (semi*np.cos(angle*pi/180.0)) / (1.0 - ecc**2.0)
	frac2 = (frac1 - 1.0) / ecc
	omega = np.arcsin(frac2)
	
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


##FUNCTION: makeimp
#Purpose: Below as a method to return a given number (n) of impact parameters, sampled uniformly between 0 and (rs + rp)/rs, where rs = star radius, rp = orbiting radius.
def makeimp(n, r1, r2):
	#Below imports basic packages
	import random as rand
	n = int(n)
	
	#Below sets the primary and orbiting radii based on which one is larger
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
	randdist = np.zeros(n)
	for b in range(0, n):
		randdist[b] = rand.random() * endrange[b]
	
	#Below generates uniform impact parameter results
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
def makeperiod(n, p):
	#Below imports basic packages
	import random as rand
	n = int(n)
	
	#Below samples random numbers between 0 and endrange
	endrange = 10 #Gives the cutoff for random sampling, so random numbers drawn from sample [0, endrange)
	randnum = np.zeros(n) #To hold random numbers
	for a in range(0, n):
		randnum[a] = rand.random() * endrange
	
	#Below generates the power law results
	perdone = randnum**p
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
	
	
	

	
	
	