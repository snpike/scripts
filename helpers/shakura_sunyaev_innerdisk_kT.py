import numpy as np

### A bunch of equations about accretion disks from Shakura & Sunyaev 1973 
### (Astronomy and Astrophysics, Vol. 24, p. 337 - 355)

def surface_temp(L38, m, R, alpha=1.0):
	### Return the surface temperature in Kelvin of the inner disk (eq. 3.7)
	# L38: unabsorbed bolometric luminosity in units of 10^38erg/s
	# m: mass of the accreting object in units of solar mass
	# R: inner disk radius in kilometers
	r = R/(9.*m)
	Mdot_crit = Mdot_crit_from_m(m)
	Mdot = Mdot_from_L38(L38)
	return (1.4e9*np.power(alpha, 2./9.)*np.power((Mdot/Mdot_crit), 8./9.)*np.power(m, -2./9.)*\
		np.power(r,-5./3.)*np.power(1.-np.power(r,-1./2.), 8./9.))


def Mdot_from_L38(L38, eta=0.1):
	### Return mass accretion rate in g/s from the luminosity
	# L38: unabsorbed bolometric luminosity in units of 10^38erg/s
	# eta: efficiency of radiation
	return (L38/(eta*np.square(3e10)))*1e38

def Mdot_crit_from_m(m, eta=0.1):
	### Return critical mass accretion rate in g/s (page 339, top left)
	# m: mass of the accretor in units of solar mass
	# eta: efficiency of radiation
	return (m/eta) * 1.14e17