import numpy as np 

# Obtained by finding the start time in MJD, MJD-OBS, using fkeyprint and adding the elapsed time, TELAPSE
obs10002013001 = [5.611328497898149E+04, 5.611431833546297E+04]
obs10002013003 = [5.614510150675926E+04, 5.614545923824074E+04]

### From Falanga et al. 2015 (Table 3)
T0 = 52846.688810
T0err = 0.000024

Porb = 3.891923160
Porberr = 0.000000066

PorbDot = (-3.541E-6/365.0) * Porb
PorbDoterr = PorbDot*np.sqrt(np.square(0.002/3.541) + np.square(Porberr/Porb))
###

Tn = []
Tnerr = []

temp = T0
n = 0
while temp < 57690.00000000:
	# Determine the time of eclipse using Equation 1 from Falanga et al
	temp = T0 + (n*Porb) + (np.square(n)*Porb*PorbDot/2.0)
	Tn.append(temp)
	# Determine error in eclipse time Tn by propagating the standard deviation errors in T0, Porb, and PorbDot
	Tnerr.append(np.sqrt(np.square(T0err) + np.square((n+(np.square(n)*PorbDot/2.0))*Porberr) + np.square(np.square(n)*PorbDot*PorbDoterr/2.0)))

	if obs10002013001[0] < temp < obs10002013001[1]:
		print('Eclipse in observation 10002013001')
	if obs10002013003[0] < temp < obs10002013003[1]:
		print('Eclipse in observation 10002013003')

	n += 1

#print(Tn)

# Eclipses immediately preceding and following observation 10002013001
eclipses = [56111.960621598875, 56115.852421396623]
eclipseIndices = []
for i in range(len(Tn)):
	if Tn[i]==56111.960621598875:
		eclipseIndices.append(i)
	if Tn[i]==56115.852421396623:
		eclipseIndices.append(i)
#print(eclipseIndices)


# Returns the phase given time in MJD
def phase(x):
	#print(eclipses)
	return (x/(eclipses[1]-eclipses[0])) - (eclipses[0]/(eclipses[1]-eclipses[0]))

# Returns the phase error by propagating standard deviation errors in Tn
def phaseErr(x):
	return ((x-eclipses[0])*np.sqrt(np.square(Tnerr[eclipseIndices[0]]/(eclipses[1]-eclipses[0])) + np.square(Tnerr[eclipseIndices[1]]))/np.square(eclipses[1]-eclipses[0]))

print('\n')
print('Obs 10002013001:')
print('Start: ' + str(phase(obs10002013001[0])) + ' +- ' + str(phaseErr(obs10002013001[0])))
print('End: ' + str(phase(obs10002013001[1])) + ' +- ' + str(phaseErr(obs10002013001[1])))
print('\n')
print(Tnerr[eclipseIndices[0]]*24.0*3600.0)

# Eclipses immediately preceding and following observation 10002013003
eclipses = [56143.095015866362, 56146.98681448854]
eclipseIndices = []
for i in range(len(Tn)):
	if Tn[i]==eclipses[0]:
		eclipseIndices.append(i)
	if Tn[i]==eclipses[1]:
		eclipseIndices.append(i)
#print(eclipseIndices)
print('Obs 10002013003:')
print('Start: ' + str(phase(obs10002013003[0])) + ' +- ' + str(phaseErr(obs10002013003[0])))
print('End: ' + str(phase(obs10002013003[1])) + ' +- ' + str(phaseErr(obs10002013003[1])))
print('\n')
print(Tnerr[eclipseIndices[0]]*24.0*3600.0)

