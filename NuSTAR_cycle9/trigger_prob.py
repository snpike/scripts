import numpy as np 
import matplotlib.pyplot as plt
import astropy.coordinates as coord 
from astropy.time import Time
from scipy.stats import poisson
import seaborn as sns

sns.set_context('paper', font_scale=2.0)
sns.set_style("whitegrid")
sns.set_palette("colorblind")

sun_angle = []
trig_num = []
years = 11
exp_time = 40
n_max = 3
known = 0
unknown = 0

with open("trigger_cat.txt") as file:
    for line in file:
        splitline = line.split()
        # if (splitline[0] not in trig_num):
        #     trig_num.append(splitline[0])
        #     if (splitline[3] =='Unknown'):
        #         unknown += 1
        #     else:
        #         known += 1
        if (splitline[1][:2] != "11") and (splitline[1][:2] != "23") and (splitline[0] not in trig_num) and (splitline[3] =='Unknown') and ('WARNING' not in line):
            trig_num.append(splitline[0])
            pos = coord.SkyCoord(float(splitline[4]), float(splitline[5]), frame='icrs', unit='deg').transform_to(coord.GCRS)
            time_str = '20' + splitline[1][:2] + '-' + splitline[1][3:5] + '-' + splitline[1][-2:] + 'T' + splitline[2]
            time = Time(time_str, format='isot', scale='utc')
            sun_pos = coord.get_sun(time)
            sun_angle.append(pos.separation(sun_pos).value)

print('Total triggers: ' + str(len(trig_num)))
print('Known triggers: ' + str(known))
print('Unkown triggers: ' + str(unknown))

sun_angle = np.array(sun_angle)

plt.figure(figsize=(9,6))
plt.hist(sun_angle, bins=25)
plt.axvspan(10, 47, facecolor='red', alpha=0.3)
# plt.axvline(10, color='red', ls = 'dashed', lw=2.0)
# plt.axvline(47, color='red', ls = 'dashed', lw=2.0)
plt.xlim((0,180))
plt.ylabel('Number of MAXI triggers')
plt.xlabel('Angular separation from the Sun (degrees)')
plt.tight_layout()
# plt.show()
plt.savefig('/Users/sean/scripts/NuSTAR_cycle9/sun_angle_dist.pdf')


too_close_mask = sun_angle <= 10
trigger_mask = (sun_angle > 10) * (sun_angle <= 47)
too_far_mask = (sun_angle > 47)
print('Number of full years considered: ' + str(years))
print('Number of objects too close to the Sun: ' + str(np.sum(too_close_mask)))
print('Number of objects which would trigger: ' + str(np.sum(trigger_mask)))
print('Number of objects too far to the Sun: ' + str(np.sum(too_far_mask)))
print('Number of Unknown objects: ' + str(len(sun_angle)))

trig_rate = np.sum(trigger_mask)/years
print('Number of triggering objects per year: ' +str(trig_rate))
trig_dist = poisson(trig_rate)
trig_prob = 1-poisson.cdf(0.9, trig_rate)
# print(trig_prob)

total_weighted_time = 0
for i in range(n_max):
    total_weighted_time += exp_time * i * trig_dist.pmf(i)
    print('Prob of at least ' + str(i) + ' trig: ' + str(1-poisson.cdf(i-0.1, trig_rate)))

print('Prob of at least ' + str(n_max) + ' trig: ' + str(1-poisson.cdf(n_max-0.1, trig_rate)))

total_weighted_time += n_max * exp_time * (1-poisson.cdf(n_max-0.1, trig_rate))

print('Total weighted exp time: ' + str(total_weighted_time))




