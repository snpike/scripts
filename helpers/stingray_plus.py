import numpy as np
import seaborn as sns
import scipy

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.modeling import functional_models, fitting
from tqdm import tqdm

import matplotlib.pyplot as plt

import stingray.events as ev
import stingray.lightcurve as lc
from stingray import io
import stingray.powerspectrum as powspec 
import stingray.crossspectrum as crossspec
import stingray.gti as sting_gti
import stingray.pulse.pulsar as plsr
from stingray import stats
from stingray.pulse.search import z_n_search

class Lightcurve_ext(lc.Lightcurve):
    # Extend the stingray Lightcurve class to also hold onto an extraction region 

    def __init__(self, time, counts, err=None, input_counts=True, gti=None, err_dist='poisson', mjdref=0, dt=None, skip_checks=True, low_memory=False, centroid=None, radius=None):
        
        super().__init__(time=time, counts=counts, err=err, input_counts=input_counts, gti=gti, err_dist=err_dist, mjdref=mjdref, dt=dt, skip_checks=skip_checks, low_memory=low_memory)
        self.centroid=centroid
        self.radius=radius
    
    @staticmethod
    def make_lightcurve(toa, dt, tseg=None, tstart=None, gti=None, mjdref=0, use_hist=False, centroid=None, radius=None):
        temp_curve = lc.Lightcurve.make_lightcurve(toa, dt, tseg=tseg, tstart=tstart, gti=gti, mjdref=mjdref, use_hist=use_hist)
        return Lightcurve_ext(temp_curve.time, temp_curve.counts, err=temp_curve.counts_err, input_counts=True, gti=gti, err_dist='poisson', mjdref=mjdref, dt=dt, centroid=centroid, radius=radius)

    def calc_area(self):
        if self.radius:
            return np.pi*np.square(self.radius)
        else:
            print('The radius has not been defined')
            return None

    def rebin(self, dt_new=None, f=None, method='sum'):
        # for x in self.time:
        #     tmp_mask = (np.abs(self.time - x) <= dt_new/2) 
        #     new_counts.append(np.mean(bkg_lc.counts[tmp_mask]))
        #     bkg_err.append(np.sqrt(np.sum(np.square(bkg_lc.counts_err[tmp_mask])))/np.sum(tmp_mask))

        rebinned_lc = super().rebin(dt_new=dt_new, f=f, method=method)
        return Lightcurve_ext(rebinned_lc.time, rebinned_lc.counts, gti=self.gti, err_dist='poisson', mjdref=self.mjdref, dt=dt_new, centroid=self.centroid, radius=self.radius)

class EventList_ext(ev.EventList):
    # Extend the stingray EventList class to also hold onto the PRIOR column, 
    # which shows livetime since last event, and X and Y position columns.

    def __init__(self, time=None, energy=None, ncounts=None, mjdref=0, dt=0, notes="", gti=None, pi=None, prior=None, x=None, y=None, xy_weights = None, centroid=None, radius=None):
        
        super().__init__(time=time, energy=energy, ncounts=ncounts, mjdref=mjdref, dt=dt, notes=notes, gti=gti, pi=pi)
        self.centroid = centroid
        self.radius = radius
        self.prior = prior
        self.x = x
        self.y = y
        self.xy_weights = np.ones(np.shape(self.time))

    def join(self, other):
        ev_temp = super().join(other)

        sorted_arg = np.argsort(np.concatenate([self.time, other.time]))

        temp_x = np.concatenate([self.x, other.x])
        temp_x = temp_x[sorted_arg]

        temp_y = np.concatenate([self.y, other.y])
        temp_y = temp_y[sorted_arg]

        temp_xy_weights = np.concatenate([self.xy_weights, other.xy_weights])
        temp_xy_weights = temp_xy_weights[sorted_arg]

        temp_prior = np.concatenate([self.prior, other.prior])
        temp_prior = temp_prior[sorted_arg]

        # temp_pi = np.concatenate([self.pi, other.pi])
        # temp_pi = temp_pi[sorted_arg]

        ev_new = EventList_ext(time=ev_temp.time, energy=ev_temp.energy, ncounts=ev_temp.ncounts, mjdref=ev_temp.mjdref, dt=ev_temp.dt, gti=ev_temp.gti, pi=ev_temp.pi, prior=temp_prior, \
            x=temp_x, y=temp_y, xy_weights = temp_xy_weights)

        return ev_new
        
    def get_times(self, PI_min = 35, PI_max = 1909):
        PI_mask = (self.pi >= PI_min) * (self.pi <= PI_max)        
        return self.time[PI_mask]
    
    def livetime_curve(self, start_time, end_time):
        if (start_time < np.min(self.time)):
            print('Invalid start time')
            return None
        elif (end_time > np.max(self.time)):
            print('Invalid end time')
            return None
        else:
            time_mask = (self.time >= start_time) * (self.time <= end_time)
            temp_times = self.time[time_mask]
            temp_prior = self.prior[time_mask]
            live_fraction = temp_prior[1:]/(temp_times[1:]-temp_times[:-1])
            return temp_times, live_fraction
        
    def make_image(self, plot=False):
        # Produce a 2D histogram of events

        H, xedges, yedges = np.histogram2d(self.x, self.y, bins=500)
        H = H.T
        if plot:
            X, Y = np.meshgrid(xedges, yedges)
            plt.figure(figsize=(9,6))
            plt.pcolormesh(X, Y, H)
            
        xcenters, ycenters = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)
        return H, xcenters, ycenters
    
    def set_xy_weights(self, centroid=[520, 460]):
        # We can weight the events based on their distance from the source.
        # It's fine to just use region filtering instead, but in the case of really high count rates this might be useful.

        H, xcenters, ycenters = self.make_image(plot=False)
        g_init = functional_models.Gaussian2D(amplitude=np.max(H), x_mean=centroid[0], y_mean=centroid[1],\
                                                  x_stddev=10.0, y_stddev=10.0)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, xcenters, ycenters, H)
        self.xy_weights = g(self.x, self.y)/g.amplitude
        return self.xy_weights
        
        
    def fold_events_ltcorr(self, *frequency_derivatives, time_intervals= None, pi_min=35, pi_max=1909, \
                           region_filter=False, centroid = None, radius=None, ref_time=None, nbin = 64, weights  = 1, gtis = None, expocorr=False, weight_pos=False):
        # Do Epoch folding to look for pulsations while also correction for livetime variations. This is important for high count rates and high pulse fractions.
        # Includes region and energy filtering, position weighting, and custom weighting.

        if centroid == None:
            centroid = self.centroid
        if radius == None:
            radius = self.radius
        if time_intervals == None:
            time_intervals = [[self.time[0], self.time[-1]]]
        if ref_time == None:
            ref_time = self.time[0]
        if gtis == None:
            gtis = self.gti
        
        time_mask = np.zeros(np.shape(self.time))
        if np.shape(time_intervals)[-1] != 2:
            print('The array of time intervals has the wrong shape')
            return None
        for interval in time_intervals:
                start_time, end_time = interval
                if (start_time < np.min(self.time)):
                    print('Invalid start time')
                    return None
                elif (end_time > np.max(self.time)):
                    print('Invalid end time')
                    return None
                time_mask = time_mask + ((self.time >= start_time) * (self.time <= end_time))
        time_mask = time_mask.astype(bool)
        pi_mask = ((self.pi > pi_min) * (self.pi < pi_max)).astype(bool)
    
        reg_mask = np.ones(np.shape(self.time)).astype(bool)
        p_weights = 1.0
        if region_filter:
            if weight_pos:
                print('Region filtering overrides position weighting.')
            reg_mask = (np.sqrt(np.square(self.x-centroid[0]) + np.square(self.y-centroid[1])) < radius).astype(bool)
        elif weight_pos:
#             print('Make sure you have called set_xy_weights')
            p_weights = self.xy_weights[time_mask*pi_mask*reg_mask]
#             print([g.x_mean, g.y_mean, g.amplitude, g.x_stddev, g.y_stddev])
       
        # The times to actually fold into a profile    
        fold_times = self.time[time_mask*pi_mask*reg_mask]
        # The phase of each folded event
        fold_phases = plsr.pulse_phase(fold_times, *frequency_derivatives)
        
        # We should use PRIOR from every event though
        temp_times = self.time[time_mask]
        temp_prior = self.prior[time_mask]
        temp_phases = plsr.pulse_phase(temp_times, *frequency_derivatives)

#         print(frequency_derivatives)
        p_derivs = plsr.p_to_f(*frequency_derivatives)
        p_t = np.zeros(np.shape(temp_times))
        # Taylor expand P(t)
        for i in range(len(p_derivs)):
            p_t = p_t + (p_derivs[i] * np.power(temp_times-ref_time, i)/scipy.special.factorial(i, exact=True))
        
        phase_bins, profile, profile_err = plsr.fold_events(fold_times, *frequency_derivatives, ref_time=ref_time, \
                                                           nbin=nbin, weights=weights * p_weights, gtis=gtis, expocorr=expocorr)
        livetime_profile = np.zeros(np.shape(profile))
        
        # During which phase bin did PRIOR start counting before each event?
        start_bins = np.floor(nbin * (temp_phases - (temp_prior/p_t)))
        # What phase bin is each event at?
        end_bins = np.floor(nbin * temp_phases)
        
        # Add 1 to every phase bin for which PRIOR was active.
        for i in range(len(start_bins)):
            start = start_bins[i]
            end = end_bins[i]
#             print(start)
            # If PRIOR started counting in a previous cycle, add 1 to each full cycle and to the partial cycles
            if start < 0:
                for j in range(int(np.floor(np.abs(start)/nbin))):
                    livetime_profile = livetime_profile + 1
                livetime_profile[int(start%nbin):] = livetime_profile[int(start%nbin):] + 1
                livetime_profile[:int(end)] = livetime_profile[:int(end)] + 1
                
            # Else, just add 1 to the portion of this cycle during which PRIOR was counting
            else:
                livetime_profile[int(start):int(end)] = livetime_profile[int(start):int(end)] + 1
        
        # livetime corresponding to the phase bin for each photon
        fold_lts = np.array([livetime_profile[int(b)] for b in np.floor(nbin * fold_phases)])
        
        z_stat = plsr.z_n(fold_phases, n=2, norm=weights * p_weights * np.max(livetime_profile)/fold_lts)
#         livetime_profile = livetime_profile/np.max(livetime_profile)
        
        return phase_bins, profile, profile_err, livetime_profile, z_stat
    
    def fold_events(self, *frequency_derivatives, time_intervals= None, pi_min=35, pi_max=1909, 
                    region_filter=False, centroid = None, radius=None, ref_time=None, nbin = 64, weights = 1, gtis = None, expocorr=False, weight_pos=False, z_n=2):
        
        # Epoch folding without livetime correction.
        # Includes region and energy filtering, position weighting, and custom weighting.

        if region_filter or weight_pos:
            if centroid == None:
                centroid = self.centroid
            if radius == None:
                radius = self.radius
        if time_intervals == None:
            time_intervals = [[self.time[0], self.time[-1]]]
        if ref_time == None:
            ref_time = self.time[0]
        if gtis == None:
            gtis = self.gti
        
        time_mask = np.zeros(np.shape(self.time))
        if np.shape(time_intervals)[-1] != 2:
            print('The array of time intervals has the wrong shape')
            return None
        for interval in time_intervals:
                start_time, end_time = interval
                if (start_time < np.min(self.time)):
                    print('Invalid start time')
                    return None
                elif (end_time > np.max(self.time)):
                    print('Invalid end time')
                    return None
                time_mask = time_mask + ((self.time >= start_time) * (self.time <= end_time))
        time_mask = time_mask.astype(bool)
        pi_mask = ((self.pi > pi_min) * (self.pi < pi_max)).astype(bool)
    
        reg_mask = np.ones(np.shape(self.time)).astype(bool)
        p_weights = 1.0
        if region_filter:
            # if weight_pos:
                # print('Region filtering overrides position weighting.')
            reg_mask = (np.sqrt(np.square(self.x-centroid[0]) + np.square(self.y-centroid[1])) < radius).astype(bool)
        elif weight_pos:
#             print('Make sure you have called set_xy_weights')
            p_weights = self.xy_weights[time_mask*pi_mask*reg_mask]
#             print([g.x_mean, g.y_mean, g.amplitude, g.x_stddev, g.y_stddev])
       
        # The times to actually fold into a profile    
        fold_times = self.time[time_mask*pi_mask*reg_mask]
        # The phase of each folded event
        fold_phases = plsr.pulse_phase(fold_times, *frequency_derivatives)
        
        phase_bins, profile, profile_err = plsr.fold_events(fold_times, *frequency_derivatives, ref_time=ref_time, \
                                                           nbin=nbin, weights=weights * p_weights, gtis=gtis, expocorr=expocorr)
        
        
        z_stat = plsr.z_n(fold_phases, n=z_n, norm=weights * p_weights)
        
        return phase_bins, profile, profile_err, z_stat

    
    def split_by_gti(self):
        # Split up the Events object into multiple Events instances based on its gti.

        split_ev = []
        for g in self.gti:
            g_mask = (self.time <= g[1]) * (self.time >= g[0])
            split_ev.append(EventList_ext(time=self.time[g_mask], mjdref=self.mjdref, \
                                          dt=self.dt, notes=self.notes, gti=[g], pi=self.pi[g_mask], prior=self.prior[g_mask], \
                                          x=self.x[g_mask], y=self.y[g_mask], xy_weights = self.xy_weights[g_mask]))
        return split_ev

    
    def split_by_time(self, bintime=100, gti = None):
        # Split up the Events object into multiple Events instances based on the input time.

        split_ev = []
        if gti is None:
            gti = self.gti
        else:
            gti = sting_gti.cross_two_gtis(self.gti, gti)
        for g in gti:
            g_len = g[1]-g[0]
            for i in range(int(np.floor(g_len/bintime))):
                g_mask = (self.time <= g[0]+((i+1)*bintime)) * (self.time >= g[0]+(i*bintime))
                split_ev.append(EventList_ext(time=self.time[g_mask], mjdref=self.mjdref, \
                                  dt=self.dt, notes=self.notes, gti=[[g[0]+(i*bintime), g[0]+((i+1)*bintime)]], pi=self.pi[g_mask], prior=self.prior[g_mask], \
                                  x=self.x[g_mask], y=self.y[g_mask], xy_weights = self.xy_weights[g_mask]))
        return split_ev
    
    def to_lc(self, dt, pi_low = 35, pi_high = 1909, centroid=None, radius = None, tstart=None, tseg=None, gti=None, buff = False, buffersize = 100.0):
        # Bin this Events instance into a Lightcurve object.
        # Unlike the built-in version, this includes region and energy filtering, and you can introduce a new gti.

        if gti is None:
            gti = self.gti

        else:
            gti = sting_gti.cross_two_gtis(self.gti, gti)

        if buff:
            buffered_gti = []
            for x,y in gti:
                if np.abs(y-x) > 2*buffersize:
                    buffered_gti.append([x+buffersize, y-buffersize])
            gti = buffered_gti

        if tstart is None and gti is not None:
            tstart = gti[0][0]
            tseg = gti[-1][1] - tstart
        
        reg_mask = np.ones(np.shape(self.time)).astype(bool)
        if (centroid is not None) and (radius is not None):
            reg_mask = (np.sqrt(np.square(self.x-centroid[0]) + np.square(self.y-centroid[1])) < radius).astype(bool)

        
        pi_mask = (self.pi > pi_low) * (self.pi <= pi_high)
        tot_mask = reg_mask * pi_mask
        return Lightcurve_ext.make_lightcurve(self.time[tot_mask], dt, tstart=tstart, gti=gti, tseg=tseg, mjdref=self.mjdref, centroid=centroid, radius=radius)


        
        
        
def nuproducts_to_stingray_lc(nuprod_lc_file, rebin = False, buff = False, rebinsize=1.0, buffersize = 100.0, rebin_method='sum'):
    # Take a light curve produced by nuproducts and turn it into a Stingray Lightcurve instance.

    nuproducts_lc = fits.open(nuprod_lc_file)
    time = nuproducts_lc[1].data['TIME'] + nuproducts_lc[1].header['TSTART']
    dt = nuproducts_lc[1].header['TIMEDEL']
    countrate = nuproducts_lc[1].data['RATE']
    error = nuproducts_lc[1].data['ERROR']
    mjdref = nuproducts_lc[1].header['MJDREFI'] + nuproducts_lc[1].header['MJDREFF']
    gti = [[x,y] for x,y in nuproducts_lc[2].data]
    centroid = [nuproducts_lc[3].data['X'][0], nuproducts_lc[3].data['Y'][0]]
    radius = nuproducts_lc[3].data['R'][0]

    # Buffer around the orbital gaps if specified.
    if buff:
        buffered_gti = []
        for x,y in gti:
            if np.abs(y-x) > 2*buffersize:
                buffered_gti.append([x+buffersize, y-buffersize])
        gti = buffered_gti
    if rebin:
        return Lightcurve_ext(time, countrate, err=error, gti=gti, mjdref=mjdref, dt = dt, \
            input_counts=False, skip_checks=True, centroid=centroid, radius=radius).rebin(dt_new=rebinsize, method=rebin_method)
    else:
        return Lightcurve_ext(time, countrate, err=error, gti=gti, mjdref=mjdref, dt = dt, input_counts=False, skip_checks=True, centroid=centroid, radius=radius)
    
    
    
def extract_events(file_A, file_B, buff=0):
    # Extracts events for FPMA and FPMB from .evt files (assumed clean).
    # The gtis are crossed to make sure they are the same.
    # Can also use stingray.io.load_events_and_gtis, but this takes in X, Y, and PRIOR.
    
    ev_files = [fits.open(file_A), fits.open(file_B)] 
    
    ev_data_A = ev_files[0][1].data
    ev_data_B = ev_files[1][1].data
        
    ev_gti_A = [[x,y] for x,y in ev_files[0][2].data]
    ev_gti_B = [[x,y] for x,y in ev_files[1][2].data]
    common_gti = sting_gti.cross_two_gtis(ev_gti_A, ev_gti_B)
    buffered_gti = []
    if buff > 0:
        for x,y in common_gti:
            if np.abs(y-x) > 2*buff:
                buffered_gti.append([x+buff, y-buff])
    else:
        buffered_gti = common_gti
        

    events = [EventList_ext(time=ev_data_A['TIME'], gti=buffered_gti, pi = ev_data_A['PI'], \
                       mjdref=ev_files[0][0].header['MJDREFI'] + ev_files[0][0].header['MJDREFF'], \
                           prior=ev_data_A['PRIOR'], x=ev_data_A['X'], y=ev_data_A['Y']), \
              EventList_ext(time=ev_data_B['TIME'], gti=buffered_gti, pi = ev_data_B['PI'], \
                       mjdref=ev_files[1][0].header['MJDREFI'] + ev_files[1][0].header['MJDREFF'], \
                           prior=ev_data_B['PRIOR'], x=ev_data_B['X'], y=ev_data_B['Y'])]
    ev_files[0].close()
    ev_files[1].close()
    return events

def split_ev_by_gti(events):
    # Split up an Events instance into multiple Events instances based on the gtis.
    # For use with built-in EventList instances, not EventList_ext, which has its own class-specific version.
    split_ev = []
    gti = events.gti
    for g in gti:
        split_ev.append(ev.EventList(time=events.time, gti=[g], pi = events.pi, \
                       mjdref=events.mjdref))
    return split_ev

def split_ev_by_time(events, bintime=100):
    # Split up an Events instance into multiple Events instances based on the input time.
    # For use with built-in EventList instances, not EventList_ext, which has its own class-specific version.
    split_ev = []
    gti = events.gti
    for g in gti:
        g_len = g[1]-g[0]
        for i in range(int(np.floor(g_len/bintime))):
            g_mask = (events.time <= g[0]+((i+1)*bintime)) * (events.time >= g[0]+(i*bintime))
            split_ev.append(ev.EventList(time=events.time[g_mask], gti=[[g[0]+(i*bintime), g[0]+((i+1)*bintime)]], \
                                         pi = events.pi[g_mask], mjdref=events.mjdref))
    return split_ev

def new_gti(events, gti):
    # Makes a new EventList instance with the specified gti. 
    # TODO: Make this work with EventList_ext
    return ev.EventList(time=events.time, gti=gti, pi = events.pi, \
                       mjdref=events.mjdref)
    

def read_QDP(file):
#     Takes in a file and returns an array with the QDP data
    data = [[]]
    i = 0
    j = 0
    for line in file:
        if i >2:
            temp = line.split()
            if temp[0] != 'NO':
                if temp[-1]=='NO':
                    temp = temp[:-1]
                data[-1].append(temp)
            else:
                data.append([])
        i += 1
    for i in range(len(data)):
        data[i] = np.array(data[i]).astype(float).T
    return data
    
def calc_pvals(power,nspec):
    # Calculates the pvalue for a given power and number of spectra, nspec.
    exp_sigma = np.sqrt(2) / np.sqrt(nspec)
    gauss = scipy.stats.norm(0, exp_sigma)
    pval = gauss.sf(power)
    return pval

def sum_lc(lc_1, lc_2):
    # Add up two Lightcurve instances. Must be simultaneous and bins must line up.

    common_gti = sting_gti.cross_two_gtis(lc_1.gti, lc_2.gti)
    lc_1.gti = common_gti
    lc_2.gti = common_gti
    lc_1.apply_gtis()
    lc_2.apply_gtis()
    if np.sum(lc_1.time) != np.sum(lc_2.time):
        print('Lightcurves don\'t line up. Exiting')
        return None

    area_ratio = lc_1.calc_area()/lc_2.calc_area()
    summed_err = np.sqrt(np.square(lc_1.counts_err) + np.square(lc_2.counts_err*area_ratio))

    # Handle edge case where counts=0. TODO: more accurate error estimation.
    summed_err[summed_err==0.0] = 1.0

    # New Lightcurve
    summed_lc = Lightcurve_ext(lc_1.time, lc_1.counts + (lc_2.counts*area_ratio), err=summed_err, gti=common_gti, mjdref=lc_1.mjdref, dt = lc_1.dt, input_counts=True, skip_checks=True, centroid=None, radius=lc_1.radius)
    return summed_lc

def bkg_subtract(src_lc, bkg_lc, bkg_bin=5):
    # Subtract the moving average of a background lightcurve from a source light curve.
    # bkg_bin is the moving average window in seconds.

    bkg_counts = []
    bkg_err = []
    outlier_mask = np.abs(bkg_lc.counts - np.mean(bkg_lc.counts)) < 5*np.std(bkg_lc.counts)
    for x in src_lc.time:
        tmp_mask = (np.abs(bkg_lc.time - x) <= bkg_bin/2) * outlier_mask
        bkg_counts.append(np.mean(bkg_lc.counts[tmp_mask]))
        bkg_err.append(np.sqrt(np.sum(np.square(bkg_lc.counts_err[tmp_mask])))/np.sum(tmp_mask))

    area_ratio = src_lc.calc_area()/bkg_lc.calc_area()
    bkg_counts = np.array(bkg_counts) * area_ratio
    bkg_err = np.array(bkg_err) * area_ratio

    return Lightcurve_ext(src_lc.time, src_lc.counts - bkg_counts, dt = src_lc.dt, err=np.sqrt(np.square(src_lc.counts_err) + np.square(bkg_err)), gti=src_lc.gti, mjdref=src_lc.mjdref, input_counts=True, skip_checks=True, centroid=src_lc.centroid, radius=src_lc.radius)
    
def efold_search_ltcorr(events, f_min, f_max, f_steps, time_intervals=None, nbin = 32, pi_min = 35, pi_max = 260, region_filter=True, weight_pos=False):
    # Scan over frequency and do livetime corrected epoch folding.

    f_arr = np.linspace(f_min, f_max, num = f_steps)
    z_prob = []
    for f in tqdm(f_arr):
        phase_bins, profile, profile_err, livetime_profile, z_stat = \
            events.fold_events_ltcorr(f, time_intervals = time_intervals, \
                                    nbin = nbin, ref_time = events.time[0], region_filter=region_filter, pi_min=pi_min, pi_max=pi_max, weight_pos=weight_pos)
        #     phase_bins_B, profile_B, profile_err_B, livetime_profile_B, z_stat_B = \
        #         events[1].fold_events_ltcorr(f, time_intervals = [[(np.min(events[0].time) + 24150), (np.min(events[0].time) + 24180)], [(np.min(events[0].time) + 37890), (np.min(events[0].time) + 37900)]], \
        #                                 nbin = nbin, ref_time = events[0].time[0], region_filter=False, pi_min=35, pi_max=260)

        #     summed_profile = np.mean(livetime_profile_A + livetime_profile_B) * (profile_A + profile_B)/(livetime_profile_A + livetime_profile_B)
        #     summed_err = np.mean(livetime_profile_A + livetime_profile_B) * np.sqrt(np.square(profile_err_A) + np.square(profile_err_B))/(livetime_profile_A + livetime_profile_B)
        #     summed_stat = stat.z2_n_probability(z_stat_A, err=summed_err)
        #     summed_real = 1.0 - plsr.fold_profile_probability(summed_stat, nbin)

        eps = stats.z2_n_probability(float(z_stat), n=2, ntrial=len(f_arr))
        z_prob.append(1.0-eps)

    z_prob = np.array(z_prob)
    return f_arr, z_prob

def efold_search(events, f_min, f_max, f_steps, ref_time=None, time_intervals=None, nbin = 32, pi_min = 35, pi_max = 260, region_filter=False, weight_pos=False, z_n=2):
    # Scan over frequency and do livetime corrected epoch folding.

    if ref_time is None:
        ref_time=events.time[0]

    f_arr = np.linspace(f_min, f_max, num = f_steps)
    z_prob = []
    z_stat_arr =[]
    for f in tqdm(f_arr):
        phase_bins, profile, profile_err, z_stat = \
            events.fold_events(f, time_intervals = time_intervals, \
                                    nbin = nbin, ref_time = ref_time, region_filter=region_filter, pi_min=pi_min, pi_max=pi_max, weight_pos=weight_pos, z_n=z_n)
        z_stat_arr.append(z_stat)
        eps = stats.z2_n_probability(float(z_stat), n=z_n, ntrial=len(f_arr))
        z_prob.append(1.0-eps)

    z_stat_arr = np.array(z_stat_arr)
    z_prob = np.array(z_prob)
    return f_arr, z_prob, z_stat_arr


def efold_search_AandB(events_A, events_B, f_min, f_max, f_steps, fdots=None, time_intervals=None, nbin = 32, pi_min = 35, pi_max = 260, return_peak = False, z_n=2):
    # Scan over frequency and do epoch folding.

    A_mask = np.sqrt(np.square(events_A.x - events_A.centroid[0]) + np.square(events_A.y - events_A.centroid[1])) <= events_A.radius
    B_mask = np.sqrt(np.square(events_B.x - events_B.centroid[0]) + np.square(events_B.y - events_B.centroid[1])) <= events_B.radius

    temp_time = np.concatenate([events_A.time[A_mask], events_B.time[B_mask]])
    sorted_arg = np.argsort(temp_time)
    temp_time = temp_time[sorted_arg]
    temp_pi = np.concatenate([events_A.pi[A_mask], events_B.pi[B_mask]])[sorted_arg]

    joined_ev = EventList_ext(time=temp_time, gti=sting_gti.cross_two_gtis(events_A.gti, events_B.gti) , pi=temp_pi)
    ref_time = joined_ev.time[0]
    f_arr = np.linspace(f_min, f_max, num = f_steps)

    # if fdots:
    #     fgrid, fdgrid, z_stats = z_n_search(joined_ev.time, f_arr, nharm=z_n, nbin=nbin, gti=joined_ev.gti, fdots=fdots, segment_size=1e6)
    # else:
    #     fgrid, z_stats = z_n_search(joined_ev.time, f_arr, nharm=z_n, nbin=nbin, gti=joined_ev.gti, fdots=fdots, segment_size=1e6)
    
    pi_mask = ((joined_ev.pi > pi_min) * (joined_ev.pi < pi_max)).astype(bool)
    # The times to actually fold into a profile    
    fold_times = joined_ev.time[pi_mask] - ref_time
    z_stats=[]
    if fdots is not None:
        f_grid, fd_grid = np.meshgrid(f_arr, fdots)
        z_stats = np.zeros(f_grid.shape)
        for x in tqdm(range(f_steps)):
            for y in range(len(fdots)):
                # The phase of each folded event
                fold_phases = plsr.pulse_phase(fold_times, *[f_grid[y,x], fd_grid[y,x]])
                z_stats[y,x] = plsr.z_n(fold_phases, n=z_n)

        z_prob = stats.z2_n_logprobability(z_stats, ntrial=(f_steps*len(fdots)), n=z_n)

        if return_peak:
            max_yx =  np.unravel_index(np.argmax(z_stats, axis=None), z_stats.shape)
            phase_bins, profile, profile_err, _ = \
                joined_ev.fold_events(*[f_grid[max_yx], fd_grid[max_yx]], time_intervals = time_intervals, \
                                        nbin = nbin, ref_time = ref_time, region_filter=False, pi_min=pi_min, pi_max=pi_max, weight_pos=False, z_n=z_n)

            return f_grid, fd_grid, z_prob, z_stats, phase_bins, profile, profile_err
        
        else:
            return f_grid, fd_grid, z_prob, z_stats

    else: 
        for f in tqdm(f_arr):
            # The phase of each folded event
            fold_phases = plsr.pulse_phase(fold_times, f)
            
            z_stat = plsr.z_n(fold_phases, n=z_n)

            # _, _, _, z_stat = \
            #     joined_ev.fold_events(f, time_intervals = time_intervals, \
            #                             nbin = nbin, ref_time = joined_ev.time[0], region_filter=False, pi_min=pi_min, pi_max=pi_max, weight_pos=False, z_n=z_n)

            z_stats.append(z_stat)
        
        z_stats = np.array(z_stats)
        z_prob = stats.z2_n_logprobability(z_stats, ntrial=len(f_arr), n=z_n)


        if return_peak:
            phase_bins, profile, profile_err, _ = \
                joined_ev.fold_events(f_arr[np.argmax(z_stats)], time_intervals = time_intervals, \
                                        nbin = nbin, ref_time = ref_time, region_filter=False, pi_min=pi_min, pi_max=pi_max, weight_pos=False, z_n=z_n)

            return f_arr, z_prob, z_stats, phase_bins, profile, profile_err
        
        else:
            return f_arr, z_prob, z_stats


def PI_to_eV(PI):
    return (PI*40) + 1600

def eV_to_PI(eV):
    return (eV-1600)/40

# def minimize_remainder(arr, min_div, max_div):
#     divisors = np.linspace(min_div, max_div, num=100)
#     remainders = []
#     for div in divisors:
#         remainders.append(np.sum(np.mod(arr, div)))
        
#     return divisors[np.argmin(remainders)]

def power_law(f, B, gamma):
    return B*np.power(f,gamma)

def Lorentzian(f, peakf, Q, A):
    # gamma = HWHM
    # peakf = centroid frequency
    gamma = peakf/(2.0 * Q)
    return (A * np.square(gamma)/(np.pi*gamma*(np.square(f-peakf) + np.square(gamma))))

def zero_center_Lorentzian(f, gamma, A):
    # gamma = HWHM
    # peakf = centroid frequency
    return (A * np.square(gamma)/(np.pi*gamma*(np.square(f) + np.square(gamma))))

# def Lorentzian_C(f, peakf, Q, A, C):
#     return Lorentzian(f, peakf, Q, A) + C

def Lorentzian_power(f, peakf, Q, A, B):
    return Lorentzian(f, peakf, Q, A) + power_law(f, B)

# def N_Lorentzian(f, *args):
#     N = int(len(args)/3)
#     peak_nu = args[:N]
#     Qs = args[N:N+N]
#     As = args[N+N:N+N+N]
#     model = np.zeros(np.shape(f))
#     for i in range(N):
#         model = model + Lorentzian(f, peak_nu[i], Qs[i], As[i])
        
#     return model

# def N_Lorentzian_C(f, *args):
#     N = int((len(args)-1)/3)
#     peak_nu = args[:N]
#     Qs = args[N:N+N]
#     As = args[N+N:N+N+N]
#     C = args[-1]
#     model = C * np.ones(np.shape(f))
#     for i in range(N):
#         model = model + Lorentzian(f, peak_nu[i], Qs[i], As[i])
        
#     return model

# def N_Lorentzian_power(f, *args):
#     N = int((len(args)-2)/3)
#     peak_nu = args[:N]
#     Qs = args[N:N+N]
#     As = args[N+N:N+N+N]
#     B, alpha = args[-2:]
#     model = power_law(f, B, alpha)
#     for i in range(N):
#         gamma = peak_nu[i]/(2.0 * Qs[i])
#         model = model + (As[i] * np.square(gamma)/(np.pi*gamma*(np.square(f-peak_nu[i]) + np.square(gamma))))
        
#     return model

def QPO_scan(cross_spec, f_min=1e-4, f_max=2000., f_bin=1000, n_lorentz = 1):
    f_mask = (cross_spec.freq > f_min) * (cross_spec.freq < f_max)
    freq_steps = np.logspace(np.log10(cross_spec.freq[f_mask][0]), np.log10(cross_spec.freq[f_mask][-1]), f_bin + 2)
    xdata = cross_spec.freq[f_mask]
    ydata = cross_spec.power[f_mask]
    sigma = cross_spec.power_err[f_mask]
    
    pl_popt, pl_pcov = scipy.optimize.curve_fit(power_law, xdata, ydata, sigma = sigma, p0 = [1e-5], \
                                                bounds=np.array([(0.0), (np.inf)]))
    print(pl_popt)
    chisq0 = np.sum(((ydata - power_law(xdata, *pl_popt)) / sigma) ** 2)
    chisq = []
    for i in tqdm(range(len(freq_steps[1:-1]))):
        f = freq_steps[i+1]
        popt, pcov = scipy.optimize.curve_fit(Lorentzian_power, xdata, ydata, sigma = sigma, p0 = [f, 2.0, 0.1, pl_popt[0]], \
                                              bounds=np.array([(f - (f-freq_steps[i])/2., f + (freq_steps[i+2] - f)/2.0), (1.0,np.inf), (0.0,np.inf), (0.0, np.inf)]).T)
        chisq.append(np.sum(((ydata - Lorentzian_power(xdata, *popt)) / sigma) ** 2))
    dof = len(xdata)-len(popt)
    return freq_steps[1:-1], chisq0, np.array(chisq), dof



# def fit_peaks(xdata, ydata, sigma, nu_peak):
    
#     popt_arr = []
#     pcov_arr = []

#     for i, p in enumerate(nu_peak):
#         f_bound = None
#         if len(nu_peak)==1:
#             f_bound = (np.min(xdata), np.max(xdata))
#         else:
#             if i == 0:
#                 f_bound = (np.min(xdata), p + (nu_peak[i+1] - p)/2)
#             elif i==len(nu_peak)-1:
#                 f_bound = (p + (nu_peak[i-1] - p)/2, np.max(xdata))
#             else:
#                 f_bound = (p + (nu_peak[i-1] - p)/2, p + (nu_peak[i+1] - p)/2)

#         par_bounds = np.array([f_bound, (1.0,np.inf), (0, np.inf), (0, np.inf), (-np.inf, 0.0)]).T
#         p0 = [p, 5.0, 0.1, 0.02, -0.5]
#         popt, pcov = scipy.optimize.curve_fit(Lorentzian_power, xdata, ydata, sigma = sigma, p0 = p0, bounds = par_bounds)
#         popt_arr.append(popt)
#         pcov_arr.append(pcov)

#     popt_arr = np.array(popt_arr)
#     pcov_arr = np.array(pcov_arr)

#     return popt_arr, pcov_arr

# def get_rms(events, centroids, radius, PI_min=35, PI_max=1210, nu_min=1e-4, nu_max=100., split_time=512, bin_time = 1/1024, plot = False):
    
#     curve_A = events[0].to_lc(dt = bin_time, pi_low=PI_min, pi_high=PI_max, centroid = centroids[0], radius = radius)
#     curve_B = events[1].to_lc(dt = bin_time, pi_low=PI_min, pi_high=PI_max, centroid = centroids[1], radius = radius)
    
#     average_CPDS = crossspec.AveragedCrossspectrum(curve_A, curve_B, segment_size=split_time, norm = 'frac')
    
#     if plot:
#         f_res = 0.1
#         plt.figure(figsize=(9,6))
#         averaged_cross_log = average_CPDS.rebin_log(f=f_res)
#         averaged_cross_log_err = average_CPDS.df*np.power(1.+f_res, range(len(averaged_cross_log.freq)))/2.
#         plt.errorbar(averaged_cross_log.freq, averaged_cross_log.power*averaged_cross_log.freq, xerr=averaged_cross_log_err, yerr=averaged_cross_log.power_err*averaged_cross_log.freq, fmt='none', lw=0.5)
#         plt.xscale('log')
#         plt.ylim((1e-4,5*np.max(averaged_cross_log.power.real)))
#         plt.yscale('log')
#         plt.xlabel('Frequency (Hz)')
#         plt.ylabel('Leahy Power')
#         plt.ylabel(r'$\mathrm{\nu P_{\nu}\ (rms/mean)^{2}}$')
#         plt.tight_layout()
#         plt.show()
# #         plt.savefig(plot_dir + 'averaged_cross_spectrum_' + str(int(split_time)) + 's.pdf')
#         plt.close()
    
#     rms_square = np.sum(average_CPDS.power[(average_CPDS.freq > nu_min) * (average_CPDS.freq < nu_max)])*average_CPDS.df
#     rms_square_err = np.sqrt(np.sum(np.square(average_CPDS.power_err[(average_CPDS.freq > nu_min) * (average_CPDS.freq < nu_max)])))*average_CPDS.df
    
#     if rms_square < 0.0:
#         uplim = True
#         rms = np.sqrt(rms_square + 2.6*rms_square_err)
#     else:
#         uplim = False
#         rms = np.sqrt(rms_square)
        
#     rms_err = np.sqrt(np.sum(np.square(average_CPDS.power_err[(average_CPDS.freq > nu_min) * (average_CPDS.freq < nu_max)])))*average_CPDS.df/(2*rms)

#     return rms, rms_err, uplim

def model_continuum_noise_zero_center(cpds, plot=True, plot_dir='/Users/sean/Desktop/', plot_name='CPDS_Pnu_continuum_zero_center.pdf', nu_max = 1e5, f_res = 0.5):
    nu_mask = cpds.freq < nu_max
    chisq0 = np.sum(((cpds.power[nu_mask]-np.mean(cpds.power[nu_mask]))/ cpds.power_err[nu_mask]) ** 2)

    popt, pcov = scipy.optimize.curve_fit(zero_center_Lorentzian, cpds.freq[nu_mask], cpds.power[nu_mask], sigma = cpds.power_err[nu_mask], \
                                          p0 = [np.min(cpds.freq[nu_mask]), 10.], bounds= [[0.0,0.0], [np.inf, np.inf]])
    chisq = np.sum(((cpds.power[nu_mask] - zero_center_Lorentzian(cpds.freq[nu_mask], *popt)) / cpds.power_err[nu_mask]) ** 2)

    if plot:
        cpds_log = cpds.rebin_log(f=f_res)
        temp_err = cpds.df*np.power(1.+f_res, range(len(cpds_log.freq)))/2.

        # plt.figure(figsize=(9,6))
        # plt.errorbar(cpds_log.freq, cpds_log.power*cpds_log.freq, xerr=temp_err, \
        #              yerr=cpds_log.power_err*cpds_log.freq, fmt='none', lw=0.5, color='black')
        # plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
        #          np.concatenate([cpds_log.power*cpds_log.freq, [(cpds_log.power*cpds_log.freq)[-1]]]), where='post', color='black', lw=0.5)
        
        # plt.plot(cpds.freq,zero_center_Lorentzian(cpds.freq, *popt)*cpds.freq, color='red', lw=1.0)
        # plt.xscale('log')
        # plt.ylim((1e-6,3*np.max(cpds_log.power.real*cpds_log.freq)))
        # plt.yscale('log')
        # plt.xlabel('Frequency (Hz)')
        # plt.ylabel(r'$\mathrm{\nu P_{\nu}\ (rms/mean)^{2}}$')
        # plt.tight_layout()
        # plt.savefig(plot_dir + plot_name)
        # plt.close()

        plt.figure(figsize=(9,6))
        plt.errorbar(cpds_log.freq, cpds_log.power, xerr=temp_err, \
                     yerr=cpds_log.power_err, fmt='none', lw=0.5, color='black')
        plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
                 np.concatenate([cpds_log.power, [cpds_log.power[-1]]]), where='post', color='black', lw=0.5)
        
        plt.plot(cpds.freq,zero_center_Lorentzian(cpds.freq, *popt), color='red', lw=1.0)
        plt.xscale('log')
        plt.ylim((1e-6,3*np.max(cpds_log.power.real)))
        plt.yscale('log')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel(r'$\mathrm{Power\ (rms/mean)^{2}/Hz}$')
        plt.tight_layout()
        plt.savefig(plot_dir + plot_name)
        plt.close()
    
    return popt, pcov, chisq0, chisq

def model_continuum_noise_pl(cpds, plot=True, plot_dir='/Users/sean/Desktop/', plot_name='CPDS_Pnu_continuum_pl.pdf', nu_max = 1e5, f_res = 0.5):
    nu_mask = cpds.freq < nu_max
    chisq0 = np.sum(((cpds.power[nu_mask]-np.mean(cpds.power[nu_mask]))/ cpds.power_err[nu_mask]) ** 2)

    popt, pcov = scipy.optimize.curve_fit(power_law, cpds.freq[nu_mask], cpds.power[nu_mask], sigma = cpds.power_err[nu_mask], \
                                          p0 = [1e-5, -2], bounds= [[0.0, -np.inf], [np.inf, 0.0]])
    chisq = np.sum(((cpds.power[nu_mask] - power_law(cpds.freq[nu_mask], *popt)) / cpds.power_err[nu_mask]) ** 2)

    if plot:
        cpds_log = cpds.rebin_log(f=f_res)
        temp_err = cpds.df*np.power(1.+f_res, range(len(cpds_log.freq)))/2.

        # plt.figure(figsize=(9,6))
        # plt.errorbar(cpds_log.freq, cpds_log.power*cpds_log.freq, xerr=temp_err, \
        #              yerr=cpds_log.power_err*cpds_log.freq, fmt='none', lw=0.5, color='black')
        # plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
        #          np.concatenate([cpds_log.power*cpds_log.freq, [(cpds_log.power*cpds_log.freq)[-1]]]), where='post', color='black', lw=0.5)
        
        # plt.plot(cpds.freq,zero_center_Lorentzian(cpds.freq, *popt)*cpds.freq, color='red', lw=1.0)
        # plt.xscale('log')
        # plt.ylim((1e-6,3*np.max(cpds_log.power.real*cpds_log.freq)))
        # plt.yscale('log')
        # plt.xlabel('Frequency (Hz)')
        # plt.ylabel(r'$\mathrm{\nu P_{\nu}\ (rms/mean)^{2}}$')
        # plt.tight_layout()
        # plt.savefig(plot_dir + plot_name)
        # plt.close()

        plt.figure(figsize=(9,6))
        plt.errorbar(cpds_log.freq, cpds_log.power, xerr=temp_err, \
                     yerr=cpds_log.power_err, fmt='none', lw=0.5, color='black')
        plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
                 np.concatenate([cpds_log.power, [cpds_log.power[-1]]]), where='post', color='black', lw=0.5)
        
        plt.plot(cpds.freq,power_law(cpds.freq, *popt), color='red', lw=1.0)
        plt.xscale('log')
        plt.ylim((1e-6,3*np.max(cpds_log.power.real)))
        plt.yscale('log')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel(r'$\mathrm{Power\ (rms/mean)^{2}/Hz}$')
        plt.tight_layout()
        plt.savefig(plot_dir + plot_name)
        plt.close()
    
    return popt, pcov, chisq0, chisq



# def model_continuum_noise(cpds, plot=True, plot_dir='/Users/sean/Desktop/', f_res = 0.5):
#     n_L = 0
#     L_args = []
#     L_bounds = []
#     chisq0 = np.sum(((cpds.power-np.mean(cpds.power))/ cpds.power_err) ** 2)
#     chisq = []
#     popt_list = []
#     pcov_list = []
#     del_chisq = -100000
#     while del_chisq < 0.0:
#         n_L += 1
#         if n_L ==1:
#             L_args = np.array([1.0, 0.1, 10.])
#             L_bounds = np.array([[[0.0, np.max(cpds.freq)], [0.0, 2.0], [0,np.inf]]])
#         else:
#             L_args = np.vstack((L_args, [1.0, 0.1, 10.]))
#             L_bounds = np.concatenate((L_bounds, [[(0.0, np.max(cpds.freq)), (0.0, 10.0), (0,np.inf)]]))
                
#         temp_p0 = L_args.T.flatten()
#         temp_bounds = [L_bounds.T.flatten()[:3*n_L], L_bounds.T.flatten()[3*n_L:]]
        
#         popt, pcov = scipy.optimize.curve_fit(N_Lorentzian, cpds.freq, cpds.power, sigma = cpds.power_err, \
#                                               p0 = temp_p0, bounds= temp_bounds)
#         temp_chisq = np.sum(((cpds.power - N_Lorentzian(cpds.freq, *popt)) / cpds.power_err) ** 2)
        
#         if n_L==1:
#             del_chisq = temp_chisq-chisq0
#         else:
#             del_chisq = temp_chisq-chisq[-1]
#         chisq.append(temp_chisq)
#         popt_list.append(popt)
#         pcov_list.append(pcov)
        
#         L_args = np.array([popt[:n_L], popt[n_L:n_L + n_L], popt[n_L+n_L:n_L+n_L+n_L]]).T

#         if plot:
#             cpds_log = cpds.rebin_log(f=f_res)
#             temp_err = cpds.df*np.power(1.+f_res, range(len(cpds_log.freq)))/2.

#             plt.figure(figsize=(9,6))
#             plt.errorbar(cpds_log.freq, cpds_log.power*cpds_log.freq, xerr=temp_err, \
#                          yerr=cpds_log.power_err*cpds_log.freq, fmt='none', lw=0.5, color='black')
#             plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
#                      np.concatenate([cpds_log.power*cpds_log.freq, [(cpds_log.power*cpds_log.freq)[-1]]]), where='post', color='black', lw=0.5)
            
#             for i in range(n_L):
#                 plt.plot(cpds.freq, Lorentzian(cpds.freq, *((L_args[i])))*cpds.freq, color='red', ls='dotted', lw=1.0)
#             plt.plot(cpds.freq,N_Lorentzian(cpds.freq, *popt)*cpds.freq, color='red', lw=1.0)
#             plt.xscale('log')
#             plt.ylim((1e-6,3*np.max(cpds_log.power.real*cpds_log.freq)))
#             plt.yscale('log')
#             plt.xlabel('Frequency (Hz)')
#             plt.ylabel(r'$\mathrm{\nu P_{\nu}\ (rms/mean)^{2}}$')
#             plt.tight_layout()
#             plt.savefig(plot_dir + 'CPDS_nuPnu_continuum_' + str(int(n_L)) + '_comps.pdf')
#             plt.close()
            
#             plt.figure(figsize=(9,6))
#             plt.errorbar(cpds_log.freq, cpds_log.power, xerr=temp_err, \
#                          yerr=cpds_log.power_err, fmt='none', lw=0.5, color='black')
#             plt.step(np.concatenate([cpds_log.freq-temp_err, [cpds_log.freq[-1]+temp_err[-1]]]), \
#                      np.concatenate([cpds_log.power, [cpds_log.power[-1]]]), where='post', color='black', lw=0.5)

#             for i in range(n_L):
#                 plt.plot(cpds.freq, Lorentzian(cpds.freq, *((L_args[i]))), color='red', ls='dotted', lw=1.0)
#             plt.plot(cpds.freq,N_Lorentzian(cpds.freq, *popt), color='red', lw=1.0)
#             plt.xscale('log')
#             plt.ylim((1e-6,3*np.max(cpds_log.power.real)))
#             plt.yscale('log')
#             plt.xlabel('Frequency (Hz)')
#             plt.ylabel(r'$\mathrm{Power\ (rms/mean)^{2}/Hz}$')
#             plt.tight_layout()
#             plt.savefig(plot_dir + 'CPDS_Pnu_continuum_' + str(int(n_L)) + '_comps.pdf')
#             plt.close()
    
#     return n_L, popt_list, pcov_list, chisq0, chisq



