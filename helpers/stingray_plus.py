import numpy as np
import seaborn as sns
import scipy

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.modeling import functional_models, fitting

import stingray.events as ev
import stingray.lightcurve as lc
from stingray import io
import stingray.powerspectrum as powspec 
import stingray.crossspectrum as crossspec
from hendrics.efsearch import dyn_folding_search, z_n_search, folding_search
import stingray.gti as sting_gti
import stingray.pulse.pulsar as plsr
# from stingray import stats


class EventList_ext(ev.EventList):
    # Extend the stingray EventList class to also hold onto the PRIOR column, which shows livetime since last event
    def __init__(self, time=None, energy=None, ncounts=None, mjdref=0, dt=0, notes="", gti=None, pi=None, prior=None, x=None, y=None, xy_weights = None):
        
        super().__init__(time=time, energy=energy, ncounts=ncounts, mjdref=mjdref, dt=dt, notes=notes, gti=gti, pi=pi)
        self.prior = prior
        self.x = x
        self.y = y
        self.xy_weights = np.ones(np.shape(self.time))
        
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
        H, xedges, yedges = np.histogram2d(self.x, self.y, bins=500)
        H = H.T
        if plot:
            X, Y = np.meshgrid(xedges, yedges)
            plt.figure(figsize=(9,6))
            plt.pcolormesh(X, Y, H)
            
        xcenters, ycenters = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)
        return H, xcenters, ycenters
    
    def set_xy_weights(self, centroid=[520, 460]):
        H, xcenters, ycenters = self.make_image(plot=False)
        g_init = functional_models.Gaussian2D(amplitude=np.max(H), x_mean=centroid[0], y_mean=centroid[1],\
                                                  x_stddev=10.0, y_stddev=10.0)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, xcenters, ycenters, H)
        self.xy_weights = g(self.x, self.y)/g.amplitude
        return self.xy_weights
        
        
    def fold_events_ltcorr(self, *frequency_derivatives, time_intervals= None, pi_min=35, pi_max=1909, \
                           region_filter=False, centroid = [520,460], radius=46.783962, ref_time=None, nbin = 64, weights  = 1, gtis = None, expocorr=False, weight_pos=False):
        
        if centroid == None:
            centroid = self.centroid
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
                    region_filter=False, centroid = [520,460], radius=46.783962, ref_time=None, nbin = 64, weights = 1, gtis = None, expocorr=False, weight_pos=False):
        
        if centroid == None:
            centroid = self.centroid
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
        
#         print(frequency_derivatives)
        p_derivs = plsr.p_to_f(*frequency_derivatives)
        p_t = np.zeros(np.shape(temp_times))
        # Taylor expand P(t)
        for i in range(len(p_derivs)):
            p_t = p_t + (p_derivs[i] * np.power(temp_times-ref_time, i)/scipy.special.factorial(i, exact=True))
        
        phase_bins, profile, profile_err = plsr.fold_events(fold_times, *frequency_derivatives, ref_time=ref_time, \
                                                           nbin=nbin, weights=weights * p_weights, gtis=gtis, expocorr=expocorr)

        
        z_stat = plsr.z_n(fold_phases, n=2, norm=weights * p_weights)
#         livetime_profile = livetime_profile/np.max(livetime_profile)
        
        return phase_bins, profile, profile_err, livetime_profile, z_stat


    
    def split_by_gti(self):
        split_ev = []
        for g in self.gti:
            g_mask = (self.time <= g[1]) * (self.time >= g[0])
            split_ev.append(EventList_ext(time=self.time[g_mask], mjdref=self.mjdref, \
                                          dt=self.dt, notes=self.notes, gti=[g], pi=self.pi[g_mask], prior=self.prior[g_mask], \
                                          x=self.x[g_mask], y=self.y[g_mask], xy_weights = self.xy_weights[g_mask]))
        return split_ev

    
    def split_by_time(self, bintime=100, gti = None):
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
    
    def to_lc(self, dt, pi_low = 35, pi_high = 1909, centroid=None, radius = None, tstart=None, tseg=None):
        if tstart is None and self.gti is not None:
            tstart = self.gti[0][0]
            tseg = self.gti[-1][1] - tstart
        
        reg_mask = np.ones(np.shape(self.time)).astype(bool)
        if (centroid is not None) and (radius is not None):
            reg_mask = (np.sqrt(np.square(self.x-centroid[0]) + np.square(self.y-centroid[1])) < radius).astype(bool)

        
        pi_mask = (self.pi > pi_low) * (self.pi < pi_high)
        tot_mask = reg_mask * pi_mask
        return lc.Lightcurve.make_lightcurve(self.time[tot_mask], dt, tstart=tstart, gti=self.gti, tseg=tseg, mjdref=self.mjdref)


        
        
        
def nuproducts_to_stingray_lc(nuprod_lc_file, rebin = False, buffer = False, rebinsize=1.0, buffersize = 100.0, rebin_method='sum'):
    nuproducts_lc = fits.open(nuprod_lc_file)
    time = nuproducts_lc[1].data['TIME'] + nuproducts_lc[1].header['TSTART']
    dt = time[1]-time[0]
    counts = nuproducts_lc[1].data['RATE'] * dt
    error = nuproducts_lc[1].data['ERROR'] * dt
    mjdref = nuproducts_lc[1].header['MJDREFI'] + nuproducts_lc[1].header['MJDREFF']
    gti = [[x,y] for x,y in nuproducts_lc[2].data]
    if buffer:
        buffered_gti = []
        for x,y in gti:
            if np.abs(y-x) > 2*buffersize:
                buffered_gti.append([x+buffersize, y-buffersize])
        gti = buffered_gti
    if rebin:
        return lc.Lightcurve(time, counts, err=error, gti=gti, mjdref=mjdref).rebin(dt_new=rebinsize, method=rebin_method)
    else:
        return lc.Lightcurve(time, counts, err=error, gti=gti, mjdref=mjdref)
    
    
    
def extract_events(file_A, file_B, buff=0):
    # Extracts events for FPMA and FPMB from .evt files (assumed clean).
    # The gtis are crossed to make sure they are the same.
    
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
    split_ev = []
    gti = events.gti
    for g in gti:
        split_ev.append(ev.EventList(time=events.time, gti=[g], pi = events.pi, \
                       mjdref=events.mjdref))
    return split_ev

def split_ev_by_time(events, bintime=100):
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
    exp_sigma = np.sqrt(2) / np.sqrt(nspec)
    gauss = scipy.stats.norm(0, exp_sigma)
    pval = gauss.sf(power)
    return pval

def sum_lc(lc_1, lc_2):
    common_gti = sting_gti.cross_two_gtis(lc_1.gti, lc_2.gti)
    lc_1.gti = common_gti
    lc_2.gti = common_gti
    lc_1.apply_gtis()
    lc_2.apply_gtis()
    summed_lc = lc.Lightcurve(lc_1.time, lc_1.counts + lc_2.counts, err=np.sqrt(np.square(lc_1.counts_err) + np.square(lc_2.counts_err)), gti=common_gti, mjdref=lc_1.mjdref)


    