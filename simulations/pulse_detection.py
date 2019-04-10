import matplotlib as mpl
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import stingray.simulator as ss
from stingray.lightcurve import Lightcurve
from stingray.events import EventList
from tqdm import tqdm
import matplotlib as mpl
import time
import seaborn as sns

def pulsar_events_mp(length, period, ctrate, pulsed_fraction, mean_obs, bkg_ctrate, detlev, nbin = 128):
    
    BUFFER_LC = Lightcurve([0, 1], [1, 1], gti=[[-0.5, 1.5]], dt=1, err_dist='gauss')
    nustar_orb = 5808
    
    dt = period / 20
    # The total length of the time series should be the number of pointings times the time per orbit. 
    # Add one orbit for buffer.
    N_orb = int(round(length/mean_obs, 0))
    tot_len = (N_orb + 1)*nustar_orb
    times = numpy.arange(dt/2, tot_len + dt/2, dt)
    cont_lc = numpy.random.poisson((ctrate * (1 + pulsed_fraction * numpy.cos(2 * numpy.pi / period * times)) * dt)) + numpy.random.poisson(bkg_ctrate*dt)
    lc = BUFFER_LC
    lc.time = times
    lc.counts = cont_lc
    # The orbital period is 5808s. Every 5808s, a continuous observation with min_obs < length < max_obs begins
    start_t = numpy.multiply(numpy.arange(N_orb), numpy.random.normal(loc = nustar_orb, scale = 120, size = N_orb))
    point_t = numpy.random.normal(loc=mean_obs, scale = mean_obs/4, size = N_orb)
    end_t = numpy.add(start_t, point_t)
    exposure = numpy.sum(point_t)
    lc.gti = numpy.column_stack((start_t, end_t))
    lc.dt = dt
    events = EventList()
    events.gti = lc.gti
    events.simulate_times(lc)
    phase = numpy.arange(0, 1, 1 / nbin)
    zsq = z_n(phase, n=2,
              norm=fold_events(events.time, 1/period, nbin=nbin)[1])
    detected = zsq > detlev
    return (detected, exposure)


def detected_pulse_fraction_mp(pf_min, pf_max, length_min, length_max, 
                            ctrate_min=1.4, ctrate_max=1.4, 
                            period_min=1, period_max=1, n_realizations=1000, 
                            ntrial=1000, results=None, nbin=128, min_mean_obs = 3000, max_mean_obs = 3250, 
                            bkg_ctrate = 0, cores = 8):
    if results is None:
        results = Table(names=["period", "countrate", "pf", 
                               "length", "mean_obs", "detected"], 
                        dtype=[float, float, float, float, float, bool])

    pfs = 10**np.random.uniform(np.log10(pf_min), 
                                np.log10(pf_max), n_realizations)
    lengths = 10**np.random.uniform(np.log10(length_min), 
                                    np.log10(length_max), n_realizations)
    periods = 10**np.random.uniform(np.log10(period_min), 
                                    np.log10(period_max), n_realizations)
    ctrates = 10**np.random.uniform(np.log10(ctrate_min), 
                                    np.log10(ctrate_max), n_realizations)
    mean_obss = 10**np.random.uniform(np.log10(min_mean_obs), 
                                    np.log10(max_mean_obs), n_realizations)
    
    c = ipp.Client()
    v = c[:]
    with v.sync_imports():
        import numpy
        from stingray.events import EventList
        from stingray.lightcurve import Lightcurve
        from stingray.pulse.pulsar import z2_n_detection_level, z_n, fold_events
        
    detlev = z2_n_detection_level(ntrial=ntrial)
    map_results = v.map(pulsar_events_mp, lengths, periods, ctrates, pfs, mean_obss, [bkg_ctrate for i in range(n_realizations)], [detlev for i in range(n_realizations)])
    c.wait_interactive()
    c.shutdown()
    for i in range(n_realizations):
        pf, period, ctrate, mean_obs = pfs[i], periods[i], ctrates[i], mean_obss[i]
        detected, length = map_results[i]
        results.add_row([period, ctrate, pf, length, mean_obs, detected])
    
    return results

