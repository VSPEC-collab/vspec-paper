"""
Flare frequency and the lightcurve.
"""

from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt

from vspec_vsm.flares import FlareGenerator
from VSPEC import ObservationModel, PhaseAnalyzer
import pypsg

import paths

pypsg.docker.set_url_and_run()

SEED = 10
outfile = paths.figures / 'flare_freq.pdf'
cfg_path = paths.scripts / 'flare_lc.yaml'

dt = 10000*u.day  # a long time.
gen = FlareGenerator(
    dist_teff_mean=9000*u.K,
    dist_teff_sigma=1000*u.K,
    dist_fwhm_mean=3*u.hr,
    dist_fwhm_logsigma=0.3,
    min_energy=1e33*u.erg,
    cluster_size=4,
    rng=np.random.default_rng(seed=SEED)
)
Es = np.logspace(np.log10(gen.min_energy.to_value(u.erg)),
                 np.log10(gen.min_energy.to_value(u.erg))+4)*u.erg

flares = gen.generate_flare_series(dt)

energies = np.array([flare.energy.to_value(u.erg) for flare in flares])

energies_ge_E = np.array([np.sum(energies >= E) for E in Es.to_value(u.erg)])

measured_freq = energies_ge_E/dt
measured_freq_err = np.sqrt(energies_ge_E)/dt

beta = gen.beta
alpha = gen.alpha

expected_log_freq = beta + alpha*np.log10(Es/u.erg)
expected_freq = 10**expected_log_freq / u.day

ratio = np.where(energies_ge_E > 0, measured_freq/expected_freq, np.nan)
ratio_err = np.where(
    energies_ge_E > 0, measured_freq_err/expected_freq, np.nan)

fig, axes = plt.subplots(2, 1,figsize=(4,5),height_ratios=[1,2])

axes[0].plot(Es, expected_freq, c='xkcd:azure', label='Expected')
axes[0].errorbar(Es, measured_freq, yerr=measured_freq_err, fmt='o',
                 color='xkcd:rose pink', label='Observed', markersize=3)
axes[0].set_xlabel('Energy (erg)')
axes[0].set_ylabel('Frequency (1/day)')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].legend()

model = ObservationModel.from_yaml(cfg_path)
model.build_planet()
model.build_spectra()

data = PhaseAnalyzer(model.directories['all_model'])
wl_pixels = [0,20,50,120]
time = data.time.to(u.day)
for i in wl_pixels:
    wl = data.wavelength[i]
    lc = data.lightcurve(
        source='star',
        pixel=i,
        normalize=0
    )
    axes[1].plot(time,lc,label=f'{wl:.1f}')
axes[1].legend()
axes[1].set_xlabel(f'time ({time.unit})')
axes[1].set_ylabel('Flux (normalized)')

fig.tight_layout()

fig.savefig(outfile)
