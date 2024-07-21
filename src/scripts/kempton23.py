from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pypsg

from VSPEC import ObservationModel,PhaseAnalyzer
from VSPEC import params
from VSPEC.config import MSH

import paths

SEED = 1214
import setup_psg;setup_psg.setup_psg()
pypsg.docker.set_url_and_run()

outfile = paths.figures / 'kempton23.pdf'

# Instrument
inst = params.InstrumentParameters.miri_lrs()

# Observation
observation = params.ObservationParameters(
    observation_time=41.0*u.hr,
    integration_time=15*u.min
)

# PSG
psg_params = params.psgParameters(
    gcm_binning=200,
    phase_binning=1,
    use_continuum_stellar=True,
    nmax=0,
    lmax=0,
    continuum=['Rayleigh', 'Refraction', 'CIA_all'],
    use_molecular_signatures=True
)

# Star and Planet
star_teff = 3250*u.K
star_rad = 0.215*u.R_sun
planet_rad = 2.742*u.R_earth
orbit_rad = 0.01490*u.AU
orbit_period = 1.58040433*u.day
planet_rot_period = orbit_period
star_rot_period = 120 * u.day
planet_mass = 8.17*u.M_earth
star_mass = 0.178*u.M_sun
inclination = 88.7*u.deg

start_time_before_eclipse = 2*u.hr
angle_before_eclipse = (2*np.pi*u.rad * start_time_before_eclipse/orbit_period).to(u.deg)
initial_phase = 0*u.deg - angle_before_eclipse

planet_params = params.PlanetParameters(
    name='GJ1214b',
    radius=planet_rad,
    gravity=params.GravityParameters('kg',planet_mass),
    semimajor_axis=orbit_rad,
    orbit_period=orbit_period,
    rotation_period=planet_rot_period,
    eccentricity=0,
    obliquity=0*u.deg,
    obliquity_direction=0*u.deg,
    init_phase=initial_phase,
    init_substellar_lon=0*u.deg
)

system_params = params.SystemParameters(
    distance=14.6427*u.pc,
    inclination=inclination,
    phase_of_periasteron=0*u.deg
)

star_dict = {
    'teff': star_teff,
    'radius': star_rad
}
planet_dict = {'semimajor_axis': orbit_rad}

gcm_dict = {
    'nlayer': 30,
    'nlon': 30,
    'nlat': 15,
    'epsilon': 6,
    'albedo': 0.3,
    'emissivity': 1.0,
    'lat_redistribution': 0.1,
    'gamma': 1.4,
    'psurf': 1*u.bar,
    'ptop': 1e-5*u.bar,
    'wind': {'U': '0 m/s','V':'0 m/s'},
    'molecules':{'CO2':0.99}
}

gcm = params.gcmParameters.from_dict({
    'star':star_dict,
    'planet':planet_dict,
    'gcm':{'vspec':gcm_dict,'mean_molec_weight':44}
})

star_kwargs = dict(
    psg_star_template='M',
    teff=star_teff,
    mass=star_mass,
    radius=star_rad,
    period=star_rot_period,
    misalignment=0*u.deg,
    misalignment_dir=0*u.deg,
    ld=params.LimbDarkeningParameters.proxima(),
    faculae=params.FaculaParameters.none(),
    flares=params.FlareParameters.none(),
    granulation=params.GranulationParameters.none(),
    grid_params=(500,1000),
    spectral_grid='default'
)

quiet_star = params.StarParameters(
    spots=params.SpotParameters.none(),
    **star_kwargs
)
spotted_star = params.StarParameters(
    spots=params.SpotParameters(
        distribution='iso',
        initial_coverage=0.2,
        area_mean=300*MSH,
        area_logsigma=0.2,
        teff_umbra=2700*u.K,
        teff_penumbra=2700*u.K,
        equillibrium_coverage=0.2,
        burn_in=0*u.s,
        growth_rate=0.0/u.day,
        decay_rate=0*MSH/u.day,
        initial_area=10*MSH
    ),
    **star_kwargs
)

# Set parameters for simulation
header_kwargs = dict(
    teff_min=2300*u.K,teff_max=3400*u.K,
    seed = SEED,
    verbose = 0
)
internal_params_kwargs = dict(
    planet=planet_params,
    system=system_params,
    obs=observation,
    gcm=gcm,
    psg=psg_params,
    inst=inst
)

params_quiet = params.InternalParameters(
    header=params.Header(data_path=Path('.vspec/gj1214_quiet'),**header_kwargs),
    star = quiet_star,
    **internal_params_kwargs
)

params_spotted = params.InternalParameters(
    header=params.Header(data_path=Path('.vspec/gj1214_spotted'),**header_kwargs),
    star = spotted_star,
    **internal_params_kwargs
)

model_quiet = ObservationModel(params=params_quiet)
model_quiet.build_planet()
model_quiet.build_spectra()

model_spotted = ObservationModel(params_spotted)
model_spotted.build_planet()
model_spotted.build_spectra()

data_spotted = PhaseAnalyzer(model_spotted.directories['all_model'])

data_quiet = PhaseAnalyzer(model_quiet.directories['all_model'])
flux_unit = u.Unit('W m-2 um-1')
def get_star(data:PhaseAnalyzer):
    i_eclipse1 = np.argmin(data.lightcurve('total',(0,-1))[:data.n_images//4])
    i_eclipse2 = np.argmin(data.lightcurve('total',(0,-1))[3*data.n_images//4:]) + 3*data.n_images//4
    time = (data.time-data.time[0]).to_value(u.hr)
    star_spec1 = data.spectrum('total',i_eclipse1).to_value(flux_unit)
    star_spec2 = data.spectrum('total',i_eclipse2).to_value(flux_unit)

    def func(t:float):
        m = (star_spec2 - star_spec1)/(time[i_eclipse2]-time[i_eclipse1])
        x = t-time[i_eclipse1]
        b = star_spec1
        y = m * x + b
        return y

    return func


def get_lc(data:PhaseAnalyzer,w1:u.Quantity,w2:u.Quantity):
    """Get the lightcurve given a bandpass"""
    wl = data.wavelength
    i_left = int(np.argwhere(wl > w1)[0])
    try:
        i_right = int(np.argwhere(wl > w2)[0])
    except IndexError:
        i_right = -1
    interp = get_star(data)
    ts = (data.time-data.time[0]).to_value(u.hr)
    star_im = np.array([interp(t) for t in ts]).T
    total_im = data.total.to_value(flux_unit)
    pl_im = total_im-star_im
    lc = 1e6*pl_im[i_left:i_right,:]/star_im[i_left:i_right,:]
    return lc.mean(axis=0)

bin_edges = np.arange(5.0,12.0,0.5)
n_ax = len(bin_edges)
fig,axes = plt.subplots(n_ax,1,figsize=(4,10),sharex=True)
for edge,ax in zip(bin_edges,axes):
    w1,w2 = edge*u.um,(edge+0.5)*u.um
    quiet_lc = get_lc(data_quiet,w1,w2)
    spotted_lc = get_lc(data_spotted,w1,w2)
    time = (data_quiet.time - data_quiet.time[0]).to(u.hr)
    ax.plot(time,(quiet_lc),c='xkcd:azure',label='No Spots')
    ax.plot(time,(spotted_lc),c='xkcd:lavender',label='Spotted')
    ax.text(0.55,0.7,f'{w1:.1f} - {w2:.1f}',transform=ax.transAxes)
    ax.set_ylim(-100,700)
fig.subplots_adjust(hspace=0,wspace=0)
axes[0].legend(loc='upper left')
axes[-1].set_xlabel('Time (hour)')
_ = axes[n_ax//2].set_ylabel('Planet Flux (ppm)')

fig.savefig(outfile)