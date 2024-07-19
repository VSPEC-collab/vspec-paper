from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pypsg

from VSPEC import ObservationModel,PhaseAnalyzer
from VSPEC import params
from VSPEC.config import MSH

import paths

outfile = paths.figures / 'moran23.pdf'

SEED = 10
PRE_TRANSIT = 8
IN_TRANSIT = 13
import setup_psg;setup_psg.setup_psg()
pypsg.docker.set_url_and_run()

# Instrument

inst = params.InstrumentParameters(
    telescope=params.SingleDishParameters.jwst(),
    bandpass=params.BandpassParameters(
        wl_blue=2.87*u.um,
        wl_red=5.14*u.um,
        resolving_power=200,
        wavelength_unit=u.um,
        flux_unit=u.Unit('W m-2 um-1')
    ),
    detector=params.DetectorParameters(
        beam_width=5*u.arcsec,
        integration_time=0.5*u.s,
        ccd=params.ccdParameters(
            pixel_sampling=8,
            read_noise=16.8*u.electron,
            dark_current=0.005*u.electron/u.s,
            throughput=0.3,
            emissivity=0.1,
            temperature=50*u.K
        )
    )
)

# Observation

observation = params.ObservationParameters(
    observation_time=3.53*u.hr,
    integration_time=8*u.min
)

# PSG

psg_kwargs = dict(
    gcm_binning=200,
    phase_binning=1,
    nmax=0,
    lmax=0,
    use_continuum_stellar=True,
    continuum=['Rayleigh', 'Refraction'],
)
psg_params = params.psgParameters(
    use_molecular_signatures=True,
    **psg_kwargs
)
psg_no_atm = params.psgParameters(
    use_molecular_signatures=False,
    **psg_kwargs
)

# Star and Planet

star_teff = 3343*u.K
star_rad = 0.339*u.R_sun
a_rstar = 11.229 # Moran+23 Table 1
rp_rstar = 0.03709 # Moran+23 Table 1
planet_rad = star_rad * rp_rstar
orbit_rad = star_rad * a_rstar
orbit_period = 1.467119*u.day
planet_rot_period = orbit_period
star_rot_period = np.inf*u.s # assume the star does not change.
planet_mass = 2.82*u.M_earth
star_mass = 0.323*u.M_sun
inclination = 89.06*u.deg
planet_norm_factor = {
    'rock_quiet': 1,
    'rock_spotted': 0.97,
    'water_quiet': 0.98,
} # We will multiply the radius by these scalars rather than fit
  # later (because I have run this code before and know how much
  # each model is off by).

print(f'Planet radius: {planet_rad.to_value(u.R_earth)} R_earth')
print(f'Semimajor axis: {orbit_rad.to_value(u.AU)} AU')

observation_angle = (2*np.pi*u.rad * observation.observation_time/orbit_period).to(u.deg)
initial_phase = 180*u.deg - 0.5*observation_angle

planet_kwargs = dict(
    name='GJ486b',
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

quiet_rock_planet = params.PlanetParameters(
    radius=planet_rad*planet_norm_factor['rock_quiet'],
    **planet_kwargs
)
quiet_water_planet = params.PlanetParameters(
    radius=planet_rad*planet_norm_factor['water_quiet'],
    **planet_kwargs
)
spotted_rock_planet = params.PlanetParameters(
    radius=planet_rad*planet_norm_factor['rock_spotted'],
    **planet_kwargs
)


system_params = params.SystemParameters(
    distance=8.07*u.pc,
    inclination=inclination,
    phase_of_periasteron=0*u.deg
)


star_dict = {
    'teff': star_teff,
    'radius': star_rad
}
planet_dict = {'semimajor_axis': orbit_rad}

gcm_dict = {
    'nlayer': 40,
    'nlon': 60,
    'nlat': 30,
    'epsilon': 100,
    'albedo': 0.3,
    'emissivity': 1.0,
    'lat_redistribution': 0.0,
    'gamma': 1.0,
    'psurf': 1*u.bar,
    'ptop': 1e-10*u.bar,
    'wind': {'U': '0 m/s','V':'0 m/s'},
    'molecules':{'H2O':0.99}
}

# Create two sets of GCM Parameters

h2o_atm = {'molecules':{'H2O':0.99}}
gcm_h2o = params.gcmParameters.from_dict({
    'star':star_dict,
    'planet':planet_dict,
    'gcm':{'vspec':dict(gcm_dict,**h2o_atm),'mean_molec_weight':18}
})
star_kwargs = dict(
    psg_star_template='M',
    teff=star_teff,
    mass=star_mass,
    radius=star_rad,
    period=star_rot_period,
    misalignment=0*u.deg,
    misalignment_dir=0*u.deg,
    ld=params.LimbDarkeningParameters.lambertian(),
    faculae=params.FaculaParameters.none(),
    flares=params.FlareParameters.none(),
    granulation=params.GranulationParameters.none(),
    grid_params=100000,
    spectral_grid='default'
)
quiet_star = params.StarParameters(
    spots=params.SpotParameters.none(),
    **star_kwargs
)
spotted_star = params.StarParameters(
    spots=params.SpotParameters(
        distribution='iso',
        initial_coverage=0.11,
        area_mean=500*MSH,
        area_logsigma=0.2,
        teff_umbra=2700*u.K,
        teff_penumbra=2700*u.K,
        equillibrium_coverage=0.0,
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
    seed = SEED
)
internal_params_kwargs = dict(
    system=system_params,
    obs=observation,
    gcm=gcm_h2o,
    inst=inst
)

# Make the three cases

params_rock_quiet = params.InternalParameters(
    header=params.Header(data_path=Path('.vspec/rock_quiet'),**header_kwargs),
    star = quiet_star,
    psg=psg_no_atm,
    planet=quiet_rock_planet,
    **internal_params_kwargs
)
params_h2o_quiet = params.InternalParameters(
    header=params.Header(data_path=Path('.vspec/h2o_quiet'),**header_kwargs),
    star = quiet_star,
    psg=psg_params,
    planet=quiet_water_planet,
    **internal_params_kwargs
)

params_rock_spotted = params.InternalParameters(
    header=params.Header(data_path=Path('.vspec/rock_spotted'),**header_kwargs),
    star = spotted_star,
    psg=psg_no_atm,
    planet=spotted_rock_planet,
    **internal_params_kwargs
)

model_rock_quiet = ObservationModel(params_rock_quiet)
model_rock_quiet.build_planet()
model_rock_quiet.build_spectra()
data_rock_quiet = PhaseAnalyzer(model_rock_quiet.directories['all_model'])

model_h2o_quiet = ObservationModel(params_h2o_quiet)
model_h2o_quiet.build_planet()
model_h2o_quiet.build_spectra()
data_h2o_quiet = PhaseAnalyzer(model_h2o_quiet.directories['all_model'])

model_rock_spotted = ObservationModel(params_rock_spotted)
model_rock_spotted.build_planet()
model_rock_spotted.build_spectra()
data_rock_spotted = PhaseAnalyzer(model_rock_spotted.directories['all_model'])

import pandas as pd

filename = paths.static / 'moran2023_fig3.txt'

df = pd.read_fwf(filename,colspecs=[(0,8),(9,14),(15,20),(21,25),(26,28)],
    header=20,names=['Reduction','Wave','Width','Depth','e_Depth'])
used_eureka = df['Reduction']=='Eureka'
moran_x = df.loc[used_eureka,'Wave'].values
moran_y = df.loc[used_eureka,'Depth'].values
moran_yerr = df.loc[used_eureka,'e_Depth'].values
moran_mean = np.mean(moran_y)

fig, ax = plt.subplots(1,1,figsize=(4,3),tight_layout=True)

for data,label,color in zip(
    [data_h2o_quiet,data_rock_spotted],
    ['H2O', 'Rock+Spots'],
    ['xkcd:azure','xkcd:golden yellow']
):
    pre_transit = PRE_TRANSIT
    in_transit = IN_TRANSIT
    unocculted_spectrum = data.spectrum('total',pre_transit)
    occulted_spectrum = data.spectrum('total',in_transit)
    lost_to_transit = unocculted_spectrum-occulted_spectrum
    transit_depth = (lost_to_transit/unocculted_spectrum).to_value(u.dimensionless_unscaled)
    ax.plot(data.wavelength,transit_depth*1e6,label=label,color=color)

ax.errorbar(moran_x,moran_y,yerr=moran_yerr,
    fmt='o',color='k',label='Moran+23 Eureka!',markersize=3)

ax.set_xlabel('Wavelength (um)')
ax.set_ylabel('Transit depth (ppm)')
ax.legend()

fig.savefig(outfile)