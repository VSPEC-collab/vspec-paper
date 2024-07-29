



from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from vspec_vsm.spots import SpotCollection, StarSpot
from vspec_vsm.faculae import FaculaCollection, Facula
from vspec_vsm.granules import Granulation
from vspec_vsm.flares import FlareGenerator
from vspec_vsm.star import Star
from vspec_vsm.config import MSH

import pypsg

import setup_psg;setup_psg.setup_psg()
pypsg.docker.set_url_and_run()

from VSPEC import ObservationModel,PhaseAnalyzer

from cartopy import crs as ccrs

import paths

outfile = paths.figures / 'surface_map_and_lc.pdf'

fig = plt.figure(figsize=(4,6))
ax1 = fig.add_subplot(2,1,1,projection=ccrs.Orthographic(0,90))
ax2 = fig.add_subplot(2,1,2)


CFG_PATH = Path(__file__).parent / 'spot_lightcurve.yaml'

model = ObservationModel.from_yaml(CFG_PATH)
model._build_star()
model._warm_up_star(facula_warmup_time=10*u.hr)



lon0 = 0*u.deg
lat0 = 0*u.deg
pl_radius = 1*u.R_earth
pl_orbit = 0.05*u.AU
inclination = 89.8*u.deg
phase = 180.4*u.deg

model.star.plot_surface(lat0, lon0, ax1)

# LC


model.build_planet()
model.build_spectra()

data = PhaseAnalyzer(model.directories['all_model'])

wl_pixels = [0,300,500,700]
time = data.time.to(u.day)
for i in wl_pixels:
    wl = data.wavelength[i]
    lc = data.lightcurve(
        source='star',
        pixel=i,
        normalize=0
    )
    ax2.plot(time,lc,label=f'{wl:.1f}')
ax2.legend()
ax2.set_xlabel(f'time ({time.unit})')
ax2.set_ylabel('Flux (normalized)')

fig.tight_layout()
fig.savefig(outfile)