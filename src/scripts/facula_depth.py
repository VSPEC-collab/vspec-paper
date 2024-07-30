"""
Facula depth and the lightcurve.

https://vspec-vsm.readthedocs.io/en/latest/auto_examples/plot_effect_facula_depth.html#sphx-glr-auto-examples-plot-effect-facula-depth-py
"""

from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import shutil

from vspec_vsm.faculae import Facula
from vspec_vsm.helpers import round_teff


import paths

fig,ax = plt.subplots(1,1,figsize=(4,3))

outfig = paths.figures / 'facula_depth.pdf'


def make_fac(_depth_over_rad: float) -> Facula:
    """Return a default facula"""
    radius = 0.005*u.R_sun
    return Facula(
        lat=0*u.deg,
        lon=0*u.deg,
        r_max=radius,
        r_init=radius,
        depth=radius*_depth_over_rad,
        # None of the below parameters will affect this
        # example, but we must set them to something
        lifetime=1*u.day,
        floor_teff_slope=0*u.K/u.km,
        floor_teff_min_rad=1*u.km,
        floor_teff_base_dteff=-100*u.K,
        wall_teff_intercept=100*u.K,
        wall_teff_slope=0*u.K/u.km
    )


rad_star = 0.15*u.R_sun
WALL_BRIGHTNESS = 1.1
FLOOR_BRIGNTNESS = 0.9


def relative_flux(_facula: Facula, angle: u.Quantity) -> float:
    """Get the contrast from a toy flux model"""
    effective_area = _facula.effective_area(angle)
    area_floor = effective_area[round_teff(_facula.floor_dteff)]
    area_wall = effective_area[round_teff(_facula.wall_dteff)]
    area_of_disk = np.pi*rad_star**2
    floor_fraction = (
        area_floor/area_of_disk).to_value(u.dimensionless_unscaled)
    wall_fraction = (area_wall/area_of_disk).to_value(u.dimensionless_unscaled)
    return 1 + floor_fraction*(FLOOR_BRIGNTNESS - 1) + wall_fraction*(WALL_BRIGHTNESS-1)

angles = np.linspace(0, 90, 50)*u.deg
depth_over_rad = np.logspace(-1, 1, 5)

for dor in depth_over_rad:
    facula = make_fac(dor)
    flux = np.array([
        relative_flux(facula, angle) for angle in angles
    ])
    x = np.concatenate([np.flip(-angles), angles])
    y = (np.concatenate([np.flip(flux), flux]) - 1)*1e6
    log_dor = np.log10(dor)
    color = cm.get_cmap('viridis')(0.5*(log_dor+1))
    ax.plot(x, y, label=f'$\\log{{\\frac{{D}}{{R}}}} = {log_dor:.1f}$', c=color)
ax.set_xlabel('angle from disk center (deg)')
ax.set_ylabel('Relative flux (ppm)')
_ = ax.legend(loc='lower right',prop={'size':8})

fig.tight_layout()

fig.savefig(outfig)

# Copy the other figure now

static_fig_name = 'facula_diagram.pdf'
static_fig_path = paths.static / static_fig_name
target_path = paths.figures / static_fig_name

shutil.copy(static_fig_path, target_path)
