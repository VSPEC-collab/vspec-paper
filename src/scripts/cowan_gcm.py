"""
https://vspec.readthedocs.io/en/latest/auto_examples/other/gcm/plot_cowan2011.html#sphx-glr-auto-examples-other-gcm-plot-cowan2011-py
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
import cartopy.crs as ccrs

import VSPEC.gcm.heat_transfer as ht

import paths

outfile = paths.figures / 'cowan_gcm.pdf'

epsilon = 2*np.pi
star_teff = 5800*u.K
albedo = 0.3
r_star = 1*u.R_sun
r_orbit = 1*u.AU

tmap = ht.TemperatureMap.from_planet(
    epsilon=epsilon,
    star_teff=star_teff,
    albedo=albedo,
    r_star=r_star,
    r_orbit=r_orbit
)
lons = np.linspace(-180,180,90,endpoint=False)*u.deg
lats = np.linspace(-90,90,46,endpoint=True)*u.deg

longrid,latgrid = np.meshgrid(lons,lats)
data = tmap.eval(lon=longrid,lat=latgrid,alpha=0)

fig = plt.figure(figsize=(4,3))
proj = ccrs.Robinson(central_longitude=0)
ax = fig.add_subplot(projection=proj)

im = ax.pcolormesh(lons.to_value(u.deg),lats.to_value(u.deg),data.to_value(u.K),cmap='gist_heat',transform=ccrs.PlateCarree(),rasterized=True)
gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
    color='grey', alpha=0.8, linestyle='--')
gl.top_xlabels = False
gl.right_ylabels = False


_=fig.colorbar(im,ax=ax,label='Surface Temperature (K)',pad=0.15,shrink=0.6)

fig.tight_layout()

fig.savefig(outfile)