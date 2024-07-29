"""
Plot the lightcurve of a transit with limb darkening.
"""

import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import VSPEC
import pypsg

import paths

LAMBERT_CONFIG = paths.scripts / 'ld_lc_lamb.yaml'
SOLAR_CONFIG = paths.scripts / 'ld_lc_solar.yaml'
OUTFILE = paths.figures / 'ld_transit.pdf'

pypsg.docker.set_url_and_run()

lam_model = VSPEC.ObservationModel.from_yaml(LAMBERT_CONFIG)
sol_model = VSPEC.ObservationModel.from_yaml(SOLAR_CONFIG)

# lam_model._build_star()
# sol_model._build_star()

lam_model.build_planet()
lam_model.build_spectra()
sol_model.build_planet()
sol_model.build_spectra()

lam_data = VSPEC.PhaseAnalyzer.from_model(lam_model)
sol_data = VSPEC.PhaseAnalyzer.from_model(sol_model)

fig = plt.figure(figsize=(8, 6))
gs = fig.add_gridspec(5, 2)
fig.subplots_adjust(hspace=0,wspace=0.5)
ax1: plt.Axes = fig.add_subplot(gs[1:4,1])

wl_index = np.argmin(np.abs(lam_data.wavelength - 1*u.um))

ax1.plot(lam_data.time.to_value(u.hr), lam_data.lightcurve('total',wl_index,normalize=0), 'r', label='Lambertian')
ax1.scatter(lam_data.time.to_value(u.hr), lam_data.lightcurve('total',wl_index,normalize=0,noise=True), c='r',marker='.',s=1,alpha=1)
ax1.plot(sol_data.time.to_value(u.hr), sol_data.lightcurve('total',wl_index,normalize=0), 'k', label='$u_1=0.5,~u_2=0.15$')
ax1.scatter(sol_data.time.to_value(u.hr), sol_data.lightcurve('total',wl_index,normalize=0,noise=True), c='k', marker='.',s=1,alpha=1)


fig.text(0.87,0.75,f'$\\lambda={lam_data.wavelength[wl_index].value:.0f} \\mu\\mathrm{{m}}$',fontsize=12)


ax1.legend(loc=(0.1,1.1))


ax1.set_xlabel('Time (hr)')
ax1.set_ylabel('Normalized Flux')


lam_ax = fig.add_subplot(gs[:2,0], projection=ccrs.Orthographic())
sol_ax = fig.add_subplot(gs[3:,0], projection=ccrs.Orthographic())

lam_model.star.plot_surface(0*u.deg,0*u.deg,ax=lam_ax,
                            orbit_radius=lam_model.params.planet.semimajor_axis,
                            radius=lam_model.params.planet.radius,
                            phase=179.3*u.deg,
                            inclination=lam_model.params.system.inclination
                            )
sol_model.star.plot_surface(0*u.deg,0*u.deg,ax=sol_ax,
                            orbit_radius=sol_model.params.planet.semimajor_axis,
                            radius=sol_model.params.planet.radius,
                            phase=179.3*u.deg,
                            inclination=sol_model.params.system.inclination
                            )

# inax = lam_ax.inset_axes([-0.99,-0.15,0.3,0.3],transform=ccrs.PlateCarree())
# tmask, _ =lam_model.star.get_transit_mask(
#     0*u.deg,0*u.deg,
#     orbit_radius=lam_model.params.planet.semimajor_axis,
#     radius=lam_model.params.planet.radius,
#     phase=179.3*u.deg,
#     inclination=lam_model.params.system.inclination,
# )
# lats,lons,_=lam_model.star.gridmaker.display_grid(500,1000,tmask)
# inax.pcolormesh(lons,lats,tmask.T)



fig.savefig(OUTFILE)