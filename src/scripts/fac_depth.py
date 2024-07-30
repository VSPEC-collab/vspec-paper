import numpy as np
import matplotlib.pyplot as plt

import paths

from vspec_vsm.helpers import calc_circ_fraction_inside_unit_circle

OUTFILE = paths.figures / 'fac_depth.png'

h=2
alpha = np.pi/4 * 0.1


x = np.linspace(-2,2,1000)
y = np.linspace(-2,2,1000)

X,Y = np.meshgrid(x,y)

R = np.sqrt((X/np.cos(alpha))**2 + Y**2)
THETA = np.arctan2(Y,X)

INSIDE = R<1

DIST_FROM_RIGHT = (np.cos(alpha)*np.sqrt(np.abs(1-Y**2))-X)

Heff = np.sin(alpha) * h
WALL = (DIST_FROM_RIGHT<Heff) & INSIDE

Z = np.ones_like(X) * np.nan
Z = np.where(INSIDE,0,Z)
Z = np.where(WALL,1,Z)

wall_frac = np.sum(WALL)/np.sum(INSIDE)
circ_frac = (1 - calc_circ_fraction_inside_unit_circle(h*np.tan(alpha),0,1))


fig,ax = plt.subplots(1,1,figsize=(4,3))
ax.set_aspect('equal')

im = ax.pcolormesh(X/np.cos(alpha),Y,Z,vmin=0,vmax=1)

ax.text(0.3,1.1,f'$f_{{\\rm wall}} = {wall_frac:.4f}$')
ax.text(0.3,1.4,f'$f_{{\\rm circ}} = {circ_frac:.4f}$')

fig.colorbar(im)
fig.savefig(OUTFILE,dpi=300)

