"""
recreate Figure 1 from Cowan & Agol (2011)
"""

import numpy as np
import matplotlib.pyplot as plt
from VSPEC.gcm import heat_transfer as ht

import paths

OUTFIG = paths.figures / 'cowan_fig1.pdf'

eps = [1e-4,0.2*np.pi,2*np.pi]
label = ['10^{-4}','2\\pi/10','2\\pi']
modes = ['ivp_reflect','ivp_reflect','bvp']
colors = ['b','k','r']

n_points = 100

fig,axes = plt.subplots(2,1,figsize=(4,6))
ax = axes[0]

for e,l,m,c in zip(eps,label,modes,colors):
    lons, tsurf = ht.get_equator_curve(e,n_points,m)
    index = np.argwhere(lons>=-0.5*np.pi)[0][0]
    lons = np.concatenate([lons[index:],lons[:index]+2*np.pi])
    tsurf = np.concatenate([tsurf[index:],tsurf[:index]])
    ax.plot(lons,tsurf,color=c,label=f'$\\epsilon = {l}$')

ax.set_xlabel('Longitude')
ax.set_ylabel('$T/T_0$')
ax.set_xticks([-0.5*np.pi,0,0.5*np.pi,np.pi, 1.5*np.pi,])
ax.set_xticklabels(['$-\\pi/2$','$0$','$\\pi/2$','$\\pi$','$3\\pi/2$'])
_=ax.legend()

ax:plt.Axes = axes[1]
alphas = [0.0,0.5,1.0]
labels = ['0','\\frac{1}{2}','1']
colors = ['xkcd:azure', 'xkcd:lavender', 'xkcd:purple']

def get_curve(alpha:float, lat:np.ndarray) -> np.ndarray:
    return (np.pi/4 * alpha + (1-alpha)*np.cos(lat))**(0.25)

n_points = 1000
for a,l,c in zip(alphas,labels,colors):
    x = np.linspace(0,np.pi/2,n_points)
    tsurf = get_curve(a,x)
    ax.plot(x,tsurf,label=f'$\\alpha = {l}$',color=c)

ax.set_xlabel('Latitude')
ax.set_ylabel('$T/T_{\\rm eq}$')
ax.set_xticks([0,0.25*np.pi,0.5*np.pi])
ax.set_xticklabels(['$0$','$\\pi/4$','$\\pi/2$'])
_=ax.legend()

fig.tight_layout()
fig.savefig(OUTFIG)
