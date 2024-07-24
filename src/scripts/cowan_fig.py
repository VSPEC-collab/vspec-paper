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

fig,ax = plt.subplots(1,1,figsize=(4,3))

for e,l,m,c in zip(eps,label,modes,colors):
    lons, tsurf = ht.get_equator_curve(e,n_points,m)
    index = np.argwhere(lons>=-0.5*np.pi)[0][0]
    lons = np.concatenate([lons[index:],lons[:index]+2*np.pi])
    tsurf = np.concatenate([tsurf[index:],tsurf[:index]])
    ax.plot(lons,tsurf,color=c,label=f'$\\epsilon = {l}$')

ax.set_xlabel('$\\Phi$')
ax.set_ylabel('$T/T_0$')
ax.set_xticks([-0.5*np.pi,0,0.5*np.pi,np.pi, 1.5*np.pi,])
ax.set_xticklabels(['$-\\pi/2$','$0$','$\\pi/2$','$\\pi$','$3\\pi/2$'])
_=ax.legend()
fig.tight_layout()
fig.savefig(OUTFIG)
