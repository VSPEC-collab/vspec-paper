"""
Compare jax and scipy speed.
"""

import jax.numpy as jnp
import numpy as np
from time import time
import matplotlib.pyplot as plt
from astropy import units as u

from GridPolator import GridSpectra

import paths

OUTFIE = paths.figures / 'jax_v_scipy.pdf'

w1 = 5 * u.um
w2 = 12 * u.um
resolving_power = 100
teffs = [2800,2900,3000,3100,3200,3300]
impl_bin = 'rust'

wl_jnp = jnp.linspace(5.0, 11.2, 100)
wl_np = np.linspace(5.0, 11.2, 100)
params_jnp = (jnp.array([2900.]),)
params_np = (np.array([2900.]),)


fig, ax  = plt.subplots(1,1,figsize=(4,3))

N=3000

g_jax = GridSpectra.from_vspec(
    w1=w1,
    w2=w2,
    resolving_power=resolving_power,
    teffs=teffs,
    impl_bin=impl_bin,
    impl_interp='jax',
    fail_on_missing=False
)
g_scipy = GridSpectra.from_vspec(
    w1=w1,
    w2=w2,
    resolving_power=resolving_power,
    teffs=teffs,
    impl_bin=impl_bin,
    impl_interp='scipy',
    fail_on_missing=False
)

dt_jax = np.zeros(N)

for i in range(N):
    start = time()
    flux_jnp = g_jax.evaluate(params_jnp, wl_jnp)
    end = time()
    dt_jax[i] = end - start

dt_scipy = np.zeros(N)

for i in range(N):
    start = time()
    flux_np = g_scipy.evaluate(params_np, wl_np)
    end = time()
    dt_scipy[i] = end - start

x = np.arange(N)

ax.plot(x, np.cumsum(dt_jax), label='JAX',c='#B96EBD')
ax.plot(x, np.cumsum(dt_scipy), label='Scipy',c='#0054A6')
ax.set_xlabel('Iteration')
ax.set_ylabel('Time (s)')
fig.tight_layout()
_=ax.legend()
fig.savefig(OUTFIE)