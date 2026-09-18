"""Azimuth scan: does ANY in-plane observer direction reproduce Ref. [1]'s
published corner tables? Runs on nucosmo; one particle pass per halo.

Uses THEIR observer-placement formula exactly (GalacticCoordinates.jl,
generate_sun_position):  p_obs = p_c + d * (cos a * v_d + sin a * (v_n x v_d)),
with v_n = AHF Ec and v_d = AHF Ea at the bin nearest 10 kpc/h (their figure
path hard-codes a = 0). If some a reproduces all three published ratio
triplets to ~0.02, their v_d was not Ea and we learn where it pointed. If no
a does, the viewpoint is excluded altogether and the residual lives in the
map construction or particle selection.

Level/ratio logic copied verbatim from figure_render.py.
"""
import json, sys, time
import numpy as np
sys.path.insert(0, '/Volumes/Data2/hestia_muru')
import stage3_views as s3
from scipy.signal import fftconvolve

D_GC, MURU_SMAX = 8.0, 15.0
X = np.linspace(-19.95, 19.95, 400)
ST_LEVELS = (0.7, 0.5, 0.3)
ALPHAS = list(range(0, 360, 10))
NAME = {"G12":"G1.1","G11":"G1.2","G22":"G2.1","G21":"G2.2","G32":"G3.1","G31":"G3.2"}  # swapped pairs

def region_axis_ratio(m, level):
    ib, il = np.where(m >= level * m.max())
    if len(ib) < 8: return float('nan')
    c = np.column_stack([X[il], X[ib]]); c = c - c.mean(axis=0)
    s = np.linalg.svd(c, compute_uv=False)
    return float(s[1] / s[0])
def b0_profile(m): return 0.5 * (m[199] + m[200])
def matched_levels(m, st):
    ps, pm = b0_profile(st/st.max()), b0_profile(m/m.max())
    return [float(pm[np.where(ps >= L)[0].max()]) for L in ST_LEVELS]

PROF = {}
for tag, halo in NAME.items():
    d = np.loadtxt(f'/Volumes/Data2/hestia_muru/profiles/profile_{tag}_AHF.txt', comments='#')
    PROF[halo] = d[d[:, 0] > 0]
def base_frame(halo, target=10.0):
    d = PROF[halo]; row = d[int(np.argmin(np.abs(d[:, 0] - target)))]
    Ea, Ec = row[13:16], row[19:22]
    v_n = Ec/np.linalg.norm(Ec)
    v_d = Ea - (Ea @ v_n)*v_n; v_d /= np.linalg.norm(v_d)
    return v_n, v_d

kern = s3.disk_kernel(); pad = kern.shape[0]//2
ext = np.linspace(-20-pad*s3.STEP, 20+pad*s3.STEP, s3.NPIX+2*pad+1)
frames = json.load(open('/Volumes/Data2/hestia_muru/muru_frames.json'))
res = {}
for name in s3.FILES:
    t0 = time.time()
    pos_c, mass, rho, star, old, z_hat, x_hat, fr = s3.load_halo(name)
    pos_c = pos_c + np.array(frames[name]["center_refined_offset_kpc"])
    dm_sel = rho > 0
    v_n, v_d = base_frame(name)
    v_p = np.cross(v_n, v_d)                       # their sin(a) direction
    res[name] = {}
    for a in ALPHAS:
        ar = np.deg2rad(a)
        v3 = np.cos(ar)*v_d + np.sin(ar)*v_p; v3 /= np.linalg.norm(v3)
        l, b, s = s3.lbs(pos_c, D_GC*v3, v_n)
        def conemass(sel):
            m = (s <= MURU_SMAX) & sel & (np.abs(l) < 23+s3.STEP) & (np.abs(b) < 23+s3.STEP)
            h, _, _ = np.histogram2d(b[m], l[m], bins=[ext, ext], weights=mass[m])
            return fftconvolve(h, kern, mode="same")[pad:-pad, pad:-pad]
        mu, st = conemass(dm_sel), conemass(old)
        lv = matched_levels(mu, st)
        res[name][str(a)] = {"dm_levels": lv,
                             "dm_ratios": [region_axis_ratio(mu/mu.max(), x) for x in lv],
                             "star_ratios": [region_axis_ratio(st/st.max(), x) for x in ST_LEVELS]}
        del l, b, s, mu, st
    print(f"{name}: {len(ALPHAS)} azimuths in {time.time()-t0:.0f} s", flush=True)
    del pos_c, mass, rho, star, old
json.dump(res, open('/Volumes/Data2/hestia_muru/azimuth_scan.json', 'w'), indent=1)
print("wrote azimuth_scan.json")
