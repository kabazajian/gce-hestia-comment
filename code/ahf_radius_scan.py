"""Scan the AHF profile radius: is there ANY bin that reproduces Ref. [1]'s
published corner tables? Runs on nucosmo. Loads each halo once and loops
radii inside, so the cost is one particle pass per halo.

Level/ratio logic copied VERBATIM from figure_render.py -- if these drift the
comparison is meaningless.
"""
import json, sys, time
import numpy as np
sys.path.insert(0, '/Volumes/Data2/hestia_muru')
import stage3_views as s3
from scipy.signal import fftconvolve

D_GC, MURU_SMAX, H = 8.0, 15.0, 0.677
X = np.linspace(-19.95, 19.95, 400)
ST_LEVELS = (0.7, 0.5, 0.3)
TARGETS = [2., 3., 4., 5., 6.77, 8., 10., 12., 15., 20., 25., 30.]
NAME = {"G12":"G1.1","G11":"G1.2","G22":"G2.1","G21":"G2.2","G32":"G3.1","G31":"G3.2"}

def region_axis_ratio(m, level):
    ib, il = np.where(m >= level * m.max())
    if len(ib) < 8: return float('nan')
    c = np.column_stack([X[il], X[ib]]); c = c - c.mean(axis=0)
    s = np.linalg.svd(c, compute_uv=False)
    return float(s[1] / s[0])

def b0_profile(m): return 0.5 * (m[199] + m[200])

def matched_levels(m, st):
    ps, pm = b0_profile(st/st.max()), b0_profile(m/m.max())
    out = []
    for L in ST_LEVELS:
        idx = np.where(ps >= L)[0].max()
        out.append(float(pm[idx]))
    return out

# AHF profiles, keyed by our halo name, converged rows only
PROF = {}
for tag, halo in NAME.items():
    d = np.loadtxt(f'/Volumes/Data2/hestia_muru/profiles/profile_{tag}_AHF.txt', comments='#')
    PROF[halo] = d[d[:, 0] > 0]

def frame_at(halo, target):
    d = PROF[halo]
    i = int(np.argmin(np.abs(d[:, 0] - target))); row = d[i]
    Ea, Ec = row[13:16], row[19:22]
    v_n = Ec/np.linalg.norm(Ec)
    v_maj = Ea - (Ea @ v_n)*v_n
    return v_n, v_maj/np.linalg.norm(v_maj), float(row[0]), float(row[11])

kern = s3.disk_kernel()
pad = kern.shape[0]//2
ext = np.linspace(-20-pad*s3.STEP, 20+pad*s3.STEP, s3.NPIX+2*pad+1)
res = {}
for name in s3.FILES:
    t0 = time.time()
    pos_c, mass, rho, star, old, z_hat, x_hat, fr = s3.load_halo(name)
    off = np.array(json.load(open('/Volumes/Data2/hestia_muru/muru_frames.json'))[name]["center_refined_offset_kpc"])
    pos_c = pos_c + off
    dm_sel = rho > 0
    res[name] = {}
    for target in TARGETS:
        v_n, maj, r_bin, b_a = frame_at(name, target)
        l, b, s = s3.lbs(pos_c, D_GC*maj, v_n)
        def conemass(sel):
            m = (s <= MURU_SMAX) & sel & (np.abs(l) < 23+s3.STEP) & (np.abs(b) < 23+s3.STEP)
            h, _, _ = np.histogram2d(b[m], l[m], bins=[ext, ext], weights=mass[m])
            return fftconvolve(h, kern, mode="same")[pad:-pad, pad:-pad]
        mu, st = conemass(dm_sel), conemass(old)
        lv = matched_levels(mu, st)
        res[name][f"{target}"] = {
            "r_bin_kpch": r_bin, "b_over_a": b_a,
            "dm_levels": lv,
            "dm_ratios": [region_axis_ratio(mu/mu.max(), x) for x in lv],
            "star_ratios": [region_axis_ratio(st/st.max(), x) for x in ST_LEVELS]}
        del l, b, s, mu, st
    print(f"{name}: {len(TARGETS)} radii in {time.time()-t0:.0f} s", flush=True)
    del pos_c, mass, rho, star, old
json.dump(res, open('/Volumes/Data2/hestia_muru/radius_scan.json','w'), indent=1)
print("wrote radius_scan.json")
