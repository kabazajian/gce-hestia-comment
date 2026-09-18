"""Viewing frames from Muru's AHF profile files (delivered 2026-09-15).

WHY. Ref. [1]'s App. A places the observer on the major axis of the inertia
tensor of particles within 10 h^-1 kpc. Their code takes that axis from the
AHF *profile* file via `get_galvectors_from_profile()`, which lives in the
private HestiaUtils package. Until 2026-09-15 we approximated it by computing
the tensor ourselves from the particles (`figure_maps.their_frame`); Muru has
now sent the profile files, so the axes can be read rather than reconstructed.

FILE FORMAT (AHF documentation pp. 170-171, `external/AHF.pdf`):
  (1)  r          right edge of radial bin, **kpc/h**; a NEGATIVE r means the
                  bin is not converged (Power et al. 2003) -- dropped here.
  (12) b, (13) c  b/a and c/a of the moment-of-inertia tensor
  (14-16) Ea      largest axis   -> observer direction
  (17-19) Eb      intermediate
  (20-22) Ec      third/smallest -> disk normal
All profile quantities are CUMULATIVE "inside sphere of radius r", so the row
at r ~= 10 kpc/h is exactly the tensor Ref. [1]'s App. A describes.

🚨 THE IN-PLANE AXIS IS NEARLY DEGENERATE. At 10 kpc/h these halos have
b/a = 0.977-0.994, so `a` and `b` are within ~1% and the direction called
"major" in the plane is poorly conditioned -- it swings tens of degrees
between adjacent radial bins. `bin_sensitivity()` quantifies that, and it
bounds how well ANY reconstruction of this frame can do. The disk normal Ec
is by contrast stable, because c/a ~= 0.7 is far from 1.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ROOT / "external" / "profiles"
TARGET_KPCH = 10.0          # their App. A radius, in h^-1 kpc
#: 🚨 THE PROFILE FILES ARE SWAPPED WITHIN EACH PAIR relative to the particle
#: files. `profile_G12_AHF.txt` is the halo in
#: `HESTIA_09_18_8192_G11_minimal_particles.arrow`, and so on for all three
#: pairs. Established 2026-09-15 by a 6x6 cross-match on the cumulative mass
#: profile M(<r) alone -- no axis convention involved: under this assignment
#: every pairing agrees to RMS |log10 M_part/M_AHF| = 0.001 over
#: r = 5-40 kpc/h, against 0.058-0.098 for the naive same-name assignment.
#: The diagonal is unambiguous (see analysis/docs/muru_ahf_frames_2026-09-15.md).
#:
#: 🚨 CORRECTED 2026-09-16: it is the PARTICLE files that are out of step with
#: the paper, not the profiles. The profile names match Table I of Ref. [1]
#: by M200c AND R200c (G11: 1.944 vs 1.94e12 Msun, 178.3 vs 178.0 kpc/h; all
#: six agree), and the azimuth scan reproduces BOTH published panels to
#: 0.01-0.02 only under the swap. So particle file G11 is the paper's G1.2,
#: etc. The mapping below (profile -> particle-file key) is still right; what
#: is wrong is the paper label attached to each particle-file key. That is
#: handled in figure_render.PAPER_LABEL, not here. The "G1.1 reproduces to
#: 0.02" that misled us was a DM-only coincidence; the stellar panel was
#: 0.10 off in the same view.
NAME = {"G12": "G1.1", "G11": "G1.2", "G22": "G2.1",
        "G21": "G2.2", "G32": "G3.1", "G31": "G3.2"}


def _rows(path):
    """Converged rows only (r > 0), as (r_kpch, b_a, c_a, Ea, Eb, Ec)."""
    d = np.loadtxt(path, comments="#")
    d = d[d[:, 0] > 0]
    return [(row[0], row[11], row[12], row[13:16], row[16:19], row[19:22])
            for row in d]


def angle_deg(u, v):
    """Angle between two axes, sign-insensitive (eigenvector sign is arbitrary)."""
    c = abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(min(1.0, c))))


def frame_at(path, target=TARGET_KPCH):
    """(v_normal, v_major, info) at the converged bin nearest `target` kpc/h."""
    rows = _rows(path)
    r = np.array([x[0] for x in rows])
    i = int(np.argmin(np.abs(r - target)))
    r_i, b_a, c_a, Ea, Eb, Ec = rows[i]
    v_n = Ec / np.linalg.norm(Ec)
    v_maj = Ea - (Ea @ v_n) * v_n                      # exactly in-plane
    v_maj /= np.linalg.norm(v_maj)
    info = {"r_kpch": float(r_i), "r_kpc_physical": float(r_i / 0.677),
            "b_over_a": float(b_a), "c_over_a": float(c_a),
            "Ea_raw": [float(v) for v in Ea],
            "tilt_Ea_off_plane_deg": angle_deg(Ea, v_maj)}
    return v_n, v_maj, info


def bin_sensitivity(path, target=TARGET_KPCH):
    """How far the axes move to the adjacent converged bins — the conditioning
    of the frame, and a floor on any reconstruction's accuracy."""
    rows = _rows(path)
    r = np.array([x[0] for x in rows])
    i = int(np.argmin(np.abs(r - target)))
    out = {}
    for j, tag in ((i - 1, "prev"), (i + 1, "next")):
        if 0 <= j < len(rows):
            out[tag] = {
                "r_kpch": float(rows[j][0]),
                "d_major_deg": angle_deg(rows[i][3], rows[j][3]),
                "d_normal_deg": angle_deg(rows[i][5], rows[j][5])}
    return out


def all_frames(target=TARGET_KPCH):
    out = {}
    for f in sorted(PROFILES.glob("profile_G*_AHF.txt")):
        key = NAME[f.name.split("_")[1]]
        v_n, v_maj, info = frame_at(f, target)
        out[key] = {"v_normal": [float(v) for v in v_n],
                    "v_major": [float(v) for v in v_maj],
                    "source": f.name, **info,
                    "bin_sensitivity_deg": bin_sensitivity(f, target)}
    return out


def main():
    fr = all_frames()
    meta = json.load(open("/private/nfs/Data2/hestia_muru/figure_maps_meta.json"))
    ours = json.load(open("/private/nfs/Data2/hestia_muru/muru_frames.json"))
    print(f"AHF frames at the converged bin nearest {TARGET_KPCH} kpc/h\n")
    print(f"{'halo':6s} {'r[kpc/h]':>9s} {'b/a':>6s} {'c/a':>6s} "
          f"{'d(major,ours)':>14s} {'d(normal,ours)':>15s} "
          f"{'bin swing: maj':>15s} {'norm':>6s}")
    for h, d in fr.items():
        dm = angle_deg(np.array(d["v_major"]), np.array(meta[h]["major_axis"]))
        dn = angle_deg(np.array(d["v_normal"]), np.array(ours[h]["z_hat"]))
        sw = d["bin_sensitivity_deg"]
        smaj = max((v["d_major_deg"] for v in sw.values()), default=float("nan"))
        snrm = max((v["d_normal_deg"] for v in sw.values()), default=float("nan"))
        print(f"{h:6s} {d['r_kpch']:9.3f} {d['b_over_a']:6.3f} {d['c_over_a']:6.3f} "
              f"{dm:14.1f} {dn:15.1f} {smaj:15.1f} {snrm:6.1f}")
    out = ROOT / "analysis" / "figures" / "ahf_frames.json"
    out.write_text(json.dumps(fr, indent=1))
    print("\nwrote", out)


if __name__ == "__main__":
    main()
