#!/usr/bin/env python3
"""Cross-check the Rust ephemeris against ERFA / astropy's `builtin` ephemeris.

Reference: astropy.coordinates solar_system_ephemeris='builtin'
  * Sun  -> erfa.epv00 (VSOP87 MHB2000 truncation)
  * Moon -> erfa.moon98 (ELP2000-82B truncation, Meeus ch.47)

Usage:  python tools/crosscheck.py
"""
import math
import os
import subprocess
import sys

import numpy as np
from astropy.coordinates import GeocentricTrueEcliptic, get_body, solar_system_ephemeris
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)

solar_system_ephemeris.set('builtin')


def arcsec(a, b):
    d = (a - b + 180.0) % 360.0 - 180.0
    return d * 3600.0


def ref_longitudes(jd_tt):
    t = Time(jd_tt, format='jd', scale='tt')
    sun = get_body('sun', t)
    moon = get_body('moon', t)
    return (
        sun.transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg,
        moon.transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg,
    )


def main():
    y0, y1 = -1975, 6025
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '--example', 'ephemeris_dump', '--',
         str(y0), str(y1)],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    sun_d, moon_d, nm_d, st_d = [], [], [], []
    sun_y, moon_y = [], []
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'EPH':
            jd = float(p[1])
            a, b = ref_longitudes(jd)
            sun_d.append(arcsec(float(p[2]), a))
            moon_d.append(arcsec(float(p[3]), b))
            sun_y.append(2000.0 + (jd - 2451545.0) / 365.25)
            moon_y.append(sun_y[-1])
        elif p[0] == 'NM':
            k = int(p[1])
            jd = float(p[2])
            nm_d.append((k, jd))
        elif p[0] == 'ST':
            st_d.append((int(p[1]), float(p[2])))

    def report(name, errs, years):
        e = np.array(errs)
        y = np.array(years)
        print(f'{name}: n={len(e)}  rms={np.sqrt((e**2).mean()):.3f}"  '
              f'max={np.abs(e).max():.3f}"')
        for lo, hi in [(-1975, -1000), (-1000, 1500), (1500, 2200), (2200, 6025)]:
            m = (y >= lo) & (y < hi)
            if m.any():
                print(f'    {lo:6}..{hi:<6} n={m.sum():5}  rms={np.sqrt((e[m]**2).mean()):7.3f}"  '
                      f'max={np.abs(e[m]).max():8.3f}"')

    report('太阳视黄经 (Rust vs ERFA)', sun_d, sun_y)
    report('月球视黄经 (Rust vs ERFA)', moon_d, moon_y)

    # 朔/节气的参考时刻：由 astropy 的经度自行求根，再与 Rust 比。
    def solve(jd, target, moon=False):
        """Newton on the apparent longitude difference (deg)."""
        for _ in range(12):
            if moon:
                a, b = ref_longitudes(jd)
                y = (a - b + 180.0) % 360.0 - 180.0
            else:
                a, _ = ref_longitudes(jd)
                y = (a - target + 180.0) % 360.0 - 180.0
            h = 0.02
            if moon:
                a1, b1 = ref_longitudes(jd + h)
                a2, b2 = ref_longitudes(jd - h)
                y1 = (a1 - b1 + 180.0) % 360.0 - 180.0
                y2 = (a2 - b2 + 180.0) % 360.0 - 180.0
            else:
                a1, _ = ref_longitudes(jd + h)
                a2, _ = ref_longitudes(jd - h)
                y1 = (a1 - target + 180.0) % 360.0 - 180.0
                y2 = (a2 - target + 180.0) % 360.0 - 180.0
            dy = (y1 - y2) / (2 * h)
            jd -= y / dy
        return jd

    nm_err = [abs(solve(jd, 0.0, moon=True) - jd) * 86400.0 for _, jd in nm_d]
    st_err = []
    for j, jd in st_d:
        target = (15.0 * j) % 360.0
        st_err.append(abs(solve(jd, target) - jd) * 86400.0)
    print(f'\n朔   时刻差: n={len(nm_err)}  中位={np.median(nm_err):.2f}s  '
          f'最大={max(nm_err):.2f}s')
    print(f'节气 时刻差: n={len(st_err)}  中位={np.median(st_err):.2f}s  '
          f'最大={max(st_err):.2f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
