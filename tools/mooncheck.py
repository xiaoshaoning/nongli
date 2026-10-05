#!/usr/bin/env python3
"""把月球位置 (λ, β, Δ) 与 ERFA/astropy 参考模型逐点比对。

参考是 `erfa.moon98`（= Meeus 第 47 章，**省略了月球平黄经的光行时改正**）。
因此本实现与它之间应有一条**系统性**黄经偏差——本脚本把它量出来，
而不是把它当成误差。

口径说明：参考的 GCRS 向量经 astropy 转到 **GeocentricMeanEcliptic**
（平黄道平分点），与 `moon_geocentric` 的约定一致；两侧同口径，参考系模型差异抵消。

用法:  python tools/mooncheck.py
"""
import os
import subprocess
import sys

import numpy as np
import astropy.units as u
import erfa
from astropy.coordinates import GCRS, GeocentricMeanEcliptic, SkyCoord
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)

AU_KM = 149597870.7
ARCSEC = 3600.0


def ref(jd):
    """erfa.moon98 的 GCRS 向量 → 平黄道平分点的 (λ, β, Δ)。"""
    t = Time(jd, format='jd', scale='tt')
    r = erfa.moon98(jd, 0.0)[0]
    sc = SkyCoord(x=r[0] * u.au, y=r[1] * u.au, z=r[2] * u.au,
                  frame=GCRS(obstime=t), representation_type='cartesian')
    e = sc.transform_to(GeocentricMeanEcliptic(equinox=t))
    return float(e.lon.deg), float(e.lat.deg), float(np.linalg.norm(r)) * AU_KM


def wrap180(x):
    return (x + 180.0) % 360.0 - 180.0


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'moon_check', '--', '-1975', '6025', '4000'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    mine, theirs, years = [], [], []
    for line in out.splitlines():
        p = line.split()
        if len(p) != 4:
            continue
        jd = float(p[0])
        mine.append((float(p[1]), float(p[2]), float(p[3])))
        theirs.append(ref(jd))
        years.append(2000.0 + (jd - 2451545.0) / 365.25)
    mine = np.array(mine)
    theirs = np.array(theirs)
    years = np.array(years)

    dlon = wrap180(mine[:, 0] - theirs[:, 0]) * ARCSEC
    dlat = (mine[:, 1] - theirs[:, 1]) * ARCSEC
    ddist = mine[:, 2] - theirs[:, 2]

    def rng(name, d, unit, lo, hi):
        m = (years >= lo) & (years < hi)
        if not m.any():
            return
        print(f'    {name} {lo:>6}..{hi:<6} n={m.sum():5}  均值={d[m].mean():+9.4f}{unit}  '
              f'rms={np.sqrt((d[m] ** 2).mean()):8.4f}{unit}  最大={np.abs(d[m]).max():9.4f}{unit}')

    print(f'样本 {len(mine)} 点，2026±4000 年')
    print('\n黄经差 (本实现 − erfa.moon98)：')
    for lo, hi in [(-1975, -1000), (-1000, 500), (500, 1500), (1500, 2200), (2200, 6025)]:
        rng('λ', dlon, '"', lo, hi)
    print(f'    {"全部":>6}            n={len(dlon):5}  均值={dlon.mean():+9.4f}"  '
          f'rms={np.sqrt((dlon ** 2).mean()):8.4f}"  最大={np.abs(dlon).max():9.4f}"')

    print('\n黄纬差：')
    for lo, hi in [(-1975, -1000), (-1000, 500), (500, 1500), (1500, 2200), (2200, 6025)]:
        rng('β', dlat, '"', lo, hi)
    print(f'    {"全部":>6}            n={len(dlat):5}  均值={dlat.mean():+9.4f}"  '
          f'rms={np.sqrt((dlat ** 2).mean()):8.4f}"  最大={np.abs(dlat).max():9.4f}"')

    print('\n距离差：')
    for lo, hi in [(-1975, -1000), (-1000, 500), (500, 1500), (1500, 2200), (2200, 6025)]:
        rng('Δ', ddist, ' km', lo, hi)
    print(f'    {"全部":>6}            n={len(ddist):5}  均值={ddist.mean():+9.4f} km  '
          f'rms={np.sqrt((ddist ** 2).mean()):8.4f} km  最大={np.abs(ddist).max():9.4f} km')

    print(f'\n黄经系统差 = {dlon.mean():+.4f}"  '
          f'(月球地心角速度 0.549"/s → 等效 {abs(dlon.mean())/0.549:.2f} s 的光行时)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
