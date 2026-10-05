#!/usr/bin/env python3
"""Count calendar-day disagreements between the Rust implementation and the
ERFA/astropy reference for 朔 (month boundaries) and 节气.

This is the number that actually matters: a timing error of a few seconds only
matters when a 朔 or 节气 falls within seconds of Beijing midnight.

The reference instant is obtained from the Rust instant by one Newton step
using the reference longitudes (vectorised), which is accurate to <<1 s.

Usage:  python tools/daydiff.py [y0] [y1]
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
sys.path.insert(0, os.path.join(ROOT, 'probe'))
import dtreference  # noqa: E402  (ΔT 参考实现，与 crates/ephemeris 保持一致)

solar_system_ephemeris.set('builtin')


# **注意（2026-10-05 起）：本工具里的月球参考现在已经比本实现差了。**
# 上面用的 `erfa.moon98` / astropy builtin 对 JPL DE421 的月球黄经误差是
# **中位 ~20″、最大 ~44″**；本实现是 **中位 0.044″、最大 0.11″**。所以月球相关的
# 残差现在量的是**参考的误差**，不是我们的。要量我们的绝对精度请用
# `tools/truthcheck.py`（对 JPL DE421 真值）。


def beijing_day(jd_tt):
    jd_tt = np.asarray(jd_tt, dtype=float)
    jd_ut = jd_tt - dtreference.delta_t_from_jd(jd_tt) / 86400.0
    return np.floor(jd_ut + 8.0 / 24.0 + 0.5)


def wrap180(x):
    return (x + 180.0) % 360.0 - 180.0


def lon(jd, body):
    t = Time(jd, format='jd', scale='tt')
    return get_body(body, t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg


def refine(jd, target, moon):
    """One Newton step on the reference longitudes."""
    h = 0.02
    if moon:
        y = wrap180(lon(jd, 'moon') - lon(jd, 'sun'))
        y1 = wrap180(lon(jd + h, 'moon') - lon(jd + h, 'sun'))
        y2 = wrap180(lon(jd - h, 'moon') - lon(jd - h, 'sun'))
    else:
        y = wrap180(lon(jd, 'sun') - target)
        y1 = wrap180(lon(jd + h, 'sun') - target)
        y2 = wrap180(lon(jd - h, 'sun') - target)
    return jd - y / ((y1 - y2) / (2 * h))


def main():
    y0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1526
    y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 2526
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'nongli', '--example', 'ephemeris_dump', '--',
         str(y0), str(y1), '1', '1'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    nm, st = [], []
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'NM':
            nm.append(float(p[2]))
        elif p[0] == 'ST':
            st.append((int(p[1]), float(p[2])))

    mj = np.array(nm)
    sj = np.array([jd for _, jd in st])
    starget = np.array([(15.0 * j) % 360.0 for j, _ in st])

    dt_nm = (refine(mj, 0.0, True) - mj) * 86400.0
    dt_st = (refine(sj, starget, False) - sj) * 86400.0

    bad_nm = np.nonzero(beijing_day(mj) != beijing_day(refine(mj, 0.0, True)))[0]
    bad_st = np.nonzero(beijing_day(sj) != beijing_day(refine(sj, starget, False)))[0]

    print(f'范围 {y0}..{y1}')
    print(f'  朔   n={len(nm):6}  |Δt| 中位={np.median(abs(dt_nm)):6.2f}s '
          f'最大={abs(dt_nm).max():7.2f}s   日界不符 {len(bad_nm)}  索引 {bad_nm[:8].tolist()}')
    print(f'  节气 n={len(st):6}  |Δt| 中位={np.median(abs(dt_st)):6.2f}s '
          f'最大={abs(dt_st).max():7.2f}s   日界不符 {len(bad_st)}  索引 {bad_st[:8].tolist()}')
    # 日界不符是"事件恰好贴近午夜"的巧合。margin 均匀分布于 [0, 43200] 秒，
    # 故单个事件的翻日概率 = |Δt|/43200。
    print(f'  (期望日界不符：朔 {abs(dt_nm).sum()/43200:.1f} 人次, '
          f'节气 {abs(dt_st).sum()/43200:.1f} 人次)')
    for i in bad_st[:5]:
        print(f'    节气 j={st[i][0]} rust={mj[0]*0:.0f} t={sj[i]:.5f} Δ={dt_st[i]:+.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
