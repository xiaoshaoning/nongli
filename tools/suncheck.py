#!/usr/bin/env python3
"""太阳位置的验收：拿 `ephemeris --example sun_check` 的输出对 astropy 比。

参考是 astropy 的 `builtin` 星历，也就是 ERFA 的 `epv00`——**VSOP2000**
（Moisson & Bretagnon 2001），不是 VSOP87。

**这份表不能用来判断我们的绝对精度。** 它只适用 1900–2100，出了这段参考自己
就在外推（文档写着 1500/2500 误差×10、1000/3000 ×60）。要看绝对精度请用
`tools/truthcheck.py`（对 JPL DE421 真值）。这里留着是为了看**长跨度**上两边
模型是怎么分离的。

用法:  python tools/suncheck.py
"""
import os
import subprocess
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from astropy.coordinates import (  # noqa: E402
    GeocentricTrueEcliptic, get_body, solar_system_ephemeris,
)
from astropy.time import Time  # noqa: E402

solar_system_ephemeris.set('builtin')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
J2000 = 2451545.0
AU_KM = 1.495978707e8
ERAS = [(-4000, -3000), (-3000, -2000), (-2000, -1000), (-1000, 0),
        (-200, 200), (1000, 2000), (3000, 4000)]


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'sun_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    rows = [l.split() for l in out.splitlines() if l.startswith('SUN')]
    Y = np.array([int(r[1]) for r in rows])
    app = np.array([float(r[2]) for r in rows])
    dist = np.array([float(r[3]) for r in rows])

    jd = J2000 + (Y - 2000.0) * 365.25
    t = Time(jd, format='jd', scale='tt')
    ref = get_body('sun', t).transform_to(GeocentricTrueEcliptic(equinox=t))
    d = ((app - ref.lon.deg + 180.0) % 360.0 - 180.0) * 3600.0
    dR = (dist - ref.distance.au) * AU_KM

    print('太阳视黄经 相对 astropy builtin（= ERFA epv00，VSOP2000）')
    print(f'{"区间":>16} {"中位":>9} {"最大":>9}   换算成节气时刻')
    for lo, hi in ERAS:
        m = (Y >= lo) & (Y < hi)
        if not m.any():
            continue
        med, mx = np.median(d[m]), np.abs(d[m]).max()
        print(f'  {lo:>5}..{hi:<6} {med:+8.3f}" {mx:8.3f}"   {med*24.5:+7.1f}s / {mx*24.5:6.1f}s')

    print(f'\n日地距离：最大差 {np.abs(dR).max():.2f} km（1 AU = {AU_KM/1e6:.4f}e6 km）')
    for lo, hi in [(-4000, -2000), (-2000, -1000), (-1000, 1000), (1000, 3000), (3000, 4000)]:
        m = (Y >= lo) & (Y < hi)
        if m.any():
            print(f'  {lo:>5}..{hi:<6} 最大 {np.abs(dR[m]).max():9.2f} km')

    print('''
读法：**先看有没有随年代的趋势，别急着当误差**。

我们这一侧的历元约定改正已经做了（`frames.rs::dynamical_to_iau2006_lon_offset_deg`，
用 VSOP87B + IAU2006 定出，并与 IAU1976−IAU2006 的岁差速率差交叉验证过）。
所以现代区间只剩 VSOP87 相对 VSOP2000 的理论差（~0.01"）。

剩下的量级随 |年−2000| 增长的那些，主要来自**参考自己**：`epv00` 只适用
1900–2100，出了这段它自己就在外推。要看绝对精度请用 `tools/truthcheck.py`
（对 JPL DE421 真值，1900–2050）。这里留着是为了看**长跨度**上两边模型怎么
分离，以及提醒：离 J2000 越远，任何一方的外推都越不可信。

为什么"朔"一直没事：它用日月**黄经之差**，两套理论各带各的岁差，作差时抵消
（`daydiff` 里朔的中位差只有 1.3 s）。
''')
    return 0


if __name__ == '__main__':
    sys.exit(main())
