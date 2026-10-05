#!/usr/bin/env python3
"""太阳位置的验收：拿 `ephemeris --example sun_check` 的输出对 astropy 比。

参考是 astropy 的 `builtin` 星历，也就是 ERFA 的 `epv00`——**VSOP2000**
（Moisson & Bretagnon 2001），不是 VSOP87。这一点很要紧，因为两者在离
J2000 远处会分道扬镳（见下）。

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
    dR = (dist - ref.distance.au) * 1.495978707e8   # km

    print('太阳视黄经 相对 astropy builtin（= ERFA epv00，VSOP2000）')
    print(f'{"区间":>16} {"中位":>9} {"最大":>9}   换算成节气时刻')
    for lo, hi in [(-4000, -3000), (-3000, -2000), (-2000, -1000), (-1000, 0),
                   (-200, 200), (1000, 2000), (3000, 4000)]:
        m = (Y >= lo) & (Y < hi)
        if not m.any():
            continue
        med, mx = np.median(d[m]), np.abs(d[m]).max()
        print(f'  {lo:>5}..{hi:<6} {med:+8.3f}" {mx:8.3f}"   {med*24.5:+7.1f}s / {mx*24.5:6.1f}s')

    print(f'\n日地距离：最大差 {np.abs(dR).max():.2f} km（1 AU = 1.496e8 km）')
    for lo, hi in [(-4000, -2000), (-2000, -1000), (-1000, 1000), (1000, 3000), (3000, 4000)]:
        m = (Y >= lo) & (Y < hi)
        if m.any():
            print(f'  {lo:>5}..{hi:<6} 最大 {np.abs(dR[m]).max():9.2f} km')

    print('''
读法：这个差**主要是历元/岁差约定的差别，不是我们的截断误差**。
它随 |年−2000| 平滑地线性增长，在 J2000 处正好过零——和 `gst06a` 那次
是同一类情形：换了模型，两者在 J2000 附近一致，走远才分家。

原因：VSOP87D 把黄经直接给在它自己的"当日动力学平黄道分点"里，岁差是
1988 年那套；astropy 把 epv00（VSOP2000）的 ICRS 坐标用 **IAU 2006/P03**
岁差转到当日黄道。两套岁差的速率差 ~0.003″/年，累积成上表的斜率。

**截断误差是另一回事，已经单独量过**：我们的表相对未截断的 VSOP87D
在 1526–2526 最大 0.0025″（见 tools/gen_vsop87.py），比上面的差小三个量级。

为什么"朔"没事：朔用的是日月**黄经之差**，两套理论各带各的岁差，作差时
抵消了（daydiff 里朔的中位差只有 1.3 s）。

两点自诚：
1. **参考到远处也不可信。** epv00 自己写着只适用 1900–2100，且 1500/2500
   误差×10、1000/3000 ×60。上表 ±3000 以外的数里既有它的外推误差，
   也有 VSOP87 的。真要那个区间，得换 JPL DE440 当参考。
2. **距离差证明"纯转动"的说法只能算近似。** 转动不改 R，而 R 在极端历元
   差到上万 km，说明那边两个理论确实各走各的。现代区间（±1000 内）距离
   只差几 km。
''')
    return 0


if __name__ == '__main__':
    sys.exit(main())
