#!/usr/bin/env python3
"""对 **JPL 星历内核**（真值）验收太阳与月球的位置。

其它 `*check.py` 比的是 ERFA / astropy 的**模型**——`epv00` 是 VSOP2000、
`moon98` 是截断的 ELP，它们各有各的历元约定与截断误差，会把我们的真误差
掩掉（这个项目已经栽过两次：一次把岁差约定之差当成截断误差，一次差点
把 moon98 的约定当成误差）。这里比的是真值。

需要的内核：`kernels/de421.bsp`（JPL DE421，覆盖 1900–2053）。
它不进版本库，取法：

    pip download skyfield-data -d /tmp/sfd --no-deps
    cd /tmp/sfd && unzip -o -q skyfield_data-*.whl -d x
    mkdir -p kernels && cp x/skyfield_data/data/de421.bsp kernels/

没有这个文件时本脚本会打印取法并正常退出（不算失败）。

用法:  python tools/truthcheck.py
"""
import os
import subprocess
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
KERNEL = os.path.join(ROOT, 'kernels', 'de421.bsp')
AU_KM = 1.495978707e8


def main():
    if not os.path.exists(KERNEL):
        print(f'缺少 {KERNEL}——跳过对 JPL 真值的验收。取法见本脚本头部注释。')
        return 0

    import erfa
    from jplephem.spk import SPK
    K = SPK.open(KERNEL)

    def earth(jd):
        return K[(0, 3)].compute(jd) + K[(3, 399)].compute(jd)

    def sun(jd):
        return K[(0, 10)].compute(jd)

    def moon(jd):
        return K[(0, 3)].compute(jd) + K[(3, 301)].compute(jd)

    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'truth_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    sun_rows = np.array([[float(x) for x in l.split()[1:]]
                         for l in out.splitlines() if l.startswith('SUN')])
    moon_rows = np.array([[float(x) for x in l.split()[1:]]
                          for l in out.splitlines() if l.startswith('MOON')])

    def wrap(x):
        return (x + 180.0) % 360.0 - 180.0

    # === 太阳 ===
    jd = sun_rows[:, 0]
    d_sun = np.empty(len(jd))
    d_dist = np.empty(len(jd))
    for i, j in enumerate(jd):
        g = sun(j) - earth(j)
        r = np.linalg.norm(g)
        w = erfa.ecm06(j, 0.0) @ (g / r)
        # 视黄经 = 当日平黄道 + 章动；再加周年光行差（我们的 sun_apparent_longitude 也含）
        lam = (np.degrees(np.arctan2(w[1], w[0]))
               + erfa.nut06a(j, 0.0)[0] / np.pi * 180.0
               - 20.4898 / (r / AU_KM) / 3600.0)
        d_sun[i] = wrap(sun_rows[i, 1] - lam) * 3600.0
        d_dist[i] = (sun_rows[i, 2] - r / AU_KM) * AU_KM

    # === 月球 ===
    jdm = moon_rows[:, 0]
    d_lon = np.empty(len(jdm))
    d_lat = np.empty(len(jdm))
    d_km = np.empty(len(jdm))
    for i, j in enumerate(jdm):
        g = moon(j) - earth(j)
        r = np.linalg.norm(g)
        w = erfa.ecm06(j, 0.0) @ (g / r)
        d_lon[i] = wrap(moon_rows[i, 1] - np.degrees(np.arctan2(w[1], w[0]))) * 3600.0
        d_lat[i] = (moon_rows[i, 2] - np.degrees(np.arcsin(w[2]))) * 3600.0
        d_km[i] = moon_rows[i, 3] - r

    y = 1900.0 + (sun_rows[:, 0] - 2415020.5) / 365.25
    print(f'对 JPL DE421（真值）· {len(jd)} 个时刻 · 1900–2050\n')
    print(f'{"":10} {"中位":>10} {"中位|差|":>11} {"最大":>10}')
    print(f'{"太阳黄经":10} {np.median(d_sun):+9.4f}" {np.median(abs(d_sun)):10.4f}" {abs(d_sun).max():9.4f}"')
    print(f'{"日地距离":10} {np.median(d_dist):+9.2f}  {np.median(abs(d_dist)):10.2f}  {abs(d_dist).max():9.2f} km')
    print(f'{"月球黄经":10} {np.median(d_lon):+9.4f}" {np.median(abs(d_lon)):10.4f}" {abs(d_lon).max():9.4f}"')
    print(f'{"月球黄纬":10} {np.median(d_lat):+9.4f}" {np.median(abs(d_lat)):10.4f}" {abs(d_lat).max():9.4f}"')
    print(f'{"月地距离":10} {np.median(d_km):+9.3f}  {np.median(abs(d_km)):10.3f}  {abs(d_km).max():9.3f} km')

    print('\n分年代（看有没有趋势——趋势意味着历元约定还没对上）：')
    for lo, hi in [(1900, 1925), (1925, 1950), (1950, 1975), (1975, 2000), (2000, 2025), (2025, 2050)]:
        m = (y >= lo) & (y < hi)
        print(f'  {lo}-{hi}: 太阳 {np.median(d_sun[m]):+8.4f}"  月球 {np.median(d_lon[m]):+9.4f}"'
              f'  月距 {np.median(d_km[m]):+8.3f} km')

    print('''
读法：
* 太阳的残差应当**没有随年代的趋势**（J2000 处过零、向两侧线性增长就是
  历元约定没对上，见 crates/ephemeris/src/frames.rs 的
  `dynamical_to_iau2006_lon_deg`）。剩下的是 VSOP87 相对 DE421 的理论差。
* 月球现在由 Meeus 第 47 章那 60 项截断主导（第 7 步换 ELP2000-82B 就是冲它）。
* DE421 只到 2053，所以这里看不到 ±4000 年；那边两个理论都在外推，
  真值得用 DE441。
''')
    return 0


if __name__ == '__main__':
    sys.exit(main())
