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

`measure()` 只算不印，`main()` 只印不算——`tools/gen_accuracy.py` 复用前者
生成 `docs/accuracy.md`（精度数字的唯一来源）。

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
C_KM_S = 299792.458
ERAS = [(1900, 1925), (1925, 1950), (1950, 1975), (1975, 2000), (2000, 2025), (2025, 2050)]


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def measure(kernel=KERNEL):
    """跑一遍全部比对。内核不存在时返回 None。"""
    if not os.path.exists(kernel):
        return None

    import erfa
    from jplephem.spk import SPK
    K = SPK.open(kernel)

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

    # === 太阳 ===
    # 主指标取**几何**黄经：DE421 的地心矢量经 erfa.ecm06 转到 IAU2006 当日平黄道，
    # 不含光行差与章动——两边同口径，验收不依赖我们自己的任何模型。
    jd = sun_rows[:, 0]
    d_geom = np.empty(len(jd))
    d_app = np.empty(len(jd))
    d_sdist = np.empty(len(jd))
    for i, j in enumerate(jd):
        g = sun(j) - earth(j)
        r = np.linalg.norm(g)
        w = erfa.ecm06(j, 0.0) @ (g / r)
        lam_geom = np.degrees(np.arctan2(w[1], w[0]))
        # 视黄经要在真值侧补上章动与周年光行差——这就复刻了 sun.rs 的口径，
        # 所以只作辅助指标（万一我们两边犯同样的错，它看不出来）。
        lam_app = (lam_geom + erfa.nut06a(j, 0.0)[0] / np.pi * 180.0
                   - 20.4898 / (r / AU_KM) / 3600.0)
        d_geom[i] = wrap(sun_rows[i, 2] - lam_geom) * 3600.0
        d_app[i] = wrap(sun_rows[i, 1] - lam_app) * 3600.0
        d_sdist[i] = (sun_rows[i, 3] - r / AU_KM) * AU_KM

    # === 月球 ===
    # 我们的 moon_geocentric 给的是**视**位置：先在 t−τ 处取几何位置（光行时），
    # 再加周年光行差——与 crates/ephemeris/src/moon.rs 逐条对应。
    jdm = moon_rows[:, 0]
    d_lon = np.empty(len(jdm))
    d_lat = np.empty(len(jdm))
    d_mdist = np.empty(len(jdm))
    for i, j in enumerate(jdm):
        r = np.linalg.norm(moon(j) - earth(j))
        js = j - r / C_KM_S / 86400.0
        g = moon(js) - earth(js)
        w = erfa.ecm06(j, 0.0) @ (g / np.linalg.norm(g))
        lam = np.degrees(np.arctan2(w[1], w[0])) % 360.0
        # 周年光行差只与**地球绕日**的速度有关，所以用日地距离而不是月地距离
        # （这里曾写错过：用月地距离 0.0026 AU 代进去会得到 7900"）。
        r_sun = np.linalg.norm(sun(j) - earth(j))
        lam -= 20.4898 / (r_sun / AU_KM) / 3600.0
        d_lon[i] = wrap(moon_rows[i, 1] - lam) * 3600.0
        d_lat[i] = (moon_rows[i, 2] - np.degrees(np.arcsin(w[2]))) * 3600.0
        d_mdist[i] = moon_rows[i, 3] - r

    # === 顺带量一下"别人当月球真值"的那些东西到底有多准 ===
    # 这一行以前是手抄在文档里的（19.8″ / 43.7″），抄错了不会有人知道。
    builtin = _builtin_moon_error(K, earth, sun)

    def stat(d):
        return [float(np.median(d)), float(np.median(abs(d))), float(abs(d).max())]

    y = 1900.0 + (sun_rows[:, 0] - 2415020.5) / 365.25
    eras = []
    for lo, hi in ERAS:
        m = (y >= lo) & (y < hi)
        eras.append((lo, hi, float(np.median(d_geom[m])),
                     float(np.median(d_lon[m])), float(np.median(d_mdist[m]))))
    return {
        'n': len(jd),
        'kernel': os.path.basename(kernel),
        'kernel_span': (float(K.segments[0].start_jd), float(K.segments[0].end_jd)),
        'sun_geom': stat(d_geom),
        'sun_app': stat(d_app),
        'sun_dist': stat(d_sdist),
        'moon_lon': stat(d_lon),
        'moon_lat': stat(d_lat),
        'moon_dist': stat(d_mdist),
        'builtin_moon': builtin,
        'eras': eras,
    }


def _builtin_moon_error(K, earth, sun, n=400):
    """astropy 的 builtin 月球（ERFA moon98 一路）对 DE421 的黄经误差，角秒。"""
    try:
        from astropy.coordinates import (
            GeocentricTrueEcliptic, get_body, solar_system_ephemeris,
        )
        from astropy.time import Time
    except ImportError:
        return None
    solar_system_ephemeris.set('builtin')
    jds = np.linspace(K.segments[0].start_jd + 100, K.segments[0].end_jd - 100, n)
    errs = []
    for j in jds:
        t = Time(j, format='jd', scale='tt')
        a = get_body('moon', t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg
        g = K[(0, 3)].compute(j) + K[(3, 301)].compute(j) - earth(j)
        r = np.linalg.norm(g)
        w = __import__('erfa').ecm06(j, 0.0) @ (g / r)
        lam = np.degrees(np.arctan2(w[1], w[0])) % 360.0
        r_sun = np.linalg.norm(sun(j) - earth(j))
        lam -= 20.4898 / (r_sun / AU_KM) / 3600.0
        errs.append(wrap(a - lam) * 3600.0)
    e = np.array(errs)
    return [float(np.median(abs(e))), float(abs(e).max())]


def main():
    m = measure()
    if m is None:
        print(f'缺少 {KERNEL}——跳过对 JPL 真值的验收。取法见本脚本头部注释。')
        return 0

    print(f'对 JPL {m["kernel"]}（真值）· {m["n"]} 个时刻 · 1900–2050\n')
    print(f'{"":10} {"中位":>10} {"中位|差|":>11} {"最大":>10}')
    for name, key, dp, unit, note in [
        ('太阳几何黄经', 'sun_geom', 4, '"', ''),
        ('太阳视黄经', 'sun_app', 4, '"', '   <- 辅助'),
        ('日地距离', 'sun_dist', 2, ' km', ''),
        ('月球黄经', 'moon_lon', 4, '"', ''),
        ('月球黄纬', 'moon_lat', 4, '"', ''),
        ('月地距离', 'moon_dist', 3, ' km', ''),
    ]:
        med, mad, mx = m[key]
        print(f'{name:12} {med:+9.{dp}f}{unit}  {mad:10.{dp}f}{unit}  {mx:9.{dp}f}{unit}{note}')

    if m['builtin_moon']:
        mad, mx = m['builtin_moon']
        print(f'\n顺带：astropy builtin 月球（= ERFA moon98 一路）对同一真值的黄经'
              f'：中位|差| {mad:.3f}"  最大 {mx:.3f}"')
        print('      —— 所以它不能当月球的真值。')

    print('\n分年代（看有没有趋势——趋势意味着历元约定还没对上）：')
    for lo, hi, sg, ml, md in m['eras']:
        print(f'  {lo}-{hi}: 太阳几何 {sg:+8.4f}"  月球 {ml:+9.4f}"  月距 {md:+8.3f} km')

    print('''
读法：
* **主指标是"太阳几何黄经"**：两边同口径、不含光行差与章动，所以验收不依赖
  我们自己的模型。"太阳视黄经"那一行要在真值侧补光行差与章动，等于复刻了
  sun.rs 的口径，只能当辅助。
* 太阳的残差应当**没有随年代的趋势**（J2000 处过零、向两侧线性增长就是
  历元约定没对上，见 crates/ephemeris/src/frames.rs）。剩下的是 VSOP87 相对
  DE421 的理论差；月球则已换成 ELP2000-82B，残差只剩 0.0x"。
* DE421 只到 2053，所以这里看不到 ±4000 年；那边两个理论都在外推，
  真值得用 DE441（本机取不到，见 README）。
''')
    return 0


if __name__ == '__main__':
    sys.exit(main())
