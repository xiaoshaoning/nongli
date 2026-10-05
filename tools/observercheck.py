#!/usr/bin/env python3
"""把恒星时、观测者几何、站心视差对 ERFA 比对。

**为什么不用 erfa.atco13 作地平坐标的参考**（两条独立的理由）：

1. 它是"恒星例程"：`px` 是**年视差**，模型假定天体距离远大于地球半径，
   因此**无法表达月球的周日视差**。实测把 px 设成月球的年视差 (0.5366") 只让
   结果动 0.00006"，而月球周日视差实际是 3422"。
2. 它自动加**年周光行差** (~20.5")；而本实现的 `horizontal()` 是给月球用的——
   月球的光行时已经在理论里，年周光行差**故意不加**。

所以按步骤分开验：恒星时对 `gmst06`/`gst06a`，观测者几何对 `gd2gc`，
站心视差用 `gd2gc` 的矢量在 Python 里独立减一遍。

用法:  python tools/observercheck.py
"""
import os
import subprocess
import sys

import numpy as np
import erfa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
sys.path.insert(0, os.path.join(ROOT, 'probe'))
import dtreference  # noqa: E402

AS = 206264.806
D2R = np.pi / 180.0
PRESSURE_HPA = 1010.0
TEMPERATURE_C = 10.0


def ut1_jd(tt_jd):
    """我们的 UT1：TT − ΔT。"""
    return tt_jd - float(dtreference.delta_t_from_jd(tt_jd)) / 86400.0


def wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'observer_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    gst, hz, obs = [], [], []
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'GST':
            gst.append(tuple(float(x) for x in p[1:]))
        elif p[0] == 'OBS':
            obs.append(tuple(float(x) for x in p[1:]))
        elif p[0] == 'HZ':
            hz.append([float(x) for x in p[1:]])

    # === 恒星时 ===
    print('=== 恒星时（对 erfa.gmst06 / erfa.gst06a）===')
    Y = np.array([g[0] for g in gst])
    jd = 2451545.0 + (Y - 2000.0) * 365.25
    ut1 = np.array([ut1_jd(x) for x in jd])
    ref_m = np.array([erfa.gmst06(ut1[i], 0.0, jd[i], 0.0) for i in range(len(Y))])
    ref_a = np.array([erfa.gst06a(ut1[i], 0.0, jd[i], 0.0) for i in range(len(Y))])
    mine_m = np.array([g[1] for g in gst])
    mine_a = np.array([g[2] for g in gst])
    em = np.abs(wrap(mine_m - ref_m / D2R)) * 3600
    ea = np.abs(wrap(mine_a - ref_a / D2R)) * 3600

    # 分点差的**定义**（erfa 的 ee00.c）：ee = dpsi*cos(epsA) + eect00。
    # 前一项用**平**黄赤交角。拿 erfa 自己的 nut06a/obl06 拼出来，各历元都该对上。
    dpsi = np.array([erfa.nut06a(jd[i], 0.0)[0] for i in range(len(Y))])
    epsa = np.array([erfa.obl06(jd[i], 0.0) for i in range(len(Y))])
    eect = np.array([erfa.eect00(jd[i], 0.0) for i in range(len(Y))])
    ref_def = ref_m + dpsi * np.cos(epsa)          # 不含补充项：本实现的口径
    ref_iau = ref_def + eect                       # 含补充项：完整的 IAU2000 分点差
    ed = np.abs(wrap(mine_a - ref_def / D2R)) * 3600
    ei = np.abs(wrap(mine_a - ref_iau / D2R)) * 3600

    for lo, hi in [(-4000, -1000), (-1000, -100), (-100, 100), (100, 1000), (1000, 4000)]:
        s = (Y >= lo) & (Y < hi)
        if s.any():
            print(f'  {lo:>6}..{hi:<6} n={s.sum():4}  GMST 最大 {em[s].max():10.5f}"  '
                  f'分点差定义最大 {ed[s].max():10.5f}"')
    print(f'  GMST 全程最大 {em.max():.5f}"，GAST 对"定义"全程最大 {ed.max():.5f}"')
    print(f'  差 IAU2000 完整分点差 {ei.max():.5f}"——这就是未含的补充项 eect00（±0.003"）。')

    # gst06a 走的是 CIO 路线（ERA - 分点原点差 eors），而 eors 来自 s06 级数。
    # erfa 自己的 s06.c 注明："s remains below 0.1 arcsecond throughout 1900-2100"
    # ——超出这段，是它的多项式外推在发散，不是我们的 GAST 错。
    ok = (Y >= 1900) & (Y <= 2100)
    print('')
    print(f'  对 erfa.gst06a（CIO 路线）：1900-2100 内最大 {ea[ok].max():.5f}"；')
    print(f'  全程最大 {ea.max():.2f}"，但它出现在 |年-2000| 很大处——'
          f'那段 erfa 的 s06 已超出自己声明的适用范围。')
    print(f'  旁证：s06 是 5 次级数，故 gst06a-gmst06 与 Δψcos(ε_A) 之差按 t^3 发散'
          f'（t=4 世纪 1.4"，t=40 世纪 2327"），而真正的分点差是有界的 ±20"。')

    # === 观测者几何 ===
    print('\n=== 观测者地心矢量（地球固连，对 erfa.gd2gc / WGS84）===')
    seen, worst = set(), 0.0
    for lat, lon, h, x, y, z in obs:
        if (lat, lon, h) in seen:
            continue
        seen.add((lat, lon, h))
        g = erfa.gd2gc(1, lon * D2R, lat * D2R, h)          # 返回米
        worst = max(worst, float(np.linalg.norm(np.array([x, y, z]) * 1000.0 - np.array(g))))
    print(f'  {len(seen)} 个站点（含赤道、南半球、高纬、负海拔），最大偏差 {worst:.6f} m')

    # === 站心视差 ===
    print('\n=== 月球站心视差（独立用 gd2gc 的矢量减一遍）===')
    moon = []
    for y, lat, lon, h, lam, bet, az_v, alt_v, az_r, alt_r, ra, dec, dist in hz:
        if dist > 1e9:
            continue
        tt = 2451545.0 + (y - 2000.0) * 365.25
        ut1 = ut1_jd(tt)
        v_tod = np.array([np.cos(dec * D2R) * np.cos(ra * D2R),
                          np.cos(dec * D2R) * np.sin(ra * D2R),
                          np.sin(dec * D2R)])
        t_gast = np.degrees(erfa.gst06a(ut1, 0.0, tt, 0.0))
        xyz = np.array(erfa.gd2gc(1, lon * D2R, lat * D2R, h)) / 1000.0   # km
        th = t_gast * D2R
        v = dist * v_tod - np.array([xyz[0] * np.cos(th) - xyz[1] * np.sin(th),
                                     xyz[0] * np.sin(th) + xyz[1] * np.cos(th),
                                     xyz[2]])
        r = np.linalg.norm(v)
        dd, aa = np.arcsin(v[2] / r), np.arctan2(v[1], v[0])
        H = (t_gast + lon) * D2R - aa
        phi = lat * D2R
        alt = np.arcsin(np.sin(phi) * np.sin(dd) + np.cos(phi) * np.cos(dd) * np.cos(H))
        az = np.degrees(np.arctan2(-np.cos(dd) * np.sin(H),
                                   np.sin(dd) * np.cos(phi)
                                   - np.cos(dd) * np.sin(phi) * np.cos(H))) % 360.0
        if alt_v > 0.0:
            moon.append((wrap(az_v - az) * 3600, (alt_v - np.degrees(alt)) * 3600))
    moon = np.array(moon)
    print(f'  样本 {len(moon)}（只取地平以上）')
    print(f'  Δaz  最大 {np.abs(moon[:, 0]).max():.4f}"')
    print(f'  Δalt 最大 {np.abs(moon[:, 1]).max():.4f}"')
    print('  参考未含周日光行差（~0.32"），我们的含——它就在 Δalt 里。')
    return 0


if __name__ == '__main__':
    sys.exit(main())
