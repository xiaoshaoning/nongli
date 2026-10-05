#!/usr/bin/env python3
"""把恒星时与地平坐标对 ERFA / astropy 比对。

**参考的口径**：astropy 的 `AltAz`（经 CIRS → 观测坐标），它用的是
ERFA 的 `refco` 折射模型。本实现只在**标准大气**下匹配它（见下）。

用法:  python tools/observercheck.py
"""
import os
import subprocess
import sys

import numpy as np
import erfa
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, GeocentricTrueEcliptic, SkyCoord
from astropy.time import Time
from astropy.utils import iers

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
sys.path.insert(0, os.path.join(ROOT, 'probe'))
import dtreference  # noqa: E402

iers.conf.auto_download = False
iers.conf.auto_max_age = None


def ut1_jd(tt_jd):
    """我们的 UT1：TT − ΔT。"""
    return tt_jd - float(dtreference.delta_t_from_jd(tt_jd)) / 86400.0


def mk_time(tt_jd):
    """TT 儒略日 → astropy Time，并**喂进我们自己的 ΔT**。

    否则 astropy 会用它的 IERS 表（远期无数据）算 UT1，恒星时差几十角秒，
    比的是 ΔT 而不是归算链。喂同一个 ΔT 才能隔离出归算误差。
    """
    ut1 = ut1_jd(tt_jd)
    # 直接把 UT1 当成 UTC 建 Time，并把 delta_ut1_utc 钉成 0，于是 astropy 内部的
    # UT1 与我们的**完全一致**，恒星时不再掺进两套 ΔT 的差别。
    # 代价是 astropy 的 TT 会偏 (ΔT−32.184−闰秒)，但参考系只是弱依赖 TT
    # （几十秒的岁差可忽略），对本比对无影响。
    t = Time(ut1, format='jd', scale='utc')
    t.delta_ut1_utc = 0.0     # 这一句才真正让 astropy 的 UT1 等于传入值
    return t
AS = 206264.806
D2R = np.pi / 180.0

# 与本实现 Atmosphere::default() 一致，也用于喂给 astropy
PRESSURE_HPA = 1010.0
TEMPERATURE_C = 10.0


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'observer_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    gst, hz = [], []
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'GST':
            gst.append(tuple(float(x) for x in p[1:]))
        elif p[0] == 'HZ':
            hz.append([float(x) for x in p[1:5]] + [float(x) for x in p[5:]])

    print('=== 恒星时（对 erfa.gmst06 / erfa.gst06a）===')
    Y = np.array([g[0] for g in gst])
    jd = 2451545.0 + (Y - 2000.0) * 365.25
    # 本实现内部由 TT 推 UT1；这里按同样规则还原一份 UT1 给 ERFA
    # 直接给 ERFA 喂我们算出的 UT1，绕开 astropy 的时间尺度转换
    ut1 = np.array([ut1_jd(x) for x in jd])
    ttv = jd
    m = np.array([erfa.gmst06(ut1[i], 0.0, ttv[i], 0.0) for i in range(len(Y))])
    a = np.array([erfa.gst06a(ut1[i], 0.0, ttv[i], 0.0) for i in range(len(Y))])
    mine_m = np.array([g[1] for g in gst])
    mine_a = np.array([g[2] for g in gst])
    wrap = lambda d: (d + 180) % 360 - 180
    em = np.abs(wrap(mine_m - m / D2R)) * 3600
    ea = np.abs(wrap(mine_a - a / D2R)) * 3600
    for lo, hi in [(-4000, -1000), (-1000, -100), (-100, 100), (100, 1000), (1000, 4000)]:
        s = (Y >= lo) & (Y < hi)
        if s.any():
            print(f'  {lo:>6}..{hi:<6} n={s.sum():4}  GMST 最大 {em[s].max():10.5f}"  '
                  f'GAST 最大 {ea[s].max():10.5f}"')
    print(f'  全部           GMST 最大 {em.max():.5f}"   GAST 最大 {ea.max():.5f}"')

    # 地平坐标**没有**在这里做批量外部对拍，原因如下。
    #
    # 试过用 astropy 的 AltAz 当参考：把 UT1 直接喂成 UTC、再钉住 delta_ut1_utc，
    # 已确认 astropy 内部的 UT1 与我们**完全一致**（0.000 s），但批量跑出来的方位角
    # 在相邻年代之间大幅摆动，与我们的结果差到几十度；而单独挑 J2000 的一个算例
    # 手工比对，两边只差 7"（方位角）与 1"（高度角）。
    #
    # 即：**这个脚手架本身不可靠，我还没查出原因**。查清之前不报它的数字——
    # 报一个自己都不信的数比不报更糟。
    #
    # 因此地平坐标目前的保证只有：
    #   * observer 的单元测试：天顶→alt=90°、赤道上升星→正东、折射量级；
    #   * J2000 单个算例与 astropy 手工比对的 7"/1"。
    # 要真正验收第 4 步的 <0.1"，需要一个可复现的参考（erfa.atioq 走 CIRS，
    # 或离线重实现一遍）。见 docs/plan.md 第 4 步。
    print()
    print('地平坐标：批量外部对拍尚未建立（脚手架不可靠，见脚本注释），未在此报告。')
    return 0


if __name__ == '__main__':
    sys.exit(main())
