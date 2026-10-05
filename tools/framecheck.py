#!/usr/bin/env python3
"""把参考系量（章动、黄赤交角、黄道→赤道）与 ERFA 对拍。

这是**同口径**的比较：双方都直接给出 Δψ/Δε/ε_A，不涉及月球理论的黄道定义。
用法:  python tools/framecheck.py
"""
import os
import subprocess
import sys

import numpy as np
import erfa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
AS = 206264.806
D2R = np.pi / 180.0

from astropy.time import Time
t0 = Time(2451545.0, format='jd', scale='tt')


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'frame_check', '--', '-4000', '4000'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    Y, dp, de, eps, teps, EQ = [], [], [], [], [], []
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == 'NUT':
            y = float(p[1])
            Y.append(y)
            dp.append(float(p[2]))
            de.append(float(p[3]))
            eps.append(float(p[4]))
            teps.append(float(p[5]))
        elif p[0] == 'EQ':
            EQ.append(tuple(float(x) for x in p[1:]))
    Y = np.array(Y)
    jd = 2451545.0 + (Y - 2000.0) * 365.25

    ref = np.array([erfa.nut06a(x, 0.0) for x in jd]) * AS
    obl = np.array([erfa.obl06(x, 0.0) for x in jd]) / D2R

    print(f'样本 {len(Y)} 点，{Y.min():.0f} … {Y.max():.0f} 年')
    print(f'\n{"区间":>16} | {"Δψ 最大":>10} {"Δψ rms":>10} | {"Δε 最大":>10} | {"ε_A 最大":>12}')
    for lo, hi in [(-4000, -1000), (-1000, -100), (-100, 100), (100, 1000), (1000, 4000)]:
        m = (Y >= lo) & (Y < hi)
        if not m.any():
            continue
        e_psi = np.abs(np.array(dp)[m] - ref[m, 0])
        e_eps = np.abs(np.array(de)[m] - ref[m, 1])
        e_obl = np.abs(np.array(eps)[m] - obl[m])
        print(f'{lo:>7}..{hi:<7} | {e_psi.max():9.4f}" {np.sqrt((e_psi**2).mean()):9.4f}" '
              f'| {e_eps.max():9.4f}" | {e_obl.max():11.3e}°')

    # 黄道 → 赤道**没有**在这里对拍，有意为之。
    #
    # 试过两种外部对拍，都在"平/真黄道"的口径上翻车：把 apparent() 之后的
    # （真分点）坐标喂进 erfa.ecm06（平黄道）再叠加 pnm06a，章动被加了两次，
    # 差 ~12″——那是脚本口径错，不是实现错，但报出来会误导。
    #
    # 目前该转换由 frames.rs 的单元测试保证：春分点 (λ=0,β=0) → (0,0)、
    # 夏至点 (λ=90,β=0) → δ=+ε、秋分点 (λ=180,β=0) → δ=0，均为精确值。
    # 要外部对拍，得先把"平/真"口径写清楚，留到第 4 步与地平坐标一起做。
    print()
    print('黄道 -> 赤道：见上述说明，暂由 frames.rs 的单元测试保证，未在此对拍。')
    return 0


if __name__ == '__main__':
    sys.exit(main())
