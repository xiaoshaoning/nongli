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

    # 黄道→赤道：与 ERFA 的 (黄道→ICRS→赤道) 复合比较。
    # 用 ERFA 自己的真交角，隔离代数与交角模型。
    print('\n黄道 → 赤道（用 ERFA 的真交角，只查球面三角代数）:')
    rm = erfa.ecm06(2451545.0, 0.0)          # ICRS → 黄道
    rmt = erfa.ecm06(2451545.0, 0.0).T       # 黄道 → ICRS
    eps_t = erfa.obl06(2451545.0, 0.0)
    for lam, bet, ra, dec in EQ:
        v = np.array([np.cos(bet * D2R) * np.cos(lam * D2R),
                      np.cos(bet * D2R) * np.sin(lam * D2R),
                      np.sin(bet * D2R)])
        w = rmt @ v                              # 黄道 → ICRS
        # ICRS → 真赤道（经岁差章动）
        pnm = erfa.pnm06a(2451545.0, 0.0)
        u = pnm @ w
        ra_ref = np.degrees(np.arctan2(u[1], u[0])) % 360
        dec_ref = np.degrees(np.arcsin(np.clip(u[2], -1, 1)))
        # 我们的 α、δ 是相对**当日真赤道**，而 ICRS→真赤道的 pnm 掺了岁差；
        # 在 J2000 处两者一致，故这里直接比。
        d_ra = (ra - ra_ref + 180) % 360 - 180
        print(f'  λ={lam:6.1f} β={bet:+6.1f}  Δα={d_ra*3600:+9.4f}"  Δδ={(dec-dec_ref)*3600:+9.4f}"')
    return 0


if __name__ == '__main__':
    sys.exit(main())
