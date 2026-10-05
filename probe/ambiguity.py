#!/usr/bin/env python3
"""农历结果的"唯一性"分析。

一次 朔 或 中气 只要落在北京时间午夜前后 U 秒内，当 ΔT 或星历有 U 秒级不确定度时，
它落在哪一天就不唯一 —— 而 朔 决定月界、中气决定闰月，于是整年的农历都可能翻掉。

本脚本对 2026±4000 的每一个朔、每一个中气，量出它离午夜还有多少秒，
从而给出"给定不确定度 U，哪些年不唯一"。

用法:  python probe/ambiguity.py
"""
import os
import subprocess
import sys

import numpy as np

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dtreference  # noqa: E402  (ΔT 参考实现，与 crates/ephemeris/src/time.rs 保持一致)


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)


def margins(jds):
    """每个时刻距离最近北京午夜还有多少秒。"""
    dt = dtreference.delta_t_from_jd(jds)
    jd_bt = jds - dt / 86400.0 + 8.0 / 24.0
    frac = (jd_bt + 0.5) % 1.0          # 0.0 = 午夜
    return np.minimum(frac, 1.0 - frac) * 86400.0


def main():
    y0, y1 = -1975, 6026
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'nongli', '--example', 'ephemeris_dump', '--',
         str(y0), str(y1), '1', '1'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    nm, zq = [], []
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'NM':
            nm.append(float(p[2]))
        elif p[0] == 'ST':
            j = int(p[1])
            if j % 2 == 0:               # 中气（含冬至）
                zq.append(float(p[2]))
    nm = np.array(nm)
    zq = np.array(zq)
    print(f'样本: 朔 {len(nm)} 个, 中气 {len(zq)} 个, 共 {len(nm)+len(zq)} 个"决定日界"的时刻')

    def year_of(jd):
        return (2000.0 + (jd - 2451545.0) / 365.25).astype(int)

    alljd = np.concatenate([nm, zq])
    m = margins(alljd)
    yrs = year_of(alljd)

    # 按年取最小余量
    order = np.argsort(yrs)
    yrs_s, m_s = yrs[order], m[order]
    uniq, idx = np.unique(yrs_s, return_index=True)
    mins = np.minimum.reduceat(m_s, idx)

    print(f'\n每年"最危险事件"的余量分位数 (秒):')
    for q in [0, 0.1, 0.5, 1, 5, 50]:
        print(f'   {q:5.1f}%   {np.percentile(mins, q):10.1f}')

    print(f'\n全区间余量最小的 15 个事件:')
    worst = np.argsort(m)[:15]
    for i in worst:
        kind = '朔' if i < len(nm) else '中气'
        print(f'   {yrs[i]:6d}  {kind}  余量 {m[i]:8.1f} s')

    print('\n各不确定度 U 下，"该年存在临界事件"的比例:')
    print(f'{"U":>10} | ' + ' '.join(f'{lo:>5}..{hi:<5}' for lo, hi in
          [(-1975, -500), (-500, 500), (500, 1526), (1526, 2100), (2100, 2526), (2526, 4000), (4000, 6026)]))
    for U in [1, 10, 60, 300, 1800, 7200]:
        row = []
        for lo, hi in [(-1975, -500), (-500, 500), (500, 1526), (1526, 2100),
                       (2100, 2526), (2526, 4000), (4000, 6026)]:
            sel = (yrs >= lo) & (yrs < hi)
            if not sel.any():
                row.append('   -  ')
                continue
            ev = m[sel]
            # 该区间内"有临界事件"的年数 / 总年数
            ys = np.unique(yrs[sel])
            bad = sum(1 for y in ys if m[sel & (yrs == y)].min() < U)
            row.append(f'{100.0*bad/len(ys):4.1f}%')
        print(f'{U:>8}s  | ' + ' '.join(f'{r:>11}' for r in row))
    return 0


if __name__ == '__main__':
    sys.exit(main())
