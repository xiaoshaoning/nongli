#!/usr/bin/env python3
"""ΔT 整体偏移对农历**结构**的影响。

ΔT 若整体偏 T 秒，所有 朔/中气 的北京日界同时平移。问题是：
  (a) 春节(正月初一)会差几天？
  (b) 更重要的：**闰月判定**会不会跟着变？

若闰月结构对 ΔT 不敏感，那么即便在 ±4000 年（ΔT 有小时级不确定度），
"哪一年、闰几月" 仍然稳定，只有具体的"日"不可定。

用法:  python probe/dtshift.py
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


DT = dtreference.delta_t_from_jd


def day(jde, shift_sec):
    return np.floor(jde - (DT(np.asarray(jde, dtype=float)) + shift_sec) / 86400.0 + 8.0 / 24.0 + 0.5)


def load():
    y0, y1 = -1975, 6026
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'nongli', '--example', 'ephemeris_dump', '--',
         str(y0), str(y1), '1', '1'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    nm, st = {}, {}
    for line in out.splitlines():
        p = line.split()
        if p[0] == 'NM':
            nm[int(p[1])] = float(p[2])
        elif p[0] == 'ST':
            st[int(p[1])] = float(p[2])
    return nm, st


def sui_structure(nm, st, year, shift):
    """返回 (岁内月数, 闰月在该岁内的索引 or None, 正月初一在该岁内的索引)。"""
    # 冬至 j：j ≡ 18 (mod 24)，落在 year 年 12 月
    jd_dec = 2451545.0 + (year - 2000.0) * 365.25
    j = int(round((0.98564736 * (jd_dec - 2451545.0) + 280.46646) / 15.0))
    while j % 24 != 18 % 24:
        j += 1
    # 取最接近 12 月下旬的那个
    cands = [j - 24, j, j + 24]
    j_ws = min(cands, key=lambda x: abs(st[x] - jd_dec))
    ws_day = day(st[j_ws], shift)

    def month_start_on_or_before(d):
        k = int(round((d - 2451550.09766) / 29.530588861))
        while day(nm[k], shift) > d:
            k -= 1
        while day(nm[k + 1], shift) <= d:
            k += 1
        return k

    k0 = month_start_on_or_before(ws_day)
    # 下一个冬至
    j2 = j_ws + 24
    ws2_day = day(st[j2], shift)
    k1 = month_start_on_or_before(ws2_day)
    nmonths = int(k1 - k0)

    if nmonths != 13:
        return nmonths, None, 2 if nmonths == 12 else None

    starts = np.array([day(nm[k0 + i], shift) for i in range(14)])
    has = [False] * 13
    for jj in range(j_ws, j_ws + 24):
        if jj % 2:
            continue
        d = day(st[jj], shift)
        if starts[0] <= d < starts[13]:
            has[int(np.searchsorted(starts, d, 'right') - 1)] = True
    leap = next((i for i in range(1, 13) if not has[i]), None)
    # 正月初一：不计闰月数到第 2 个月
    zy = None
    num, idx, cnt = 11, 0, 0
    for idx in range(13):
        if idx > 0:
            if idx == leap:
                pass
            else:
                num = num % 12 + 1
        if num == 1:
            zy = idx
            break
    return nmonths, leap, zy


def main():
    nm, st = load()
    print(f'样本: {len(nm)} 个朔, {len(st)} 个节气')

    shifts = [0, 10, 60, 600, 3600, 21600]
    years = list(range(-1974, 6025, 7))
    base = {y: sui_structure(nm, st, y, 0) for y in years}

    print(f'\n{"ΔT 偏移":>10} | ' + ' '.join(f'{lo:>6}..{hi:<6}' for lo, hi in
          [(-1974, -500), (-500, 1000), (1000, 2100), (2100, 4000), (4000, 6025)]))
    for s in shifts:
        row = []
        for lo, hi in [(-1974, -500), (-500, 1000), (1000, 2100), (2100, 4000), (4000, 6025)]:
            ys = [y for y in years if lo <= y < hi]
            diff = sum(1 for y in ys if sui_structure(nm, st, y, s) != base[y])
            row.append(f'{100.0*diff/len(ys):5.1f}%')
        print(f'{s:>8}s  | ' + ' '.join(f'{r:>13}' for r in row))

    print('\n(百分比 = 该区间内"岁内月数 / 闰月索引 / 正月索引"发生变化的年份占比)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
