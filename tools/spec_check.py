#!/usr/bin/env python3
"""GB/T 33661-2017 第 **5.2** 条的验收：

    朔和节气的北京时间计算精度应达到 1 s。

这是整份规范里**唯一一条定量要求**，此前只用角秒量过精度，没换成"秒"。
真值是 JPL DE421（覆盖 1900–2053）。两边用同一套 ΔT，所以作差反映的是
**模型/理论**的误差；ΔT 自身的不确定度另见 README 那张表。

用法:  python tools/spec_check.py
需要:  kernels/de421.bsp（取法见 tools/truthcheck.py 头部）
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
R2D = 180.0 / np.pi


def wrap180(x):
    return (x + 180.0) % 360.0 - 180.0


def main():
    if not os.path.exists(KERNEL):
        print(f'缺少 {KERNEL} —— 跳过规范 5.2 的验收。取法见 tools/truthcheck.py 头部。')
        return 0

    import erfa
    from jplephem.spk import SPK
    K = SPK.open(KERNEL)

    def earth(j):
        return K[(0, 3)].compute(j) + K[(3, 399)].compute(j)

    def sun(j):
        return K[(0, 10)].compute(j)

    def moon(j):
        return K[(0, 3)].compute(j) + K[(3, 301)].compute(j)

    # ---- 真值侧的视黄经。口径与 crates/ephemeris/src/{sun,moon}.rs 逐条对应。
    def sun_lon(j, apparent=True):
        g = sun(j) - earth(j)
        r = np.linalg.norm(g)
        w = erfa.ecm06(j, 0.0) @ (g / r)
        lam = np.degrees(np.arctan2(w[1], w[0])) % 360.0
        lam -= 20.4898 / (r / AU_KM) / 3600.0        # 周年光行差，与 sun.rs 同号
        if apparent:
            lam += erfa.nut06a(j, 0.0)[0] * R2D
        return lam % 360.0

    def moon_lon(j, apparent=True):
        r = np.linalg.norm(moon(j) - earth(j))
        js = j - r / C_KM_S / 86400.0                # 光行时
        g = moon(js) - earth(js)
        w = erfa.ecm06(j, 0.0) @ (g / np.linalg.norm(g))
        lam = np.degrees(np.arctan2(w[1], w[0])) % 360.0
        r_sun = np.linalg.norm(sun(j) - earth(j))
        lam -= 20.4898 / (r_sun / AU_KM) / 3600.0
        if apparent:
            lam += erfa.nut06a(j, 0.0)[0] * R2D
        return lam % 360.0

    def solve(f, guess):
        """对 f(jd)=0 做牛顿迭代（数值导数）。"""
        jd = guess
        for _ in range(12):
            h = 0.02
            y0, y1, y2 = f(jd), f(jd + h), f(jd - h)
            d = (y1 - y2) / (2 * h)
            step = y0 / d
            jd -= step
            if abs(step) < 1e-10:
                break
        return jd

    # ---- 取本实现的输出
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'nongli',
         '--example', 'spec_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    nm = [(int(p[1]), float(p[2])) for p in
          (l.split() for l in out.splitlines()) if p and p[0] == 'NEWMOON']
    tm = [(int(p[1]), float(p[2])) for p in
          (l.split() for l in out.splitlines()) if p and p[0] == 'TERM']

    d_nm, d_tm = [], []
    for _, jd in nm:
        f = lambda j: wrap180(moon_lon(j, False) - sun_lon(j, False))  # noqa: E731
        j = solve(f, jd)
        d_nm.append((jd - j) * 86400.0)
    for j, jd in tm:
        target = (15.0 * j) % 360.0
        f = lambda x: wrap180(sun_lon(x, True) - target)  # noqa: E731
        j2 = solve(f, jd)
        d_tm.append((jd - j2) * 86400.0)

    d_nm, d_tm = np.array(d_nm), np.array(d_tm)
    print(f'GB/T 33661-2017 §5.2 验收（真值 JPL DE421，1900–2050）\n')
    print(f'{"":10} {"个数":>6} {"中位":>9} {"中位|差|":>11} {"最大":>9}   判据')
    for name, d in [('朔', d_nm), ('节气', d_tm)]:
        ok = '✅' if np.abs(d).max() <= 1.0 else '❌'
        print(f'{name:10} {len(d):>6} {np.median(d):+8.3f}s {np.median(abs(d)):10.3f}s '
              f'{np.abs(d).max():8.3f}s   ≤1 s {ok}')

    print('\n按年代（看有没有趋势）：')
    yrs = np.array([1900.0 + (jd - 2415020.5) / 365.25 for _, jd in nm])
    yrs_t = np.array([1900.0 + (jd - 2415020.5) / 365.25 for _, jd in tm])
    for lo, hi in [(1900, 1930), (1930, 1960), (1960, 1990), (1990, 2020), (2020, 2050)]:
        m, mt = (yrs >= lo) & (yrs < hi), (yrs_t >= lo) & (yrs_t < hi)
        print(f'  {lo}-{hi}: 朔 中位 {np.median(d_nm[m]):+7.3f}s 最大 {np.abs(d_nm[m]).max():6.3f}s'
              f'   节气 中位 {np.median(d_tm[mt]):+7.3f}s 最大 {np.abs(d_tm[mt]).max():6.3f}s')

    print('''
说明：
* 两边用**同一套 ΔT**，所以上面量的是模型/理论的误差——正是 §5.2 要的。
  ΔT 自身的不确定度（±4000 年可达小时量级）是另一回事，README 有专门一节。
* 判据按规范取 **1 s**（未设余量）。实测余量见上表。
''')
    return 0 if max(np.abs(d_nm).max(), np.abs(d_tm).max()) <= 1.0 else 1


if __name__ == '__main__':
    sys.exit(main())
