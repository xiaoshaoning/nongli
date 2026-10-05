#!/usr/bin/env python3
"""把月相与"用向量独立算一遍"的结果比对。

独立路径：取 ERFA 的月球与地球（→太阳）地心矢量，按**余弦定理在月球顶点**算
相位角，再得照亮比例。与 Rust 侧用 Meeus 的经度公式是两条不同的算法。

用法:  python tools/moonphasecheck.py
"""
import os
import subprocess
import sys

import numpy as np
import erfa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
AU_KM = 149597870.7


def main():
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'moon', '--example', 'phase_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout

    jd, k, i_deg, chi, dist, diam, alt = [], [], [], [], [], [], []
    for line in out.splitlines():
        p = line.split()
        if len(p) != 7:
            continue
        jd.append(float(p[0]))
        k.append(float(p[1]))
        i_deg.append(float(p[2]))
        chi.append(float(p[3]))
        dist.append(float(p[4]))
        diam.append(float(p[5]))
        alt.append(float(p[6]))
    jd = np.array(jd)

    # 独立一遍：余弦定理，顶点在月球
    ks, is_, ds = [], [], []
    for x in jd:
        rm = np.array(erfa.moon98(x, 0.0)[0])                 # 地心月球 (AU, GCRS)
        rs = -np.array(erfa.epv00(x, 0.0)[0][0])              # 地心太阳 (AU)：pvh[0]
        v_me = -rm                                            # 月 → 地
        v_ms = rs - rm                                        # 月 → 日
        cos_i = float(v_me @ v_ms) / (np.linalg.norm(v_me) * np.linalg.norm(v_ms))
        i = np.arccos(np.clip(cos_i, -1, 1))
        is_.append(np.degrees(i))
        ks.append((1 + cos_i) / 2)
        ds.append(np.linalg.norm(rm) * AU_KM)
    ks = np.array(ks)
    is_ = np.array(is_)
    ds = np.array(ds)

    print(f'样本 {len(jd)} 天（{jd.min():.1f} … {jd.max():.1f} JD，约 20 年）')
    e = np.abs(np.array(k) - ks)
    print(f'\n照亮比例: 最大偏差 {e.max():.3e}   中位 {np.median(e):.3e}')
    print(f'  （最差处照亮 {np.array(k)[e.argmax()]:.6f} vs {ks[e.argmax()]:.6f}）')
    ei = np.abs(np.array(i_deg) - is_)
    print(f'相位角  : 最大偏差 {ei.max() * 3600:.3f}"')
    ed = np.abs(np.array(dist) - ds)
    print(f'地心距离: 最大偏差 {ed.max():.6f} km')
    # 视直径应由距离决定：2*asin(R/d)
    R = 1737.4
    dd = np.abs(np.array(diam) - 2 * np.degrees(np.arcsin(R / np.array(dist))))
    print(f'视直径自洽 (对 2asin(R/Δ)): 最大偏差 {dd.max() * 3600:.6f}"')
    return 0


if __name__ == '__main__':
    sys.exit(main())
