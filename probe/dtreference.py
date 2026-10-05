#!/usr/bin/env python3
"""ΔT 参考实现 —— **必须与 src/jd.rs 保持一致**。

工具/探针脚本共用这一份，避免各脚本各写一套导致比对出现假的差异。
改 src/jd.rs 的 ΔT 时，这里要同步改。

* 1973–2025：IERS 实测值（内建小表，与 src/jd.rs 的 DELTA_T_OBSERVED 同源同值）
* 其余年份：Espenak & Meeus (2006) 分段多项式，在表端做**增量衔接**以保证连续

接口刻意分成两个名字，避免"把儒略日当年份"这类静默错误：

* :func:`delta_t_from_year` —— 参数是**小数年**
* :func:`delta_t_from_jd`   —— 参数是**儒略日 (TT)**
"""

import numpy as np
import os

# 与 src/jd.rs 的 DELTA_T_OBSERVED 完全一致
DELTA_T_OBSERVED = [
    (1973.0, 43.38), (1974.0, 44.48), (1975.0, 45.48), (1976.0, 46.46), (1977.0, 47.52), (1978.0, 48.53),
    (1979.0, 49.59), (1980.0, 50.54), (1981.0, 51.38), (1982.0, 52.17), (1983.0, 52.96), (1984.0, 53.79),
    (1985.0, 54.34), (1986.0, 54.87), (1987.0, 55.32), (1988.0, 55.82), (1989.0, 56.30), (1990.0, 56.86),
    (1991.0, 57.57), (1992.0, 58.31), (1993.0, 59.12), (1994.0, 59.98), (1995.0, 60.79), (1996.0, 61.63),
    (1997.0, 62.30), (1998.0, 62.97), (1999.0, 63.47), (2000.0, 63.83), (2001.0, 64.09), (2002.0, 64.30),
    (2003.0, 64.47), (2004.0, 64.57), (2005.0, 64.69), (2006.0, 64.85), (2007.0, 65.15), (2008.0, 65.46),
    (2009.0, 65.78), (2010.0, 66.07), (2011.0, 66.32), (2012.0, 66.60), (2013.0, 66.91), (2014.0, 67.28),
    (2015.0, 67.64), (2016.0, 68.10), (2017.0, 68.59), (2018.0, 68.97), (2019.0, 69.22), (2020.0, 69.36),
    (2021.0, 69.36), (2022.0, 69.29), (2023.0, 69.20), (2024.0, 69.18), (2025.0, 69.14),
]

_T0, _T1 = DELTA_T_OBSERVED[0][0], DELTA_T_OBSERVED[-1][0]


def observed(y):
    y = min(max(y, _T0), _T1)
    i = min(int(y - _T0), len(DELTA_T_OBSERVED) - 2)
    (y0, v0), (y1, v1) = DELTA_T_OBSERVED[i], DELTA_T_OBSERVED[i + 1]
    return v0 + (v1 - v0) * (y - y0) / (y1 - y0)


def series(y):
    """Espenak & Meeus (NASA, 2006) 的绝对值。"""
    if y < -500.0:
        u = (y - 1820.0) / 100.0
        return -20 + 32 * u * u
    if y < 500.0:
        u = y / 100.0
        return (10583.6 - 1014.41 * u + 33.78311 * u**2 - 5.952053 * u**3
                - 0.1798452 * u**4 + 0.022174192 * u**5 + 0.0090316521 * u**6)
    if y < 1600.0:
        u = (y - 1000.0) / 100.0
        return (1574.2 - 556.01 * u + 71.23472 * u**2 + 0.319781 * u**3
                - 0.8503463 * u**4 - 0.005050998 * u**5 + 0.0083572073 * u**6)
    if y < 1700.0:
        t = y - 1600.0
        return 120 - 0.9808 * t - 0.01532 * t * t + t**3 / 7129.0
    if y < 1800.0:
        t = y - 1700.0
        return 8.83 + 0.1603 * t - 0.0059285 * t * t + 0.00013336 * t**3 - t**4 / 1174000.0
    if y < 1860.0:
        t = y - 1800.0
        return (13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3
                - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6
                + 0.000000000875 * t**7)
    if y < 1900.0:
        t = y - 1860.0
        return (7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3
                - 0.0004473624 * t**4 + t**5 / 233174.0)
    if y < 1920.0:
        t = y - 1900.0
        return -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4
    if y < 1941.0:
        t = y - 1920.0
        return 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3
    if y < 1961.0:
        t = y - 1950.0
        return 29.07 + 0.407 * t - t * t / 233.0 + t**3 / 2547.0
    if y < 1986.0:
        t = y - 1975.0
        return 45.45 + 1.067 * t - t * t / 260.0 - t**3 / 718.0
    if y < 2005.0:
        t = y - 2000.0
        return (63.86 + 0.3345 * t - 0.060374 * t**2 + 0.0017275 * t**3
                + 0.000651814 * t**4 + 0.00002373599 * t**5)
    if y < 2050.0:
        t = y - 2000.0
        return 62.92 + 0.32217 * t + 0.005589 * t * t
    if y < 2150.0:
        return -20 + 32 * ((y - 1820.0) / 100.0) ** 2 - 0.5628 * (2150.0 - y)
    u = (y - 1820.0) / 100.0
    return -20 + 32 * u * u


def _scalar_year(y):
    if _T0 <= y <= _T1:
        return observed(y)
    anchor = _T0 if y < _T0 else _T1
    return observed(anchor) + (series(y) - series(anchor))


def delta_t_from_year(decimal_year):
    """**小数年** → ΔT (秒)。与 src/jd.rs::delta_t_seconds 同逻辑。

    注意参数是"年"不是"儒略日"。要传儒略日用 :func:`delta_t_from_jd`。
    """
    if isinstance(decimal_year, np.ndarray):
        return np.vectorize(_scalar_year)(decimal_year)
    return _scalar_year(decimal_year)


def delta_t_from_jd(jd):
    """**儒略日 (TT)** → ΔT (秒)。标量或 numpy 数组均可。"""
    y = 2000.0 + (np.asarray(jd, dtype=float) - 2451545.0) / 365.25
    out = delta_t_from_year(y)
    return out if isinstance(out, np.ndarray) else float(out)


if __name__ == '__main__':
    import subprocess

    # 1) 与 IERS 实测值自检
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    worst = 0.0
    for y in range(1973, 2026):
        t = Time(f'{y}-01-01', scale='utc')
        obs = float((t.tt.jd - t.ut1.jd) * 86400.0)
        worst = max(worst, abs(delta_t_from_year(float(y)) - obs))
    # J2000 = 2000.0 年，用日输入的接口必须一致
    assert abs(delta_t_from_jd(2451545.0) - delta_t_from_year(2000.0)) < 1e-9
    print(f'delta_t_from_jd(2451545.0) = {float(delta_t_from_jd(2451545.0)):.3f} s  (= 2000 年)')
    print(f'与 IERS 实测的最大偏差: {worst:.3f} s')
    print(f'E-M 多项式若单独使用，最大偏差: '
          f'{max(abs(series(y) - delta_t_from_year(float(y))) for y in range(1973, 2026)):.2f} s')

    # 2) 与 src/jd.rs 逐点核对（跨语言同一份知识，必须能机器校验）
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '--example', 'dt_dump'],
        cwd=root, capture_output=True, text=True, check=True,
    ).stdout
    n, bad = 0, []
    for line in out.splitlines():
        tag, x, want = line.split()
        x, want = float(x), float(want)
        if tag == 'DT':
            got = float(delta_t_from_year(x))
        else:
            got = float(delta_t_from_year(2000.0 + (x - 2451545.0) / 365.25))
        n += 1
        if abs(got - want) > 1e-6:
            bad.append((tag, x, want, got))
    print(f'与 src/jd.rs 核对 {n} 点，不一致 {len(bad)} 点')
    for tag, x, want, got in bad[:10]:
        print(f'  BAD {tag} {x}: rust={want:.6f} py={got:.6f}')
    raise SystemExit(1 if bad else 0)
