#!/usr/bin/env python3
"""Generate crates/ephemeris/src/frames_tables.rs.

两件事，都用**拟合**而非凭记忆抄系数，并且都能对 ERFA 验证：

1. **IAU2006 岁差（Fukushima–Williams 角）**：对 `erfa.pfw06` 拟合 5 次多项式。
   这些角本就是 5 次多项式，故拟合即是复原，系数可对到 1e-10 角秒。

2. **章动 Δψ / Δε**：对 `erfa.nut06a`（IAU2006/2000A）做**正交匹配追踪**，
   用五个基本引数 (l, l', F, D, Ω) 的整数组合作基。
   IAU2000A 有约 1300 项，凭记忆抄是不可能的；拟合出来的是同一族基函数，
   残差可直接量，项数按需要的精度定。

三个基本引数多项式取自 IERS Conventions (2003) 的 IAU2000A 定义。

Run:  python tools/gen_frames.py
"""
import os
import sys

import numpy as np
import erfa

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, os.pardir, 'crates', 'ephemeris', 'src', 'frames_tables.rs')
AS2R = np.pi / (180.0 * 3600.0)   # 角秒 → 弧度
D2R = np.pi / 180.0

# ---------------------------------------------------------------- 岁差

PREC_T = 100.0        # 拟合区间 ±100 儒略世纪
PREC_N = 2001


def fit_precession():
    T = np.linspace(-PREC_T, PREC_T, PREC_N)
    jd = 2451545.0 + T * 36525.0
    raw = np.array([erfa.pfw06(x, 0.0) for x in jd])  # gamb, phib, psib, epsa (rad)
    names = ['GAMB', 'PHIB', 'PSIB', 'EPSA']
    out = {}
    for k, name in enumerate(names):
        # 用缩放变量拟合以改善条件数，再换回 T 的幂基
        c = np.polyfit(T / PREC_T, raw[:, k], 5)
        # c 是 (T/100) 的降幂系数 → 换成 T 的幂基
        cs = np.array([c[i] / PREC_T ** (5 - i) for i in range(6)])
        out[name] = cs
    return out


def check_precession(coefs):
    T = np.linspace(-200, 200, 2001)
    jd = 2451545.0 + T * 36525.0
    ref = np.array([erfa.pfw06(x, 0.0) for x in jd])
    worst = 0.0
    for k, name in enumerate(['GAMB', 'PHIB', 'PSIB', 'EPSA']):
        d = np.abs(np.polyval(coefs[name], T) - ref[:, k]).max() / AS2R
        worst = max(worst, d)
    return worst


# ---------------------------------------------------------------- 章动基

def fundamental_arguments(T):
    """IAU2000A 的五个基本引数（弧度）。T 为儒略世纪 (TT from J2000)。

    系数单位为角秒，取自 IERS Conventions (2003)。
    """
    l = (134.96340251 * 3600 + 1717915923.2178 * T + 31.8792 * T**2
         + 0.051635 * T**3 - 0.00024470 * T**4) * AS2R
    lp = (357.52910918 * 3600 + 129596581.0481 * T - 0.5532 * T**2
          + 0.000136 * T**3 - 0.00001149 * T**4) * AS2R
    f = (93.27209062 * 3600 + 1739527262.8478 * T - 12.7512 * T**2
         - 0.001037 * T**3 + 0.00000417 * T**4) * AS2R
    d = (297.85019547 * 3600 + 1602961601.2090 * T - 6.3706 * T**2
         + 0.006593 * T**3 - 0.00003169 * T**4) * AS2R
    om = (125.04455501 * 3600 - 6962890.5431 * T + 7.4722 * T**2
          + 0.007702 * T**3 - 0.00005939 * T**4) * AS2R
    return l, lp, f, d, om


def candidate_pool(limit):
    """(i, j, k, m, n) 的候选，每个引数组只保留一个符号。

    `sin(−θ) = −sin θ`、`cos(−θ) = cos θ`，所以 `c` 与 `−c` 张成同一个平面；
    两个都放进基里会让正规方程奇异。只保留首个非零分量为正的那个。
    """
    rng = range(-limit, limit + 1)
    out = []
    for i in rng:
        for j in rng:
            for k in rng:
                for m in rng:
                    for n in rng:
                        c = (i, j, k, m, n)
                        if c == (0, 0, 0, 0, 0):
                            continue
                        first = next(v for v in c if v != 0)
                        if first > 0:
                            out.append(c)
    return out


def angle_of(cands, args):
    """候选基的引数矩阵，形状 (n_cand, n_samp)。"""
    l, lp, f, d, om = args
    th = np.empty((len(cands), len(l)))
    for c, (i, j, k, m, n) in enumerate(cands):
        th[c] = i * l + j * lp + k * f + m * d + n * om
    return th


def select_and_fit(y, args, cands, n_terms, T, chunk=150, min_sep=1e-3):
    """先按投影强度排序取前 `n_terms` 项，再对它们整体最小二乘。

    每项带 4 个系数：`A·sin + B·cos + T·(A'·sin + B'·cos)`——T 依赖不可省，
    否则外推到几千年会差到角秒（主项 A' 是 −174.7 mas/世纪）。

    在近正交的基上单次投影已约等于最小二乘系数，故"先排序再拟合"比逐轮贪心
    快得多，也不会像贪心那样去拟合采样网格的噪声。

    **采样必须避开混叠**：基里含 27 天的月球周期，采样间隔必须 ≲1 天，
    否则高频项会折到低频，拟合出来的系数在区间外会爆掉。
    """
    ncol = 4
    score = np.empty(len(cands))
    for s0 in range(0, len(cands), chunk):
        cc = cands[s0:s0 + chunk]
        th = angle_of(cc, args)
        Sn, Cs = np.sin(th), np.cos(th)
        score[s0:s0 + chunk] = ((Sn @ y) ** 2 + (Cs @ y) ** 2
                                + (T * Sn @ y) ** 2 + (T * Cs @ y) ** 2)
    # 只保留频率互相分开的候选。频率挨得太近（小于采样总跨度的分辨率）的两项，
    # 系数无法分辨，最小二乘会给出巨大的、互相抵消的值，一外推就爆。
    # 频率 = θ 的 T 线性系数，单位 cycles/年。
    rate = np.array([(i * 1717915923.2178 + j * 129596581.0481 + k * 1739527262.8478
                      + m * 1602961601.2090 - n * 6962890.5431) / (1296000.0 * 100.0)
                     for (i, j, k, m, n) in cands])
    order = np.argsort(score)[::-1]
    keep = []
    for idx in order:
        if all(abs(rate[idx] - rate[j]) > min_sep for j in keep):
            keep.append(idx)
        if len(keep) >= n_terms:
            break
    sel = [cands[i] for i in keep]
    th = angle_of(sel, args)

    # 用正规方程：列数 ~4·n_terms，用 lstsq 做 SVD 太贵
    nsamp, n = len(y), 4 * len(sel)
    AtA = np.zeros((n, n))
    Aty = np.zeros(n)
    for s0 in range(0, nsamp, 4000):
        sl = slice(s0, min(s0 + 4000, nsamp))
        Sn, Cs = np.sin(th[:, sl]), np.cos(th[:, sl])
        Ac = np.empty((len(T[sl]), n))
        Ac[:, 0::4], Ac[:, 1::4] = Sn.T, Cs.T
        Ac[:, 2::4] = (Sn * T[sl]).T
        Ac[:, 3::4] = (Cs * T[sl]).T
        AtA += Ac.T @ Ac
        Aty += Ac.T @ y[sl]
    # 用 lstsq 而非 solve：即便仍有轻度共线也不会炸
    coef = np.linalg.lstsq(AtA, Aty, rcond=None)[0]
    resid = np.empty(nsamp)
    for s0 in range(0, nsamp, 4000):
        sl = slice(s0, min(s0 + 4000, nsamp))
        Sn, Cs = np.sin(th[:, sl]), np.cos(th[:, sl])
        Ac = np.empty((len(T[sl]), n))
        Ac[:, 0::4], Ac[:, 1::4] = Sn.T, Cs.T
        Ac[:, 2::4] = (Sn * T[sl]).T
        Ac[:, 3::4] = (Cs * T[sl]).T
        resid[sl] = y[sl] - Ac @ coef
    return sel, coef, resid


def evaluate(res, T):
    """在给定 T 上求拟合值。"""
    l, lp, f, d, om = fundamental_arguments(T)
    out = {}
    for name in ('DPSI', 'DEPS'):
        sel, coef = res[name]
        th = angle_of(sel, (l, lp, f, d, om))
        y = np.zeros_like(T)
        for t in range(len(sel)):
            y += (coef[4*t] + coef[4*t+2]*T) * np.sin(th[t])                + (coef[4*t+1] + coef[4*t+3]*T) * np.cos(th[t])
        out[name] = y
    return out


def sample_windows(span_years, n_windows, win_days, dt=1.0):
    """在长跨度上铺若干**短窗**，窗内细步长。

    两个约束互相冲突，只有这样才能同时满足：
      * **不混叠**要求局部步长 ≲1 天（基里有 27 天的月球周期）；
      * **可分辨**要求总跨度足够长，否则频率挨得近的两项无法分开，
        系数会互相抵消，一旦外推就爆掉（实测外推 10 世纪能差到 3.6e3 角秒）。
    总跨度取 ±span_years，窗内 1 天步长。
    """
    # 中心必须**抖动**：等距排列会以 1/间距 的频率产生混叠（间距 40 年 →
    # 0.025 cycles/年），让本来可分辩的两个频率重新纠缠在一起。
    rng = np.random.default_rng(20260217)
    centers = np.sort(rng.uniform(-span_years, span_years, n_windows))
    offs = np.arange(0.0, win_days, dt)
    days = np.concatenate([c * 365.25 + offs for c in centers])
    return days / 36525.0


def fit_nutation(n_terms, limit, span_years):
    """在 ±span_years 年内按窗采样拟合。"""
    T = sample_windows(span_years, n_windows=400, win_days=200.0)
    jd = 2451545.0 + T * 36525.0
    dp, de = np.array([erfa.nut06a(x, 0.0) for x in jd]).T
    args = fundamental_arguments(T)
    cands = candidate_pool(limit)
    print(f'  候选基 {len(cands)} 个，样本 {len(T)} 点', flush=True)
    res = {}
    for name, y in [('DPSI', dp), ('DEPS', de)]:
        sel, coef, resid = select_and_fit(y, args, cands, n_terms, T)
        rms = np.sqrt((resid ** 2).mean()) / AS2R
        mx = np.abs(resid).max() / AS2R
        print(f'  {name}: 用 {len(sel)} 项，拟合区间内残差 rms={rms:.5f}"  最大={mx:.5f}"')
        res[name] = (sel, coef)
    return res


def nutation_error(res, T):
    jd = 2451545.0 + T * 36525.0
    dp, de = np.array([erfa.nut06a(x, 0.0) for x in jd]).T
    fit = evaluate(res, T)
    return {'DPSI': (fit['DPSI'] - dp) / AS2R, 'DEPS': (fit['DEPS'] - de) / AS2R}


# ---------------------------------------------------------------- 输出

def emit(prec, res):
    L = []
    L.append('// @generated by tools/gen_frames.py -- do not edit by hand.')
    L.append('//')
    L.append('// 黄赤交角：IAU2006 Fukushima–Williams 的 EPSA 角。**对 erfa.pfw06 拟合**得到')
    L.append('// 5 次多项式，故系数是被复原的而不是凭记忆抄的。单位角秒，T 为儒略世纪，降幂。')
    L.append('//')
    L.append('// 另三个 FW 角（γ̄、φ̄、ψ̄）已在生成时拟合并验证（±200 世纪内差 2.3e-8″），')
    L.append('// 但没有消费者，故不输出——只有把"当日"化到 ICRS 时才需要。见 src/frames.rs。')
    L.append('//')
    L.append('// 章动**不在本文件里**：缩写式写在 src/frames.rs。')
    L.append('// 本脚本里有一段"对 erfa.nut06a 拟合更好章动级数"的代码，**结论是失败**：')
    L.append('// 拟合出的 200 项级数在 ±1 世纪（0.2 天步长、36.5 万点）是 Δψ 0.75″/rms 0.22″，')
    L.append('// 而现有的 4 项缩写式是 0.34″/0.12″——更差。精度在约 0.3″ 处饱和，加项数无用。')
    L.append('// 那段代码保留，是为了记录结论，并让将来拿到 IAU2000A 原始表时能对拍。')
    L.append('// 详见 docs/plan.md 第 3 步。')
    L.append('')
    # 只输出 EPSA（黄赤交角）；另三个 FW 角无消费者，见 src/frames.rs。
    for name in ['EPSA']:
        c = prec[name]
        L.append(f'/// `{name}` 的多项式系数，角秒，T 的降幂 (T⁵ … T⁰)。')
        L.append(f'pub static {name}: [f64; 6] = [')
        for v in c:
            L.append(f'    {float(v / AS2R)!r},')
        L.append('];')
        L.append('')
    # 章动**不输出**：拟合出的级数在 ±1 世纪反而比缩写式差（见文件头与 plan）。
    # 这段拟合代码保留，是为了记录结论、并让将来拿到 IAU2000A 原始表时能对拍。
    with open(OUT, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(L))
    print('wrote', OUT)


def curve(limit, span_years, counts):
    """精度-项数曲线：排序只做一次，对不同的项数分别拟合评估。"""
    T = sample_windows(span_years, n_windows=400, win_days=200.0)
    jd = 2451545.0 + T * 36525.0
    dp, de = np.array([erfa.nut06a(x, 0.0) for x in jd]).T
    args = fundamental_arguments(T)
    cands = candidate_pool(limit)
    print()
    print('精度 vs 项数（拟合跨度 ±%d 年，候选 %d）' % (span_years, len(cands)))
    hdr = '  项数   ' + ''.join(f'±{h:<7}' for h in [1, 2, 10, 40])
    print(hdr)
    for n in counts:
        row = f'  {n:5}   '
        res = {}
        for name, y in [('DPSI', dp), ('DEPS', de)]:
            sel, coef, _ = select_and_fit(y, args, cands, n, T)
            res[name] = (sel, coef)
        for h in [1, 2, 10, 40]:
            Tt = np.linspace(-h, h, 401)
            e = nutation_error(res, Tt)
            row += f'{np.abs(e["DPSI"]).max():<8.4f}'
        print(row)
    return res


def main():
    print('岁差：对 erfa.pfw06 拟合 5 次多项式')
    prec = fit_precession()
    print(f'  T∈[−200,200] 上最大偏差 {check_precession(prec):.3e}"')
    print()
    print('章动：正交匹配追踪 erfa.nut06a')
    # 失败尝试的记录；结果不输出。详见文件头。
    curve(limit=2, span_years=2000.0, counts=[200])
    res = None
    print()
    emit(prec, res)
    return 0


if __name__ == '__main__':
    sys.exit(main())
