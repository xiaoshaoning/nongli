#!/usr/bin/env python3
"""`Elp2000-82b.js` 的逐行 Python 对照实现。

用途：排查我们自己的 Rust 实现与参考实现的分歧。刻意**照抄**参考实现的写法
（包括对 5 阶/2 阶引数的分工），不做任何"改进"——它就是判据。
"""
import json
import sys
import numpy as np

JSON = sys.argv[1] if len(sys.argv) > 1 else \
    'third_party/elp2000-82b/ELP2000-82b.json'
d = json.load(open(JSON, encoding='utf-8'))

ARG = [
    [785939.8782, 1732559343.3328, -6.870, 0.006604, -0.00003169],
    [300071.6518, 14643420.3304, -38.2639, -0.045047, 0.00021301],
    [450160.3265, -6967919.8851, 6.3593, 0.007625, -0.00003586],
    [361679.1880, 129597742.3016, -0.0202, 0.000009, 0.00000015],
    [370574.4136, 1161.2283, 0.5327, -0.000138, 0.0],
]
PLANET = [
    [908103.25986, 538101628.68898], [655127.28305, 210664136.43355],
    [361679.22059, 129597742.2758], [1279559.78866, 68905077.59284],
    [123665.34212, 10925660.42861], [180278.89694, 4399609.65932],
    [1130598.01841, 1542481.19393], [1095655.19575, 786550.32074],
]
RAD = np.pi / 648000.0
CORR = {
    'deltanu': 0.55604 * RAD / ARG[0][1],
    'deltaE': 0.01789 * RAD,
    'deltaGamma': -0.08066 * RAD,
    'deltanp': -0.06424 * RAD / ARG[0][1],
    'deltaep': -0.12879 * RAD,
}
DTASM = 2.0 * 0.002571881335 / (3.0 * 0.074801329518)
AM = 0.074801329518
PREC = 5029.0966 - 0.0316


def polys(t, n):
    return [sum(ARG[i][j] * t ** j for j in range(n)) for i in range(5)]


def delargs(p):
    return {'D': p[0] - p[3] + 648000.0, 'LP': p[3] - p[4],
            'L': p[0] - p[1], 'F': p[0] - p[2]}


def main_sin(delargs, data, cosine=False):
    s = 0.0
    for r in data:
        A = r['A']
        if cosine:
            A += -(2.0 / 3.0) * r['A'] * CORR['deltanu']
        A += (r['B1'] + DTASM * r['B5']) * (CORR['deltanp'] - AM * CORR['deltanu'])
        A += r['B2'] * CORR['deltaGamma'] + r['B3'] * CORR['deltaE'] + r['B4'] * CORR['deltaep']
        arg = (r['i1'] * delargs['D'] + r['i2'] * delargs['LP']
               + r['i3'] * delargs['L'] + r['i4'] * delargs['F']) * RAD
        s += A * (np.cos(arg) if cosine else np.sin(arg))
    return s


def nonplan(prec, delargs, data):
    s = 0.0
    for r in data:
        s += r['A'] * np.sin((r['i1'] * prec + r['i2'] * delargs['D'] + r['i3'] * delargs['LP']
                              + r['i4'] * delargs['L'] + r['i5'] * delargs['F']) * RAD
                             + r['phi'] * np.pi / 180.0)
    return s


def plan1(pl, delargs, data):
    s = 0.0
    for r in data:
        s += r['A'] * np.sin((r['i1'] * pl[0] + r['i2'] * pl[1] + r['i3'] * pl[2] + r['i4'] * pl[3]
                              + r['i5'] * pl[4] + r['i6'] * pl[5] + r['i7'] * pl[6] + r['i8'] * pl[7]
                              + r['i9'] * delargs['D'] + r['i10'] * delargs['L']
                              + r['i11'] * delargs['F']) * RAD + r['phi'] * np.pi / 180.0)
    return s


def plan2(pl, delargs, data):
    s = 0.0
    for r in data:
        s += r['A'] * np.sin((r['i1'] * pl[0] + r['i2'] * pl[1] + r['i3'] * pl[2] + r['i4'] * pl[3]
                              + r['i5'] * pl[4] + r['i6'] * pl[5] + r['i7'] * pl[6]
                              + r['i8'] * delargs['D'] + r['i9'] * delargs['LP']
                              + r['i10'] * delargs['L'] + r['i11'] * delargs['F']) * RAD
                             + r['phi'] * np.pi / 180.0)
    return s


def spherical(JT):
    T = (JT - 2451545.0) / 36525.0
    pf = polys(T, 5)
    pl_ = polys(T, 2)
    df = delargs(pf)
    dl = delargs(pl_)
    zeta = ARG[0][0] + T * (ARG[0][1] + PREC)
    pargs = [p[0] + p[1] * T for p in PLANET]

    lon = main_sin(df, d['ELP01'])
    lat = main_sin(df, d['ELP02'])
    dist = main_sin(df, d['ELP03'], cosine=True)
    lon += nonplan(zeta, dl, d['ELP04']); lat += nonplan(zeta, dl, d['ELP05']); dist += nonplan(zeta, dl, d['ELP06'])
    lon += nonplan(zeta, dl, d['ELP07']) * T; lat += nonplan(zeta, dl, d['ELP08']) * T; dist += nonplan(zeta, dl, d['ELP09']) * T
    lon += plan1(pargs, dl, d['ELP10']); lat += plan1(pargs, dl, d['ELP11']); dist += plan1(pargs, dl, d['ELP12'])
    lon += plan1(pargs, dl, d['ELP13']) * T; lat += plan1(pargs, dl, d['ELP14']) * T; dist += plan1(pargs, dl, d['ELP15']) * T
    lon += plan2(pargs, dl, d['ELP16']); lat += plan2(pargs, dl, d['ELP17']); dist += plan2(pargs, dl, d['ELP18'])
    lon += plan2(pargs, dl, d['ELP19']) * T; lat += plan2(pargs, dl, d['ELP20']) * T; dist += plan2(pargs, dl, d['ELP21']) * T
    lon += nonplan(0.0, dl, d['ELP22']); lat += nonplan(0.0, dl, d['ELP23']); dist += nonplan(0.0, dl, d['ELP24'])
    lon += nonplan(0.0, dl, d['ELP25']) * T; lat += nonplan(0.0, dl, d['ELP26']) * T; dist += nonplan(0.0, dl, d['ELP27']) * T
    lon += nonplan(0.0, dl, d['ELP28']); lat += nonplan(0.0, dl, d['ELP29']); dist += nonplan(0.0, dl, d['ELP30'])
    lon += nonplan(0.0, dl, d['ELP31']); lat += nonplan(0.0, dl, d['ELP32']); dist += nonplan(0.0, dl, d['ELP33'])
    lon += nonplan(0.0, dl, d['ELP34']) * T * T; lat += nonplan(0.0, dl, d['ELP35']) * T * T
    dist += nonplan(0.0, dl, d['ELP36']) * T * T
    lon += pf[0]
    return lon * RAD, lat * RAD, dist, (main_sin(df, d['ELP01']),)


if __name__ == '__main__':
    for jd in [2448724.5, 2451545.0]:
        lon, lat, dist, (m1,) = spherical(jd)
        print(f'JD {jd}:  参考实现 λ={np.degrees(lon) % 360:.9f}  '
              f'β={np.degrees(lat):.9f}  Δ={dist:.4f} km   (ELP01 主问题和 = {m1:.4f}")')
