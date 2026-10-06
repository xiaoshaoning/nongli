#!/usr/bin/env python3
"""对 **NASA《五千年日月食目录》** 验收日食/月食预测。

真值目录见 `third_party/eclipse-catalog/`（Espenak & Meeus 的 NASA GSFC
目录，全世界通用的标准参考）。本脚本：

1. 解析那份目录，**并用每页自己的统计表校验解析结果**——解析出的各类计数
   必须与统计表逐项相同，否则直接报错（这一步抓到过两个解析坑，见该目录 README）。
2. 跑 `crates/ephemeris/examples/eclipse_check.rs`，取本实现对 1901–2100
   每一个朔/望的判断。
3. 按**日期**配对，比较类型。

判据是"日期与类型都对得上"。类型里日食只到**中心食/偏食**这一层：
全食/环食/全环食的区分取决于观测者站在食带的哪一段，目录给的是"最大食点处"
的分类，本实现不追那一层（见 `crates/ephemeris/src/eclipse.rs` 的说明）。

用法:  python tools/eclipsecheck.py
"""
import datetime
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, os.pardir)
sys.path.insert(0, os.path.join(ROOT, 'probe'))
import dtreference  # noqa: E402

CAT = os.path.join(ROOT, 'third_party', 'eclipse-catalog')
FILES = ['SE1901-2000', 'SE2001-2100', 'LE1901-2000', 'LE2001-2100']
# 目录的类型字母 → 本实现的类型字母
SOLAR_MAP = {'T': 'C', 'A': 'C', 'H': 'C', 'P': 'P'}
LUNAR_KIND = {'N': 'N', 'P': 'P', 'T': 'T'}


def strip(seg):
    t = re.sub(r'<[^>]+>', ' ', seg)
    t = re.sub(r'&nbsp;?', ' ', t)
    return re.sub(r'\s+', ' ', t).strip()


def parse(fn):
    """→ (solar, lunar)，各为 {日期: 类型字母}。"""
    s = open(os.path.join(CAT, fn + '.html'), encoding='latin-1').read()
    if fn.startswith('SE'):
        # 日食：类型是行内 hh:mm:ss 之后的第 5 个词，可能是两字符（Pe / T- / H3）
        out = {}
        for m in re.finditer(r'5MCSEmap/\d{4}-\d{4}/(\d{4}-\d{2}-\d{2})\.gif"[^>]*>\s*\d+</a>', s):
            f = strip(s[m.end():m.end() + 300]).split()
            i = next((k for k, t in enumerate(f) if re.fullmatch(r'\d{2}:\d{2}:\d{2}', t)), None)
            if i is not None and len(f) > i + 4 and re.fullmatch(r'[A-Za-z][A-Za-z0-9+\-]*', f[i + 4]):
                out[m.group(1)] = SOLAR_MAP.get(f[i + 4][0].upper())
        return out, {}
    # 月食：类型就在文件名里，而且用的是 5MCLEmap
    out = {}
    for a, b, c, k in re.findall(r'5MCLEmap/\d{4}-\d{4}/LE(\d{4})-(\d{2})-(\d{2})([A-Za-z])\.gif', s):
        out[f'{a}-{b}-{c}'] = LUNAR_KIND[k.upper()]
    return {}, out


def stats(fn):
    s = open(os.path.join(CAT, fn + '.html'), encoding='latin-1').read()
    return {m.group(1): int(m.group(2)) for m in re.finditer(
        r'<td>\s*(Total|Annular|Partial|Hybrid|Penumbral)\s*</td>\s*<td>\s*[A-Z]?\s*</td>\s*<td>\s*(\d+)\s*</td>',
        s)}


def main():
    if not os.path.isdir(CAT):
        print(f'缺少 {CAT} —— 跳过。')
        return 0

    # ---- 1. 解析目录，并用页面自己的统计表校验 ----
    cat_s, cat_l = {}, {}
    for f in FILES:
        s, l = parse(f)
        cat_s.update(s)
        cat_l.update(l)
        want = stats(f)
        got = {}
        for k in (s or l).values():
            got[k] = got.get(k, 0) + 1
        # 目录用 T/A/H/P/N，本实现用 C/P/N/T，换算后再比
        rev = {'T': 'T', 'A': 'A', 'H': 'H', 'P': 'P', 'N': 'N'}
        exp = {}
        for name, n in want.items():
            exp[{'Total': 'T', 'Annular': 'A', 'Hybrid': 'H',
                 'Partial': 'P', 'Penumbral': 'N'}[name]] = n
        if f.startswith('SE'):
            # 日食侧只映射到 C/P，统计表要合并 T+A+H
            merged = {'C': exp.get('T', 0) + exp.get('A', 0) + exp.get('H', 0),
                      'P': exp.get('P', 0)}
            ok = got == merged
        else:
            ok = got == exp
        print(f'  {f}: 解析 {len(s) + len(l):4} 条  类型 {got}  '
              f'{"与统计表一致" if ok else "与统计表不符 " + str(exp)}')
        if not ok:
            print('  解析与页面统计不符，停止。')
            return 1

    print(f'\n目录：日食 {len(cat_s)} 次、月食 {len(cat_l)} 次（1901–2100）')

    # ---- 2. 本实现的判断 ----
    out = subprocess.run(
        ['cargo', 'run', '--release', '--quiet', '-p', 'ephemeris',
         '--example', 'eclipse_check'],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    # Rust 侧已给出世界时日期（用它自己的 ymd_from_jdn），这里不再抄历法换算
    our_s, our_l = {}, {}
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        tgt = our_s if p[0] == 'SOLAR' else our_l
        tgt.setdefault(p[1], []).append(p[2])

    # ---- 3. 按日期配对 ----
    def shift(date, days):
        y, m, d = (int(x) for x in date.split('-'))
        return (datetime.date(y, m, d) + datetime.timedelta(days=days)).isoformat()

    def cmp(cat, our, name):
        hit = miss = extra = wrong = 0
        bad = []
        used = set()
        for date, kind in sorted(cat.items()):
            # 允许 ±1 天：目录给的是 UT 日期，且最大食时刻的定义略有差别
            cand = [(d, our[d]) for d in (date, shift(date, 1), shift(date, -1))
                    if d in our and d not in used]
            if not cand:
                miss += 1
                bad.append(f'漏 {date} 应为 {kind}')
                continue
            d, kinds = cand[0]
            used.add(d)
            hit += 1
            if kind not in kinds:
                wrong += 1
                bad.append(f'类型不符 {date}: 目录 {kind}，本实现 {"/".join(kinds)}')
        for d in our:
            if d not in used:
                extra += 1
                bad.append(f'多 {d} 本实现 {"/".join(our[d])}')
        print(f'\n{name}：目录 {len(cat)}  配对 {hit}  漏 {miss}  多 {extra}  类型不符 {wrong}')
        for b in bad[:12]:
            print('   ', b)
        if len(bad) > 12:
            print(f'    …还有 {len(bad) - 12} 条')
        return len(cat), hit, miss, extra, wrong

    cmp(cat_s, our_s, '日食（中心食 C / 偏食 P）')
    cmp(cat_l, our_l, '月食（N 半影 / P 偏 / T 全）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
