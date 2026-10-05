#!/usr/bin/env python3
"""查文档与代码的脱节。两类：

1. **死引用** —— 文档里提到的仓库内路径是否存在。`git rm` 一个文件之后，
   文档里还留着它的名字，没有任何东西会告诉你。本项目已经踩过
   （`tools/gen_tables.py`、`src/tables.rs`）。

2. **数字漂移** —— 精度数字手工抄到别处之后必然过期。现在确切数字只由
   `tools/gen_accuracy.py` 写进 `docs/accuracy.md`，其余地方只引**量级**；
   本脚本把那些量级与生成值对一遍。已经漂过好几次，最严重的一次
   `crates/ephemeris/src/lib.rs` 的"精度（必读）"落后整整一个第 7 步。

`docs/plan.md` 按**历史记录**豁免（它记的是"当时量到什么"，本来就是过去的数字）。
`kernels/` 下的文件不进版本库，也豁免。

用法:  python tools/doccheck.py     （有问题时退出码 1）
"""
import glob
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# 被检查的文档
DOCS = ['README.md', 'docs/compliance.md', 'docs/plan.md', 'docs/ambiguity.md']
# `docs/plan.md` 是**执行日志**，两条检查都豁免它：
# 它记的是"当时做了什么"，必然要提到后来被删掉的文件（"我们删了 X"本身是有效的
# 历史），数字也是"当时量到什么"。面向使用者的文档（README / compliance /
# accuracy）才要求路径存活、数字与生成值一致。
HISTORICAL = {'docs/plan.md'}

# 别处引用的"量级"：(文件, 该文件里应当出现的字面量模板, tracked 键, 小数位)
MAGNITUDE_CLAIMS = [
    ('README.md', '太阳视黄经 ~{}″', 'sun_app_med_abs', 2),
    ('README.md', '月球视黄经 ~{}″', 'moon_lon_med_abs', 2),
    ('README.md', '月地距离 ~{} km', 'moon_dist_med_abs', 2),
    ('README.md', '朔 ~{} s', 'spec_new_moon_max', 1),
    ('README.md', '节气 ~{} s', 'spec_term_max', 1),
    ('docs/compliance.md', '朔 ~{} s', 'spec_new_moon_max', 1),
    ('docs/compliance.md', '节气 ~{} s', 'spec_term_max', 1),
    ('docs/compliance.md', '太阳视黄经 ~{}″', 'sun_app_med_abs', 2),
    ('docs/compliance.md', '月球 ~{}″', 'moon_lon_med_abs', 2),
]

# 仓库内路径的样子
PATH_RE = re.compile(r'`((?:crates|tools|probe|docs|third_party)/[A-Za-z0-9_./{},*-]+)`')


def read(p):
    return io.open(os.path.join(ROOT, p), encoding='utf-8').read()


def load_tracked():
    """从 docs/accuracy.md 里读机器可读的清单。"""
    s = read('docs/accuracy.md')
    m = re.search(r'<!-- tracked-claims\n(.*?)\n-->', s, re.S)
    if not m:
        raise SystemExit('docs/accuracy.md 里没有 tracked-claims 块——'
                         '请先跑 tools/gen_accuracy.py')
    out = {}
    for line in m.group(1).strip().split('\n'):
        k, v = line.split('=')
        out[k] = float(v)
    return out


def expand(pat):
    """把 `{a,b}.rs` 展开成若干具体路径。"""
    if '{' not in pat:
        return [pat]
    pre, rest = pat.split('{', 1)
    inner, post = rest.split('}', 1)
    return [pre + x.strip() + post for x in inner.split(',') if x.strip()]


def main():
    tracked = load_tracked()
    bad = []

    # ---- 1. 死引用 ----
    checked = 0
    for d in DOCS + ['docs/accuracy.md']:
        if d in HISTORICAL or not os.path.exists(os.path.join(ROOT, d)):
            continue
        for pat in sorted(set(PATH_RE.findall(read(d)))):
            for p in expand(pat):
                if p.startswith('kernels/'):      # 不进版本库
                    continue
                checked += 1
                if '*' in p:
                    # 通配（如 tools/gen_*.py）：至少要匹配到一个，否则也是死引用
                    if not glob.glob(os.path.join(ROOT, p)):
                        bad.append(f'{d}: 通配 {p} 一个文件都没匹配到')
                elif not (os.path.exists(os.path.join(ROOT, p))
                          or os.path.isdir(os.path.join(ROOT, p))):
                    bad.append(f'{d}: 引用了不存在的路径 {p}')

    # ---- 2. 量级与生成值一致 ----
    for f, template, key, dp in MAGNITUDE_CLAIMS:
        want = template.format(f'{tracked[key]:.{dp}f}')
        if want not in read(f):
            bad.append(f'{f}: 找不到按生成值渲染的 "{want}"'
                       f'（{key} = {tracked[key]}）')

    # ---- 3. 别处不该出现确切数字（除历史记录文档）----
    for d in DOCS:
        if d in HISTORICAL:
            continue
        s = read(d)
        for key, v in tracked.items():
            # 生成值本身在别处出现 = 手抄了确切数字
            if f'{v:g}' in s:
                bad.append(f'{d}: 出现了确切数字 {v:g}（{key}）——'
                           f'只写量级并链到 docs/accuracy.md')

    if bad:
        print(f'文档与代码脱节 {len(bad)} 处（检查了 {checked} 个路径引用）：')
        for b in sorted(set(bad)):
            print('  ✗', b)
        return 1
    print(f'文档检查: 干净（{checked} 个路径引用；'
          f'{len(MAGNITUDE_CLAIMS)} 条量级与生成值一致）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
