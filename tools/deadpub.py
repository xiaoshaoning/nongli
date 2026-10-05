#!/usr/bin/env python3
"""找出两类**接口面太宽**的公开项。

`rustc` 只对私有项报 dead_code；`pub` 项永远不会被警告，于是很容易留下
"没人用的公开接口"。它白白加宽 crate 的接口面，正是 Ousterhout 说的
**Overexposure**。这个检查已经抓到过多次真货。

两类：

1. **死 API**——某 `pub` 项在所有非测试代码里出现次数 ≤ 1（只有定义本身）。
2. **只在本 crate 内用的 `pub`**——出现了很多次，但一次都不在"该 crate 的库目标"
   之外。这类同样该收成 `pub(crate)`，只是引用次数多，第 1 条抓不到。

关于第 2 条的边界：`src/bin/**`、`examples/**`、`tests/**` 在 Rust 里各是**独立的
crate**，它们用到的 `pub` 项确实必须是 `pub`，算正当消费者。所以判断"库目标之外"
时要把它们算进去，只把 `crates/X/src/**`（除 `bin/`）里、只有 crate X 自己的
文件引用到的项判为过宽。

`#[cfg(test)]` 之后的部分不计——测试引用不算消费者。

用法:  python tools/deadpub.py     （有问题时退出码 1）
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
ALL = sorted(ROOT.glob('crates/**/*.rs'))
ITEM = re.compile(r'^pub (?:fn|struct|enum|const|static|trait|type) ([A-Za-z_0-9]+)', re.M)

# 第 2 条只查这几种。`struct`/`enum`/`trait` 常出现在**公开函数的签名**里
# （`pub fn apparent() -> ApparentEcliptic`），名字检查看不出这一点，
# 会把它们全误报。`fn`/`const`/`static`/`type` 的名字几乎只会因为
# "有人直接调用或读取"才出现在别的文件里。
WIDE = re.compile(r'^pub (?:fn|const|static|type) ([A-Za-z_0-9]+)', re.M)


def strip_tests(text):
    i = text.find('#[cfg(test)]')
    return text[:i] if i >= 0 else text


def lib_crate(path):
    """若 path 属于某个 crate 的**库目标**，返回 crate 名；否则 None。"""
    parts = path.relative_to(ROOT).parts        # crates/<name>/src/...
    if len(parts) < 4 or parts[0] != 'crates' or parts[2] != 'src':
        return None
    if parts[3] in ('bin', 'main.rs'):          # bin 是独立 crate，不是库目标
        return None
    return parts[1]


def main():
    texts = {f: strip_tests(f.read_text(encoding='utf-8')) for f in ALL}
    homes, wideable = {}, set()
    for f in ALL:
        if lib_crate(f) is None:
            continue
        for m in ITEM.finditer(texts[f]):
            homes.setdefault(m.group(1), f)
        for m in WIDE.finditer(texts[f]):
            wideable.add(m.group(1))

    def count(sym, files):
        pat = re.compile(r'\b' + re.escape(sym) + r'\b')
        return sum(len(pat.findall(texts[f])) for f in files)

    dead, wide = [], []
    for sym, home in homes.items():
        if count(sym, ALL) <= 1:
            dead.append((sym, home.relative_to(ROOT)))
        elif sym in wideable and count(
                sym, [f for f in ALL if lib_crate(f) != lib_crate(home)]) == 0:
            wide.append((sym, lib_crate(home), home.relative_to(ROOT)))

    if not dead and not wide:
        print(f'接口面: 干净（检查了 {len(homes)} 个公开项）')
        return 0
    if dead:
        print(f'死 API（{len(dead)}）:')
        for sym, home in sorted(dead):
            print(f'  {sym:34} {home}')
    if wide:
        print(f'只在本 crate 内用的 pub（{len(wide)}）——收成 pub(crate):')
        for sym, krate, home in sorted(wide):
            print(f'  {sym:34} [{krate}] {home}')
    print('\n清掉，或等真有消费者再加。别拿"以后要用"或"已经验证过"当理由。')
    return 1


if __name__ == '__main__':
    sys.exit(main())
