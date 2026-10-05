#!/usr/bin/env python3
"""找出**没有消费者的公开项**（死 API）。

`rustc` 只对私有项报 dead_code；`pub` 项永远不会被警告，于是很容易留下
"没人用的公开接口"——它白白加宽 crate 的接口面，是 Ousterhout 说的
Overexposure。这个检查已经抓到过两次真货。

判据：某 `pub` 项在**所有非测试代码**里的出现次数 ≤ 1（只有定义本身）。
`#[cfg(test)]` 之后的部分不计——测试引用不算消费者。

用法:  python tools/deadpub.py
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FILES = sorted(ROOT.glob('crates/*/src/**/*.rs'))


def strip_tests(text):
    i = text.find('#[cfg(test)]')
    return text[:i] if i >= 0 else text


def main():
    texts = {f: strip_tests(f.read_text(encoding='utf-8')) for f in FILES}
    items = {}
    for f, txt in texts.items():
        for m in re.finditer(r'^pub (?:fn|struct|enum|const|static) ([A-Za-z_0-9]+)', txt, re.M):
            items.setdefault(m.group(1), f)

    dead = []
    for sym, home in items.items():
        pat = re.compile(r'\b' + re.escape(sym) + r'\b')
        n = sum(len(pat.findall(txt)) for txt in texts.values())
        if n <= 1:                      # 只有定义那一处
            dead.append((sym, home.relative_to(ROOT)))

    if not dead:
        print(f'死 API: 无（检查了 {len(items)} 个公开项）')
        return 0
    print(f'死 API（{len(dead)} / {len(items)} 个公开项）:')
    for sym, home in sorted(dead):
        print(f'  {sym:34} {home}')
    print('\n清掉，或等真有消费者再加。别拿"以后要用"或"已经验证过"当理由。')
    return 1


if __name__ == '__main__':
    sys.exit(main())
