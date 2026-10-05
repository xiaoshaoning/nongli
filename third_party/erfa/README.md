# third_party/erfa

## 为什么有这份代码

`nut00a.c` 是 **IAU 2000A 章动模型**（MHB2000 日月章动 + 行星章动）的参考实现，
含 678 条日月项与 687 条行星项。本仓库原先只用 4 项缩写式章动（误差 0.34″），
而把这份权威表**转录成 Rust** 需要原始数据，故留下出处。

`tools/gen_nutation.py` 读它、生成 `crates/ephemeris/src/nutation_tables.rs`。
**它不参与编译**，只是系数来源，便于离线复现。

## 出处

* 项目：ERFA（Essential Routines for Fundamental Astronomy）
* 仓库：https://github.com/liberfa/erfa
* 文件：`src/nut00a.c`
* 下载于 2026-10-05，取自 `master`

## 许可

ERFA 采用 BSD 3-clause（见同目录 `LICENSE`）。ERFA 由 IAU SOFA 派生，
但不受 SOFA "只能原样分发" 的限制。分发条件第 1 条要求保留版权声明——
`LICENSE` 与 `nut00a.c` 文件头都原样保留。

本仓库自身的代码是 MIT；两者不冲突。
