# third_party/vsop87

## 为什么有这份东西

`VSOP87D.ear` 是 **VSOP87 地球日心黄道坐标（当日动力学平黄道分点）的完整系数**，
2425 项（L 1080、B 348、R 997）。

本仓库原来用的是 **Meeus《Astronomical Algorithms》附录 III** 的转抄版，只有
195 项。两份的差别不是"少留了几项"，而是 Meeus 那份只保留了 L0 的 64/559 项、
且系数经过重排——实测它相对未截断的 VSOP87D：

| | 1526–2526 |
|---|---|
| Meeus 转抄版（129 个 L 项） | 中位 0.117″、最大 **0.511″** |
| 按同样项数取最大振幅 | 中位 0.061″、最大 0.249″ |
| 本仓库现在用的（818 个 L 项） | 中位 0.0006″、最大 **0.0025″** |

**在 J2000 那一刻它就已经差 −0.24″** —— 这是转抄与重排带来的系统差，
不是"截断得粗一点"。所以换成原始系数。

## 出处

* 文件：`VSOP87D.ear`，CDS  Strasbourg 天文数据中心，目录
  [`VI/81`](https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/81)
* 理论：Bretagnon P. & Francou G. (1988), *Planetary theories in rectangular and
  spherical variables: VSOP87 solutions*, A&A 202, 309
* SHA-256：`8b160c859136d467f2be7fc29efa8a9652e95516dfbde00e4c739d7ddc90ca91`
* 2442 行，版本行 `VSOP87 VERSION D4 EARTH`

CDS 的目录数据可自由再分发（引用上述文章）。

## 文件格式

每个变量段以一行表头开头，指明变量号与 τ 幂次：

```
 VSOP87 VERSION D4    EARTH     VARIABLE 1 (LBR)       *T**0    559 TERMS ...
 4310    1  0  0  0 ... 0  0.00000000000  1.75347045673  1.75347045673 0.00000000000  0.00000000000
```

变量 1 = L（日心黄经）、2 = B（日心黄纬）、3 = R（向径）。
数据行的最后五个浮点是 `S, K, A, B, C`，**本仓库只用后三个**：每项是

```
A · cos(B + C · τ)        τ = (JDE − 2451545) / 365250   儒略千年数
```

`A` 的单位：L、B 为**弧度**，R 为 **AU**。（`S`、`K` 是同一项在另一组基本引数
下的 sin/cos 幅度，本仓库用不到。）

## 生成的 Rust

`tools/gen_vsop87.py` → `crates/ephemeris/src/vsop87_tables.rs`。
截断阈值与实测误差见该脚本的文档字符串。这些文件**不参与编译**，
只是系数来源。
