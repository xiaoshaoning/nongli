# third_party/elp2000-82b

## 为什么有这份东西

这是 **ELP2000-82B 月球理论的完整系数**，用来替换本仓库原先的 Meeus 60 项截断级数。

| | 现状（Meeus 第 47 章） | ELP2000-82B |
|---|---|---|
| 系数项数 | 60 | **3402** |
| 绝对精度（对 ELP/MPP02，SOFA 文档） | RMS 2.9″、最坏 18.3″ | 作者标注截断误差 **0.001″ / 0.001 km** |

这是第 7 步（`docs/plan.md`）唯一的阻塞项，2026-10-05 网络恢复后取得。

## 出处

* 仓库：https://github.com/vsr83/ELP2000-82B
* 许可：**MIT**（见 `LICENSE`）——与 MIT 的本仓库兼容
* 算法依据：Chapront-Touzé, Chapront & Francou, *Lunar solution ELP version
  ELP 2000-82B* (1985)，并采用 Chapront, Chapront-Touzé & Francou (2002, A&A)
  基于 LLR 重新确定的轨道参数（见 `Elp2000-82b.js` 顶部注释）

## 文件

| 文件 | 用途 |
|---|---|
| `ELP2000-82b.json` | 全部 36 张系数表（ELP01–ELP36），共 3402 项 |
| `Elp2000-82b.js` | 参考实现（JS），本仓库移植时照此写 Rust |
| `create_json_elp2000.js` | 说明 JSON 是如何从原始 ELP 数据文件生成的 |
| `horizons_results_monthly_1900-2100.js` | **JPL Horizons 真值**，1900–2100 逐月，`[JD_TDB, x, y, z]` km |
| `elp2000ref_results_monthly_1900-2100.js` | 该仓库自己算出的结果，可用于交叉核对 |

`horizons_results_monthly_1900-2100.js` 是独立于本实现的**真值**——用它来验收
移植是否正确，比自己跟自己比有意义得多。

## 表的编码

* **ELP01–ELP03**（主问题）：`{i1..i4, A, B1..B6}`，项为
  `A·sin(i1·D + i2·l' + i3·l + i4·F)`（ELP03 是 `cos`，算距离）。
  `B1..B6` 是 DE200/LE200 拟合改正要用的偏导系数。
* **ELP04–ELP09、ELP22–ELP36**（非行星扰动 / 潮汐 / 月体形状 / 相对论）：
  `{i1..i5, phi, A, B}`，项为 `A·sin(i1·ζ + i2·D + i3·l' + i4·l + i5·F + φ)`。
* **ELP10–ELP21**（行星扰动两张表）：`{i1..i11, phi, A, B}`，前 8 个乘子是行星平黄经，
  后 3 个是 `D, l, F`。

角度一律是**角秒**（`rad = deg × π/648000`），`φ` 是**度**。

## 参考实现有一处历元错误（2026-10-05 实测）

`Elp2000-82b.js` 的函数叫 `elp2000SphericalOfDate`，但它给的实际是**混合口径**：

* **黄经以 J2000 平春分点起算**（ELP 的 W1 是恒星黄经，不含岁差）——
  与 DE421 按黄道 J2000 读差 0.19″，按当日黄道读差 388″（1992 年）；
* **黄纬已经是当日黄道**——与 DE421 按当日黄道读差 0.000″，按 J2000 读差到 46″。

它自带的 `horizons_results_monthly_1900-2100.js` 用的是同一个历元，所以两边
"自比自"地通过了（那份数据也确实是黄道 J2000：x 分量与 ICRS 逐位相同）。

本仓库的用法：**只补黄经岁差**，黄纬不动。见 `docs/plan.md` 第 7 步。
逐行复刻那份 JS 的 Python 脚本在 `probe/elp_ref.py`。

## 注意

这些文件**不参与编译**，只是系数来源与验收依据。Rust 侧的
`crates/ephemeris/src/moon_tables.rs` 由 `tools/gen_elp2000.py` 生成。
