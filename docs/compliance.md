# 对 GB/T 33661-2017 的符合性说明

标准里只有两条**技术性**要求需要逐条对照：**5.1**（计算模型）与 **5.2**（精度）。
其余（第 4 章的编排规则、第 6 章的表示方法、第 7 章的颁行要求）是可直接验证的
规则，见 `README.md` 与 `crates/nongli/src/`。

---

## 5.2 朔和节气的北京时间计算精度应达到 1 s

**实测达标，余量 1.7–4.5 倍。** 判据与量法：

* 真值：**JPL DE421**（覆盖 1900–2053），不是 ERFA/astropy 的模型——后者本身就
  有 20″ 量级的月球误差（见下文）。
* 量的是**时间**误差（秒），不是角误差——标准要的是前者。
* 两侧用**同一套 ΔT**，所以量到的是模型/理论的误差，正是 5.2 所指
  （ΔT 自身的不确定度是另一回事，见 `README.md` 的"适用范围"一节）。
* 工具：`python tools/spec_check.py`（数据源 `crates/nongli/examples/spec_check.rs`）。

数字见 **[`docs/accuracy.md`](accuracy.md)**（唯一来源，`tools/gen_accuracy.py` 生成）：
朔 ~0.2 s、节气 ~0.6 s，判据 1 s。

> 标准 5.2 括注"不计及编算时尚未正式发布的闰秒"。本实现不处理闰秒，
> 时间尺度用 TT ↔ UT1（经 ΔT），与闰秒无关。

---

## 5.1 太阳和月球的位置按国际地球自转和参考系服务规范规定的模型计算

"国际地球自转和参考系服务规范" = **IERS Conventions**。逐条对照如下。

### 规范实际规定了什么

查 IERS Conventions 2010 **第 5 章**（*Transformation between the International
Terrestrial Reference System and the Geocentric Celestial Reference System*，
2012-08-10 版，`https://iers-conventions.obspm.fr/`）：

* 第 5 章管的是 **ITRS ↔ GCRS 的归算**，即**参考系与地球定向**：
  `[GCRS] = Q(t)·R(t)·W(t)·[ITRS]`。它规定的是岁差-章动模型、天球中间极（CIP）、
  天球中间原点（CIO）、地球自转角（ERA）与 UT1 的定义。
* 它**不规定**行星或月球的位置理论。全章没有出现 JPL/DE 星历；唯一带
  "ephemeris" 的字样是旧名"天球历书极 / 天球历书原点"（CIP/CIO 的旧称）。
* 它明确以 IAU 2000/2006 决议为准（原文："comply with the IAU 2000/2006
  resolutions"），其中 **IAU 2000 决议 B1.6** 写着 IAU 1976 岁差模型
  与 IAU 1980 章动理论**应被取代**，**IAU 2006 决议 B1** 又把岁差换成 P03。

### 本实现对照

| IERS 规定的量 | 本实现 | 位置 |
|---|---|---|
| 章动模型 IAU 2000A（B1.6） | **IAU 2000A 完整表**，678 日月项 + 687 行星项，取自 ERFA `nut00a.c` | `nutation_tables.rs` |
| 岁差模型 P03（IAU 2006 决议 B1） | FW 角（γ̄、φ̄、ψ̄、ε_A）**对 `erfa.pfw06` 拟合**，±200 世纪内 2.3e-8″ | `frames_tables.rs` |
| 黄赤交角 IAU 2006 | ε_A，同上拟合，对 `erfa.obl06` 差 1.8e-6″ | `frames_tables.rs` |
| 框架偏置（frame bias） | 已含在 φ̄ 里（`pmat06 = fw2m(pfw06)`，不另乘） | `frames.rs` |
| 地球自转角 ERA / GMST | ERA + IAU 2006 GMST 多项式 | `observer.rs` |
| 分点差（equation of the equinoxes） | `Δψ·cos ε_A`（ERFA `ee00.c` 的定义；未含 ±0.003″ 的 `eect00`） | `observer.rs` |
| 天球参考系 ICRS | 归算矩阵与 `erfa.ecm06` **逐位相同**（已对拍） | `frames.rs` |
| 基本引数（l、l′、F、D、Ω…） | IERS 2003 表达式，取自 ERFA `fa*.c` | `gen_nutation.py` |
| TDB 与 TT | 用 **TT** 代替 TDB | 见下 |

关于最后一条：第 5 章自己写了（§5.5.2 附近）

> the largest term in the difference TDB−TT being 1.7 ms · sin l …, the resulting
> error … is periodic … **which is significantly below the required
> microarcsecond accuracy**.

本实现同样直接以 TT 作为 ELP2000-82B 的自变量（ELP 原始定义用 TDB），
差别 ≤1.7 ms，对月球相当于 0.001″。

### 规范没有规定、我们自行选择的部分

**太阳与月球的位置理论。** 本实现用 **VSOP87D**（行星）与 **ELP2000-82B**（月球）
的解析级数，而不是 JPL 数值星历。这不是 IERS 规定的量，因此不构成 5.1 的不符合；
对它的要求落在 5.2 上（上面的表已经量过）。

两者的差别有多大，对 JPL DE421 实测（`tools/truthcheck.py`，1900–2050）：

数字见 **[`docs/accuracy.md`](accuracy.md)**：太阳视黄经 ~0.01″、月球 ~0.04″。

顺带一提：IERS Conventions 第 5 章的参考文献里就有
**Chapront-Touzé & Chapront (1983), "The lunar ephemeris ELP 2000"**，
且第 5 章定义 J2000 黄道时用的 d₀ 取自 **Chapront et al. (2002)**——正是
本实现 ELP2000-82B 系数所用的那篇 LLR 重定论文。

### 一处曾经的不符合，已修

在本仓库 2026-10-05 那次改动之前，太阳的黄经是给在 **VSOP87D 自带的
"当日动力学平黄道分点"** 里，那套岁差是 IAU 1976 时代的。这**正是不符合 5.1**
的地方（IAU 2000 决议 B1.6 已明确取代它），量级为 0.003″/年、±4000 年达 130″。
现已把日月都换算到 IAU 2006 的当日平黄道，
见 `docs/plan.md` 的"消掉那 0.003″/年的斜率"与"第 7 步"两节。

---

## 其余条款

| 条款 | 要求 | 实现 |
|---|---|---|
| 4.1 | 以北京时间为标准时间 | 全部日界判定走 UT1 + 8 h（`nongli::beijing_jdn`） |
| 4.2 | 朔日为农历月的第一个农历日 | `calendar.rs` 以 `new_moon` 定月界 |
| 4.3 | 含冬至的月为十一月 | `calendar.rs` 以 `winter_solstice_jdn` 定位 |
| 4.4 | 两个十一月间有 13 个月则置闰 | 无中气规则，`calendar.rs` 有测试 |
| 4.5 | 十一月之后第 2 个月为年首 | 同上 |
| 6.1.1 | 干支纪年，参照 1984-02-02 起为甲子年 | `Calendar::ganzhi_year`，有测试 |
| 6.1.2 | 生肖纪年，与干支同年循环 | `Calendar::zodiac` |
| 6.2 | 数序纪月，闰月冠"闰"字 | `Calendar::month_name` |
| 6.3.1 | 数序纪日（初一…三十） | `Calendar::day_name` |
| 6.3.2 | 干支纪日，参照 1949-10-01 为甲子日 | `Calendar::ganzhi_day`，有测试 |
| 6.4 | "农历" + 年名 + 月名 + 日名 | `impl Display for LunarDate`：`农历丙午年正月初一` |
| 7.2 | 含公历对照与二十四节气 | `nongli --year <年>` 给全年公历/农历对照并标出节气；`--terms <年>` 单列二十四节气 |

（附录 A 把二十四个节气与太阳黄经一一对应，`nongli/src/terms.rs` 按其定义：
序号 j 的节气即太阳视黄经 = 15j (mod 360)°。）
