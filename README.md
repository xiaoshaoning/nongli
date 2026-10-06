# nongli —— 中国农历编算 (GB/T 33661-2017)

按 **GB/T 33661-2017《农历的编算和颁行》** 用 Rust 实现的农历计算库与命令行工具。
给定普通日历（公历）的任意日期时刻，算出对应的农历表示：农历年/月/日名称、干支纪年、
生肖纪年、干支纪日、二十四节气、闰月。

底下的 `ephemeris` crate 是通用的：太阳与月球的地心位置、参考系归算、站心地平坐标、
月相与月球天平动、**日食与月食的发生与类型**。它不知道任何农历概念。

```
$ nongli 2026-02-17
公历 2026-02-17 (星期二)
农历丙午年正月初一   [丙午马年 正月初一]
干支  丙午年 壬戌日
下一节气  雨水 (1 天后)

$ nongli 2033-12-22
公历 2033-12-22 (星期四)
农历癸丑年闰十一月初一   [癸丑牛年 闰十一月初一]
干支  癸丑年 丁未日
```

## 快速开始

```rust
use nongli::{Calendar, DateTime, LunarDate};

let d = LunarDate::from_date(2026, 2, 17, Calendar::Gregorian);
assert_eq!((d.year, d.month, d.day), (2026, 1, 1));
assert_eq!(d.to_string(), "农历丙午年正月初一");
assert_eq!(d.ganzhi_year(), "丙午");
assert_eq!(d.zodiac(), "马");
assert_eq!(d.ganzhi_day(), "壬戌");

// 带时刻 (北京时间)：农历日以北京时间为界
let nye = DateTime::new(2026, 2, 17, 23, 59, 59.0);
let new = DateTime::new(2026, 2, 18, 0, 0, 0.0);
assert_eq!(LunarDate::from_datetime(nye, Calendar::Gregorian).day, 1);
assert_eq!(LunarDate::from_datetime(new, Calendar::Gregorian).day, 2);

// UTC 输入：北京时间 = UTC + 8h
let d = LunarDate::from_utc(DateTime::new(2026, 2, 17, 7, 30, 0.0));
assert_eq!(d.day, 1); // 北京 15:30
let d = LunarDate::from_utc(DateTime::new(2026, 2, 17, 17, 0, 0.0));
assert_eq!(d.day, 2); // 北京次日 01:00

// 节气
assert_eq!(LunarDate::from_date(2024, 4, 4, Calendar::Gregorian).solar_term(), Some("清明"));
```

命令行：

```
nongli <日期>[T时刻] ...     换算农历（默认北京时间）
nongli --year  <年>          打印全年公历/农历对照 + 节气
nongli --terms <年>          打印全年二十四节气
nongli --sui   <年>          打印一个"岁"的各月与置闰判定
     --utc                   输入按 UTC 解释
     --julian                输入按儒略历解释
```

## 实现的编排规则

| 标准条款 | 内容 | 代码位置 |
|---|---|---|
| 4.1 | 以北京时间为标准时间（东经 120° 标准时） | `jd::BEIJING_OFFSET_DAYS` |
| 4.2 | 朔日为农历月的第一个农历日 | `calendar::new_moon_on_or_before` |
| 4.3 | 包含节气冬至的农历月为十一月 | `calendar::winter_solstice_jdn` |
| 4.4 | 某十一月到下一个十一月（不含）之间有 13 个月时置闰，取最先出现且不含中气的月为闰月 | `calendar::lunar_from_jdn` |
| 4.5 | 十一月之后第 2 个（不计闰月）农历月为年的起始月 | `calendar::month_number` |
| 5.1 | 太阳、月球位置按 IERS 规范模型计算 | `astro`（VSOP87D + ELP2000-82B） |
| 6.1.1 | 干支纪年，1984 农历年为甲子年 | `LunarDate::ganzhi_year` |
| 6.1.2 | 生肖纪年，同上参考年为鼠年 | `LunarDate::zodiac` |
| 6.2 | 数序纪月（正月…十二月，闰月加"闰"字） | `LunarDate::month_name` |
| 6.3.1 | 数序纪日（初一…初十、十一…二十、廿一…廿九、三十） | `names::day_name` |
| 6.3.2 | 干支纪日，1949-10-01 的农历日为甲子日 | `LunarDate::ganzhi_day` |
| 6.4 | 农历日期表示方法 | `impl Display for LunarDate` |
| 附录 A | 二十四节气名称与相应太阳地心视黄经 | `names::TERM_NAMES` |

## 计算模型

系数表全部由 `tools/gen_*.py` 生成到 `crates/ephemeris/src/*_tables.rs`，
**生成器的文档字符串里写了来源、截断阈值与实测误差**。三套级数都用**原始系数**，
不用教科书上的转抄版（转抄版实测在 J2000 处就有 −0.24″ 的系统差）。

| 量 | 模型 | 来源 | 备注 |
|---|---|---|---|
| 行星 | **VSOP87D** 地球级数 | CDS VI/81 `VSOP87D.ear` | 1080/348/997 项，按实测阈值截断 |
| 月球 | **ELP2000-82B** | `vsr83/ELP2000-82B`（MIT） | 36 张表 3402 项；2002 年 LLR 重定参数 |
| 章动 | **IAU 2000A 完整模型** | ERFA `nut00a.c`（BSD） | 678 日月项 + 687 行星项 |
| 岁差 | **IAU 2006 (P03)** | 对 `erfa.pfw06` 拟合 FW 角 | ±200 世纪内 2.3e-8″ |
| 黄赤交角 | IAU 2006 | 同上（EPSA 角） | 对 `erfa.obl06` 差 1.8e-6″ |
| 恒星时 | ERA + IAU 2006 GMST | — | 分点差按 ERFA `ee00.c` 的定义 |
| ΔT | 1973–2025 IERS 实测值；其余 Espenak & Meeus (2006) 分段多项式 | 内建小表 | 表端做增量衔接保证连续 |

两处口径要留意（都写进了接口注释）：

* **参考面**：太阳与月球都换算到 **IAU 2006 的当日平黄道分点**。VSOP87D 与
  ELP2000-82B 自带的是 1980 年代那套（差 0.003″/年），要显式改正。
* **观测效应**：月球输出是**视**位置（几何 + 光行时 + 周年光行差）；太阳的
  视黄经与几何黄经是两个分开的入口（"朔"用后者）。

几何：朔（日月**视**黄经相等）与节气（太阳视黄经 = 15° 整数倍）都用角度差上的
牛顿迭代求解。

## 适用范围：在多大范围内没有歧义

**完整分析（含全部实测数据与复现步骤）见 [`docs/ambiguity.md`](docs/ambiguity.md)。**
摘要：

**编排规则本身没有歧义**——GB/T 的规则是确定性的。歧义只来自一个地方：
一个**朔**或**中气**若恰好落在北京时间午夜前后，它的“哪一天”就取决于 ΔT 与星历
模型的精度。而朔决定月界、中气决定闰月。

实测（`probe/ambiguity.py`，2026±4000 共 194980 个朔与中气）：

| 距午夜余量 | 1% 的年份 | 5% | 50% |
|---|---|---|---|
| 该年最“危险”事件的余量 | ≤ 14.5 s | ≤ 84.1 s | ≤ 1163 s |

余量近似均匀分布，故“某年存在临界事件”的概率 ≈ U / 1770（U = 不确定度秒数）：

| 不确定度 U | 1 s | 10 s | 60 s | 300 s | 1 h |
|---|---|---|---|---|---|
| 该年出现临界事件的概率 | 0.1% | 0.7% | 3.5% | 16% | 66% |

**但“哪一天”比“哪一年/闰几月”脆弱得多。** 把 ΔT 整体平移（`probe/dtshift.py`）：

| ΔT 整体偏移 | 月数/闰月索引/正月索引发生变化的年份占比（全区间） |
|---|---|
| ≤ 60 s | **0.0%** |
| 600 s (10 min) | ≤ 0.7% |
| 3600 s (1 h) | 1.3% – 4.2% |
| 21600 s (6 h) | 12% – 16% |

因为整体平移同时移动了朔和中气的日界，而闰月判定只比较两者的**相对**关系。
换句话说：**即便在 ±4000 年，“哪一年、闰几月”仍然基本可信，不可信的只是“初几”。**

### 结论

| 区间 | ΔT 不确定度（量级） | 日边界 | 年/月/闰月 |
|---|---|---|---|
| 1973–2025 | < 0.1 s（IERS 实测，已内建） | 唯一 | 唯一 |
| 1620–2025 | ≲ 2 s | 唯一 | 唯一 |
| 2026 ± 50 | ≲ 5 s | 基本唯一 | 唯一 |
| 2026 ± 100 | ~10 s | 百年内约 1 次可能 ±1 天 | 唯一 |
| 2026 ± 500 | ~1–2 min | 约 5% 的年份有 ±1 天风险 | 几乎总对 |
| 2026 ± 1000 | ~5 min | 约 17% 的年份有 ±1 天风险 | 几乎总对 |
| 2026 ± 2000 | ~20 min | 不可唯一定日 | 个别年份闰月会变 |
| 2026 ± 4000 | ~1 h | 不可唯一定日 | 约 95% 的年份仍正确 |

过去比未来可靠：**1926–2025 的 ΔT 是实测值，不是外推**；往未来则 ΔT 本身
不可预测，2100 年各家投影已相差数十秒，4000 年后相差小时量级。

因此：

* 要“拿出去不被挑”的结果：**2026 ± 100 年**。
* 要“年月闰月总对、日偶尔差一天”：**2026 ± 500 年**。
* 更大范围仍可计算，但应视为**推算**：请把结果当作“某日 ±1 天”，
  并记住闰月也有百分之几的年份可能变。

### 另外两个与天文无关的歧义

* **公历 vs 儒略历**：标准未规定 1582 年以前“公历”的含义。库默认外推公历，
  可用 `Calendar::Julian` 切换——两者给出的农历日会差 10 天左右。
* **日的起点**：农历日以**北京时间**零点为界（3.17）。用 UTC 输入时必须先换算，
  库提供 `LunarDate::from_utc`。

## 精度

**真值参照是 JPL 星历，不是 ERFA / astropy 的模型。** 后者自身就有可观的误差：
astropy builtin 的月球对 JPL 是 20″ 量级，而本实现是 0.04″ 量级。

* 太阳视黄经 ~0.01″、月球视黄经 ~0.04″、月地距离 ~0.02 km（对 JPL 真值）
* **GB/T 33661-2017 §5.2**（朔和节气的北京时间精度应达到 1 s）：达标，
  朔 ~0.2 s、节气 ~0.6 s，余量 1.7–4.5 倍

> 上面只给**量级**。确切数字与完整表格在
> **[`docs/accuracy.md`](docs/accuracy.md)** —— 那是精度数字的**唯一来源**，
> 由 `tools/gen_accuracy.py` 生成。本节原先手抄确切数字，已经漂过好几次
> （最严重的一次落后整整一个第 7 步）。`tools/doccheck.py` 会校验这里的
> 量级与生成值一致。

符合性逐条对照见 [`docs/compliance.md`](docs/compliance.md)。

### 复现方式

```bash
python tools/truthcheck.py   # 对 JPL 真值验太阳与月球（需 kernels/de421.bsp）
python tools/spec_check.py   # §5.2 的 1 s 判据
python tools/gen_accuracy.py # 重新生成 docs/accuracy.md（精度数字的唯一来源）
python tools/doccheck.py     # 查文档与代码的脱节（死路径、数字漂移）
python tools/eclipsecheck.py # 日食/月食对 NASA《五千年日月食目录》
python tools/suncheck.py     # 太阳对 astropy（可看长跨度上两个模型怎么分离）
python tools/mooncheck.py    # 月球 λ/β/Δ 对 erfa.moon98
python tools/framecheck.py   # 章动/黄赤交角 对 ERFA 逐点比对
python tools/observercheck.py # 恒星时、站心视差 对 ERFA
python tools/moonphasecheck.py # 照亮比例/相位角 对 ERFA 向量的独立实现
python tools/daydiff.py 1526 2526  # 统计日界不符的个数
python tools/deadpub.py      # 没有消费者的公开项（rustc 不会警告 pub）
python probe/ambiguity.py    # 朔/中气距午夜的余量分布
python probe/dtshift.py      # ΔT 偏移对闰月结构的影响
python probe/dtreference.py  # ΔT 参考实现自检
python tools/gen_vsop87.py   # 重新生成太阳系数表
python tools/gen_elp2000.py  # 重新生成月球系数表
python tools/gen_nutation.py / gen_frames.py / gen_ecliptic_frame.py
```

`kernels/de421.bsp`（JPL DE421，覆盖 1900–2053）不进版本库，取法见
`tools/truthcheck.py` 头部。

### 长跨度上的模型差异

`tools/suncheck.py` 保留了一张对 astropy builtin 的表，用来观察**两个模型在长跨度上
如何分离**：太阳在 1526–2526 差 ~0.7″、3000–4000 差 ~1.4″，而那是 ERFA `epv00`
只在 1900–2100 有效造成的，不是我们的误差。

## 适用范围与限制

（歧义范围的分层结论见上文“适用范围：在多大范围内没有歧义”。）

* **近现代（2026 ± 500 年）**：可靠。已用 1949–2035 年全部春节日期、2017/2020/2023/2025
  年闰月、以及 2033 年闰十一月等公开结果核对（见 `cargo test`）。
* **2026 ± 4000 年**：程序照常给出结果，但可信度递减，原因有两个，**都与代码无关**：
  1. **ΔT 本身**在 ±4000 年有小时级的不确定度（无可观测数据可依），会直接平移
     北京时间的日界。见上文“适用范围”一节的概率表。
  2. **理论的适用区间**：VSOP87D 与 ELP2000-82B 都是对 1900–2100 前后拟合的，
     星历本身在 ±4000 年是外推。本仓库只对 **1900–2053**（JPL DE421 的覆盖）
     做过真值验证；再往外要验需 DE440/DE441。
  因此远年结果应视为推算而非定论。这一限制对所有实现（含紫金山天文台）都成立。
* 标准未规定公历的历史用法。库默认**外推公历**（1582-10-15 之前也按公历规则外推），
  可选 `Calendar::Julian` 用儒略历。
* 闰秒不处理（标准 5.2 明确“不计及编算时尚未正式发布的闰秒”）。
* **不要拿 `erfa.moon98` 或 astropy 的 builtin 当月球的真值。** 它们对 JPL 的误差
  是 **20″ 量级**，比本实现大三个数量级。`tools/mooncheck.py`、`moonphasecheck.py`、
  `daydiff.py` 里的月球残差量的是**参考的**误差（确切数字见 `docs/accuracy.md`）。

## 目录

Cargo workspace，按**知识**而非执行顺序划分：

```
crates/ephemeris/     通用星历（零依赖，不知道任何农历概念）
  src/time.rs         儒略日、公历/儒略历、ΔT。不预设时区
  src/sun.rs          VSOP87D 地球级数、太阳地心视黄经
  src/moon.rs         ELP2000-82B 完整级数 —— 理论之"缝"
  src/frames.rs       Ecliptic/Equatorial 值类型、岁差、黄赤交角、章动、黄道→赤道
  src/eclipse.rs      日食/月食的**是否发生**与类型（地心几何，不做当地见食）
  src/observer.rs     恒星时、站心视差、大气折射、地平坐标
  src/angle.rs        角度归化、角度-时刻求根
  src/vsop87_tables.rs / elp2000_tables.rs     自动生成的级数系数
  src/nutation_tables.rs / frames_tables.rs / ecliptic_frame_tables.rs
  examples/{dt_dump,sun_check,moon_check,frame_check,observer_check,truth_check,eclipse_check}.rs
                      供 tools/ 核对的输出源（truth_check 对 JPL 真值、
                      eclipse_check 对 NASA 食典）
crates/moon/          月相与月球位置 app
  src/lib.rs          相位角、照亮比例、月龄、视直径、亮边方位、出没中天
  src/render.rs       把月相画成字符图形（晨昏线椭圆，不含纹理）
  src/bin/moon.rs     命令行
crates/nongli/        农历：编排规则 + 命名 + 命令行
  src/calendar.rs     编排规则（4.2–4.5）与农历日期
  src/terms.rs        二十四节气（"冬至 = 第 0 个节气"这一约定）
  src/names.rs        干支、生肖、月名、日名、节气名
  src/bin/nongli.rs   命令行
  examples/ephemeris_dump.rs  星历导出，供校验脚本使用
tools/                系数生成与校验脚本
probe/                系数来源、ΔT 参考实现、唯一性实测
docs/                 GB/T 原文与 OCR、唯一性分析、实施计划、符合性说明
```

依赖方向单向：`nongli → ephemeris`。`nongli` 只重导出其公开签名里出现的
`Calendar` / `DateTime`，其余星历接口请直接依赖 `ephemeris`。
