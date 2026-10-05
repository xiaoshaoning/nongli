# 计划：workspace 重构 + 月相/月球位置

> 状态：**第 0 步（决策）已完成；第 1 步起未开始。**
> 已确认的决定：
> 1. 拆 **Cargo workspace**（`ephemeris` / `nongli` / `moon`）
> 2. 目标包含**站心 + 地平坐标**（从某地看）
> 3. **先做数字版**，验证正确后再画月亮，且要"画得对"
> 4. 要"天文专业级别"，愿意加更多 ELP2000 项

---

## 0. 必须先记住的精度事实

这一节是整个计划的地基。**在动手前请先读完，避免误判进度。**

### 0.1 当前的 0.088″ 不度量理论精度

`README.md` 与 `docs/ambiguity.md` 里报的“月球视黄经对 ERFA 差 0.088″（1500–2200）”，
比对对象是 `astropy.get_body('moon')`。它证明了：

* 系数没抄错
* 光行时处理与 astropy 一致
* 参考系链（黄道↔赤道）正确

它**不能**证明理论准 —— 因为两边共用同一套 Meeus 截断级数，**截断误差在作差时抵消了**。

### 0.2 一个重要发现：本实现含月球光行时，`erfa.moon98` 不含

SOFA 文档原文：

> This function is a full implementation of the algorithm published by
> Meeus except that **the light-time correction to the Moon's mean longitude
> has been omitted**.

本仓库按 Meeus 原文转录（能把例 47.a 复现到 1e-5），因此**含**该改正。实测：

| 比较 | 差 |
|---|---|
| 本实现 − `erfa.moon98` | **−0.70″（系统性）** |
| 本实现 − `get_body`（含光行时） | +0.05″ |
| `get_body` − `erfa.moon98` | −0.72″（= 光行时效应）|

**对新 crate 的两条直接后果：**

1. 理论的“缝”必须写明返回的是**几何**位置还是**含光行时**的位置。否则将来换成
   ELP/MPP02 或 DE440 时会引入 **0.7″ 的不连续**，而且很难查。
2. 任何对 `erfa.moon98` 的比对都会带一条 −0.7″ 的系统偏差，**不要误读为误差**。
   要比它，先明确把光行时关掉。

### 0.3 理论本体的绝对精度（SOFA 官方文档）

`erfa.moon98` 与 **ELP/MPP02** 在 1950–2100 的比对：

| 量 | RMS | 最坏 |
|---|---|---|
| 地心方向 | 2.9″ | 18.3″ |
| 位置 | 6.1 km | 31.7 km |
| 速度 | 36 mm/s | 172 mm/s |

**这就是本仓库当前的绝对精度上限。** 对农历日界无所谓（月球 0.55″/秒，3″ ≈ 5.5 s 时间），
对"专业级"不够——掩星、月食接触时刻要的是 0.001″ 量级。

### 0.4 归算链可以做到专业级，且全程可验证

| | 能否专业级 | 依据 |
|---|---|---|
| 归算链：岁差 IAU2006、章动 IAU2000A 级、恒星时、站心视差、光行差、大气折射、地平坐标 | **能** | ERFA 提供全部参考例程，可逐项对拍 |
| 月球理论本体 | **不能**（缺数据，见第 2 节） | 当前 2.9″ RMS 是天花板 |

### 0.5 因此有一条排序原则

> **在理论被卡在 3″ 时，把归算链磨到 0.001″ 并不会让最终结果好于 3″。**

所以归算链的目标是"正确、可验证、结构干净"，**不追位数**；瓶颈是理论本体。

---

## 1. 阻塞项（唯一）

**需要高精度月球理论数据，本机没有，网络不通。** 已核实：

* Meeus 第 47 章的 60 项已是该截断级数的**全部**，没有"更多项"可加；再往上需要
  ELP2000-82B 原始系数表（约 3.7 万项）
* `curl https://github.com` → 超时
* 本机唯一 JPL 内核 `skyfield/tests/data/de441-1969.bsp`：长弧段只有太阳（`2→299`）
  与水星（`1→199`）；**地球段 `3→399` 与月球段 `3→301` 均只有 0–4 天**
* 装有 `astropy` / `skyfield` / `jplephem` / `poliastro`，但都只提供 `builtin`（= moon98）

**不会做的事**：凭记忆写几千个系数。那会产出看起来专业、无法验证、很可能更差的东西。

### 请提供下列任一（按推荐度）

- [ ] **`de440.bsp`（约 114 MB）** ← 最推荐：一份文件同时给出"真值"与"拟合目标"，
      可让结果既专业又自包含、可验证
- [ ] ELP/MPP02 系数文件（IMCCE）
- [ ] 完整 VSOP87D 系数文件（`VSOP87D.ear` 等）
- [ ] 可用的网络出口
- [ ] 以上都拿不到 → 接受 **3″ / 6 km**，在文档中写明这是可替换的缝

---

## 2. 目标架构

```
nongli/                        workspace root
├─ Cargo.toml                  [workspace] members = ["crates/*"]
├─ crates/
│  ├─ ephemeris/               通用天文（零依赖，不知道任何农历概念）✅
│  │  └─ src/{lib,time,sun,moon,frames,angle,tables}.rs   （observer 见第 4 步）
│  ├─ nongli/                  农历：GB/T 第 4/6 章规则 + CLI ✅
│  │  └─ src/{lib,calendar,names,terms,bin/nongli.rs}
│  └─ moon/                    月相与月球位置 app   （第 5 步再建）
├─ tools/                      生成器与校验脚本（跨 crate）
├─ probe/
├─ docs/
└─ README.md
```

依赖方向（必须单向，不得回指）：

```
nongli ──> ephemeris
moon   ──> ephemeris          （若以后要显示农历日期，再加 moon ──> nongli）
```

### 划界原则：按**知识**划分，不按文件

`ephemeris` **不得知道**任何农历概念。第 1 步之前 `src/astro.rs` 里混着下面这些
农历专用的符号，现已按下表归位：

| 第 1 步之前的符号 | 为什么是农历专用 | 去向 |
|---|---|---|
| ✅ `WINTER_SOLSTICE_LONGITUDE` | “冬至”是农历概念（附录 A 第 1 行） | `nongli/src/terms.rs` |
| ✅ `term_index` | “冬至是第 0 个节气”是农历约定 | `nongli/src/terms.rs` |
| ✅ `solar_term_jde` | “节气”是农历概念 | `nongli/src/terms.rs` |
| ✅ `solar_term_index_near` | 同上 | `nongli/src/terms.rs` |
| ✅ `DEG_PER_TERM` | 15° 只用在节气里 | `nongli/src/terms.rs` |
| ✅ `MEAN_LONGITUDE_AT_J2000_DEG` / `MEAN_MOTION_DEG_PER_DAY` | 是**太阳本身**的量，不是农历 | **留在 `ephemeris/src/sun.rs`** |
| ✅ `solve_sun_longitude` | “太阳黄经 = 给定值”**本身是通用的** | 重命名为 `sun_longitude_at`，留在 `ephemeris/src/sun.rs` |
| ✅ `new_moon_index_near` | 月相 app 也要用（算月龄） | **留在 `ephemeris/src/moon.rs`** |

时区硬编码也已消除（第 1 步）：`BEIJING_OFFSET_DAYS` 从 `ephemeris` 移除，
换成 `tt_to_jdn(jd_tt, utc_offset_hours)` / `jdn_in_offset(jd_ut1, offset_hours)`；
`nongli::BEIJING_OFFSET_HOURS = 8.0` 是**农历的政策**，留在 `nongli`。

模块划分与接口设计见下面两节。

### 模块划分：按**知识**，不按教科书章节，也不按执行顺序

| 知识 | 模块 | 现有来源 |
|---|---|---|
| 历法、儒略日、ΔT、`Instant` | `time.rs` | 现 `jd.rs` |
| 太阳理论（VSOP87D、视黄经、光行差） | `sun.rs` | 现 `astro.rs` 太阳部分 |
| **月球理论（ELP 截断、λ/β/Δ、光行时）** | `moon.rs` | 现 `astro.rs` 月球部分 |
| 参考系归算（岁差 IAU2006、章动、黄赤交角、黄道↔赤道） | `frames.rs` | 新增 |
| 观测者（大地坐标、恒星时、视差、折射、时角→地平） | `observer.rs` | 新增 |
| 级数系数（生成物） | `tables.rs` | 现同名 |

比早先草案的 8 个文件少，理由（“better together or apart”）：

* `jd` 与 `time_scales` 是同一份知识 → 合并为 `time.rs`（原名 `jd` 太窄：里面还有 ΔT 与 `DateTime`）
* `precession` + `nutation` + `obliquity` 三者互锁（没有单独使用的场景）→ 一个 `frames.rs`
* `sidereal` 只被观测者用 → 并入 `observer.rs`

### 公开接口：设计两次

三个消费者真正要的东西：

* `nongli`：太阳视黄经（及其反解）、朔时刻
* `moon`：月球的位置与形状
* 不变量：两个 crate 都不应该学会“归算配方”（光行时→章动→岁差→视差→折射的**顺序与取舍**）

**设计 A（分层库，即早先草案暗示的）**：每个阶段一个 `pub fn`，分布在 6 个模块。

> 接口 ≈ 实现。配方泄漏给每个调用方；~15 个入口点；调用方必须知道 TT/UT1 之分、
> 帧的定义、哪些改正可选。**否决**（Overexposure + Information Leakage）。

**设计 B（一次深调用）**：`moon::observe(instant, observer) -> Observation`（一个 god struct）。

> 接口极小，配方隐藏。但拿不到中间量（地心黄道位置本身就是合理的科学输出），
> 且 struct 混了多个坐标系。**否决**（不够“somewhat general-purpose”）。

**设计 C（值类型 + 分阶段转换 + 一步到位）—— 选它**

```rust
// crates/ephemeris/src/lib.rs

/// 一个时刻。TT 与 UT1 都由它携带，避免把两种时间尺度搞混。
///
/// **第 2 步才落**（第 1 步只做到参数化时区）：它此刻没有生产者也没有消费者，
/// 提前建出来是脚手架。见第 1 步末尾"实施中偏离计划的两处"。
pub struct Instant { /* tt_jd, ut1_jd 私有 */ }
impl Instant {
    /// 由 TT 儒略日构造（UT1 内部由 ΔT 推出）。
    pub fn from_tt(tt_jd: f64) -> Self;
    /// 由某历法下的地方民用时构造。`utc_offset_hours` 如 +8.0。
    pub fn from_civil(dt: DateTime, cal: Calendar, utc_offset_hours: f64) -> Self;
    pub fn tt_jd(&self) -> f64;
}

pub struct Observer { pub lat_deg: f64, pub lon_deg: f64, pub height_m: f64 }

/// 地心黄道坐标，参考**当日平**分点。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Ecliptic { pub lon_deg: f64, pub lat_deg: f64, pub distance_km: f64 }
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Equatorial { pub ra_deg: f64, pub dec_deg: f64, pub distance_km: f64 }
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Horizontal { pub azimuth_deg: f64, pub altitude_deg: f64 }

// ---- 理论的“缝” ----
/// 月球地心**黄道几何**位置（含光行时，不含章动；平黄道分点）。
/// 这里就是将来换成 ELP/MPP02 或 DE440 时唯一要改的地方。
pub fn moon_geocentric(t: Instant) -> Ecliptic;
/// 太阳地心视黄经（真分点，含光行差与章动）。农历只需要这个。
pub fn sun_apparent_longitude(t: Instant) -> f64;
/// 解 `sun_apparent_longitude == target`，`near` 需落在该解 ±7 天内。
pub fn sun_longitude_at(target_deg: f64, near: Instant) -> Instant;
/// 第 k 次朔。
pub fn new_moon(k: i64) -> Instant;

// ---- 归算（配方藏在这里）----
impl Ecliptic {
    /// 平黄道 → **视**黄道（真分点）：加章动。
    pub fn apparent(self, t: Instant) -> Ecliptic;
    /// 黄道 → 视赤道。
    pub fn equatorial(self, t: Instant) -> Equatorial;
}
impl Equatorial {
    /// 视赤道 → 站心地平。`refraction` 为真时计入大气折射。
    pub fn horizontal(self, t: Instant, o: Observer, refraction: bool) -> Horizontal;
}

// ---- 常见情况一步到位 ----
/// 从某地看到的月球：地平坐标 + 赤道坐标 + 黄道坐标 + 相位。
pub fn observe_moon(t: Instant, o: Observer, refraction: bool) -> MoonObservation;
```

接口大小：**4 个函数 + 3 个值与 3 个方法 + 2 个 struct**。
比设计 A 的 ~15 个入口小，但保留了取中间量的能力。

#### 为什么 `Instant` 要把 TT 与 UT1 藏起来

“把 TT 当成 UT1 传进去”是这类库最常见的 bug，而月球 0.55″/秒 的角速度会让它
立刻变成几十角秒的错。把两者绑在一个不透明值里，并让 `from_civil` 内部处理 ΔT，
是这个错误**定义掉不存在**（principle 7），而不是多加一个参数靠调用方自觉。

#### 缝的契约（**必须写进接口注释**）

由 §0.2 的发现，`moon_geocentric` 必须写明：

| 项 | 取值 |
|---|---|
| 参考点 | 地心 |
| 参考面/点 | **平**黄道、**平**分点（*不是*真分点）|
| 光行时 | **已含**（与 Meeus 一致，与 `erfa.moon98` 不同）|
| 章动 | **不含**（由 `Ecliptic::apparent` 加）|
| 单位 | 度、度、km |

> 若不写明这五项，将来换理论时会引入 0.7″ 的不连续，而且很难查。

> 暂不引入 `trait`：当前只有一个实现，trait 是空转。等真有第二个实现再抽。

---

## 3. 验证参考（ERFA）

| 用途 | ERFA 例程 |
|---|---|
| 岁差-章动矩阵 | `pnm06a` |
| 章动 Δψ/Δε | `nut06a` |
| 平黄赤交角 | `obl06` |
| 平恒星时 GMST | `gmst06` |
| 真恒星时 GAST | `gst06a` |
| **站心地平坐标（含光行差/视差/折射）** | `atco13` |
| 大地 → 地心 | `gd2gc` |
| 折射模型 | `refco` |
| 时角 ↔ 地平 | `hd2ae` / `ae2hd` |
| ΔT | `dtdb` |
| 月球（**仅可用于查转录；含 −0.7″ 系统偏差，见 §0.2**） | `moon98` |
| 太阳（**仅自洽**，非真值） | `epv00` |

**关键手法**：第 4 步验证归算链时，要把**我们自己算出的月球位置**喂给 `atco13`，
再比它输出的 az/el 与我们的 az/el。这样理论误差两边抵消，隔离出归算误差。

---

## 4. 清单

### 第 1 步 — 拆 workspace（不依赖任何阻塞项）✅ 已完成

- [x] 建 workspace：`Cargo.toml` 加 `[workspace]`，`members = ["crates/*"]`
- [x] `crates/ephemeris`：`jd.rs` → `time.rs`、`astro.rs` 拆为 `angle/sun/moon/frames`、`tables.rs`
- [x] `crates/nongli`：`calendar.rs`、`names.rs`、新建 `terms.rs`、`bin/nongli.rs`
- [x] 时区参数化：`tt_to_beijing_jdn` → `tt_to_jdn(jd_tt, utc_offset_hours)` +
      `jdn_in_offset`；`BEIJING_OFFSET_HOURS` 移到 `nongli`（它是农历的政策）
- [x] `tools/gen_tables.py` 输出改到 `crates/ephemeris/src/tables.rs`
- [x] `tools/*.py`、`probe/*.py` 的 `cargo run` 加 `-p <crate>`；`examples/` 各归其位
- [x] README / `docs/ambiguity.md` / `docs/plan.md` 路径同步

**验收**
- [x] `cargo test --workspace` 全绿：**ephemeris 9 + nongli 16 = 25 个单测 + 1 doctest**，
      重​​构前的 21 个一个不少（逐个核对过），新增 4 个
- [x] `nongli` CLI 对参考集输出**逐字节相同**（sha256 比对，12 个调用场景）
- [x] `ephemeris` 的**代码**中农历概念 0 命中；仅 `lib.rs` 保留两行显式声明边界
- [x] 仍然零依赖（`cargo tree` 只有 `ephemeris` 与 `nongli` 两个节点）
- [x] 依赖方向单向 `nongli → ephemeris`
- [x] `probe/dtreference.py` 与 `crates/ephemeris/src/time.rs` 的 **256 点核对 0 不一致**

**实施中偏离计划的两处（已记录理由）**

1. **没有引入 `Instant`。** 计划第 1 步写了要落 `Instant`，但它此刻没有生产者也没有
   消费者：现有代码里"TT 儒略日 → 北京日"只有一个调用点，参数化足以解决；而
   `from_civil` 要等到第 5 步。此刻建它属于脚手架。**改到第 2 步**，
   与它的生产者（`moon_geocentric`）和消费者一起落。
2. **没有建 `crates/moon` 空壳。** 空 crate 是纯脚手架，等第 5 步有内容再建。

**同时清掉的三处小泄漏**（原计划未列）：`frames.rs` / `time.rs` 里用"对农历无影响"
论证设计取舍的三处注释，改为按精度本身论证。

### 第 2 步 — 月球的 λ、β、Δ（缝就位）

- [ ] `tools/gen_tables.py` 加回月球黄纬表（`probe/moon_test.py` 的 `TB` 60 项）
- [ ] 实现 `moon_geocentric(t) -> Ecliptic`：求和黄经、**黄纬**、距离三列
      （`MOON_LON` 第 6 列 Σr 现在被 `_cr` 丢掉，是现成的）
- [ ] 按 §2 “缝的契约” 五项写接口注释（地心 / 平黄道平黄经 / 已含光行时 / 不含章动 / 度·度·km）
- [ ] 加一个 `#[cfg(test)]` 开关把光行时关掉，以便与 `erfa.moon98` 对拍

**验收**
- [ ] Δ 与 `erfa.moon98` 的 `|r|` 差 < 0.1 km（同一级数，属转录检查）
- [ ] β 与 `erfa.moon98` 的纬度差 < 0.01″
- [ ] **关掉光行时后**，λ 与 `erfa.moon98` 差 < 0.1″；开着时差应在 **−0.70″±0.05″**
      （若不符合，说明 §0.2 的结论在这个实现里不成立）
- [ ] Meeus 例 47.a 的黄经仍逐位一致（回归）

### 第 3 步 — 岁差 / 章动 / 黄赤交角 / 赤道坐标

- [ ] IAU2006 岁差（Capitaine 或 Fukushima-Williams 多项式）
- [ ] 章动：**拟合 `erfa.nut06a`** 得到紧凑级数（合法且可验证，不涉及凭记忆编数据）
- [ ] 平黄赤交角 `obl06` 的对等形式
- [ ] 黄道 → 赤道：`(λ, β) → (α, δ)`

**验收**
- [ ] 章动 Δψ/Δε 对 `nut06a`：**< 0.01″**（现状为 0.9″ 缩写式）
- [ ] 黄赤交角对 `obl06`：< 0.001″
- [ ] 岁差-章动矩阵对 `pnm06a`：1900–2100 内 < 0.01″
- [ ] 赤经赤纬对 `astropy` 的 `GeocentricTrueEcliptic→ICRS` 链 < 0.1″（扣除理论误差后）

### 第 4 步 — 恒星时 / 站心 / 地平坐标 ← **你要的「站心 + 地平坐标」**

- [ ] GMST（`gmst06` 的对等形式）+ GAST（含 Δψ·cos ε）
- [ ] 观测者：大地坐标 (φ, λ, h) → 地心直角坐标
- [ ] 站心视差（月球水平视差 ~57′，不做就错一个月亮直径）
- [ ] 大气折射（`refco` 模型）
- [ ] 时角 → 地平坐标 (A, h)
- [ ] 可选：周日光行差

**验收**
- [ ] GMST / GAST 对 `gmst06` / `gst06a`：< 0.001″
- [ ] **喂同一个位置给 `atco13`**，az/el 差：无折射 < 0.1″，含折射 < 1″
- [ ] 至少覆盖：南半球、赤道、高纬、负海拔各一例；±1000 年边界各一例

### 第 5 步 — 月相（数字版收尾）

- [ ] 相位角、照亮比例、月龄（自上次朔）
- [ ] 视直径、距离
- [ ] 亮边方位角（画月牙朝向，为第 6 步铺路）
- [ ] 出没 / 中天时刻（需迭代；月球移动快，要小心）

**验收**
- [ ] 在朔/上弦/望/下弦四个几何点上：照亮比例分别为 0 / 0.5 / 1 / 0.5（±1e-6）
- [ ] 月龄与 `new_moon_jde` 一致
- [ ] 出没时刻用两套独立方法（高度根 + 迭代）互相验证
- [ ] CLI `moon` 能对任意 (时刻, 地点) 输出全部数字，并对上述各项自检

### 第 6 步 — 天平动 + 月面渲染（"画得对"）

- [ ] 光学天平动（经度方向 + 纬度方向）
- [ ] 物理天平动（另一套级数，Meeus 第 53 章量级）
- [ ] 月面坐标 (l, b) → 月缘方位
- [ ] 渲染：先 Unicode 半块字符（零依赖），再考虑 SVG

**验收**
- [ ] 天平动角度对已知值（如 1987-11-24 等 Meeus 算例）
- [ ] 渲染图在朔/望/上下弦四种情况下肉眼与物理预期一致
- [ ] 若零依赖被破坏，需先记录理由

### 第 7 步 — 更换理论本体（**被第 1 节阻塞**）

- [ ] 拿到数据后：实现 ELP/MPP02，或从 `de440.bsp` 拟合紧凑级数
- [ ] 只替换 `geocentric_ecliptic`，**归算链与 app 一行不动**
- [ ] 用 `de440.bsp` 做真值验证

**验收**
- [ ] 月球地心方向对 DE440：**< 0.01″**（当前 2.9″ RMS）
- [ ] README 更新精度表，并删掉"这是可替换的缝"的免责说明

---

## 5. 明确不做

- **不凭记忆补系数**。没有可验证的数据就不动理论本体。
- **不为了位数去磨归算链**。理论卡在 3″ 时，0.001″ 的归算没有意义（见 0.4）。
- **不提前引入 trait**。等第二个实现出现再抽。
- **不引入运行时依赖**（`chrono`/`ratatui`/`egui`）。时间处理本仓库已有；
  渲染先用零依赖方案。破坏零依赖需单独记录理由。
- **不追 2026±4000 之外的范围**。`docs/ambiguity.md` 已说明 ΔT 才是那里的瓶颈。

---

## 6. 风险

| 风险 | 影响 | 应对 |
|---|---|---|
| 拿不到 ELP/DE440 数据 | 永远停在 2.9″ RMS | 把缝留好；文档写明真实精度，不吹 |
| 第 3 步章动拟合项数过多 | 表体积、编译时间 | 先定目标 0.01″，项数超 300 就重新评估 |
| 第 6 步物理天平动比预期大 | "画得对"拖期 | 第 5 步交付后再评估；可先只做光学天平动 |
| workspace 重构改动面大 | 回归风险 | 第 1 步的验收里"CLI 输出逐字节不变"就是为此 |
