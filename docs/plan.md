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
| 月球（**仅可用于查转录；含 −0.744″ 系统偏差，见 §0.2**） | `moon98` |
| 月球几何量对拍脚本 | `tools/mooncheck.py` |
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

### 第 2 步 — 月球的 λ、β、Δ（缝就位）✅ 已完成

- [x] `tools/gen_tables.py` 加回月球黄纬表（表 47.B 的 60 项）
- [x] `moon_geocentric(t: Instant) -> Ecliptic`：黄经、**黄纬**、距离三列
      （用上了 `MOON_LON` 第 6 列 Σr，先前被 `_cr` 丢掉）
- [x] 接口注释写明"缝的契约"五项（见下）
- [x] `Instant` 落地；`_jde` 后缀全部消失（`new_moon_jde` → `new_moon` 等）
- [x] 收口公开面：`moon_longitude_mean_equinox` 并入 `moon_geocentric`；
      `earth_heliocentric` 降为 `pub(crate)`

**验收（`tools/mooncheck.py`，4001 点）**

现代化窗口 1500–2200（**决定性的那一档**）：

| 量 | 均值 | rms | 最大 |
|---|---|---|---|
| λ − `erfa.moon98` | −0.7440″ | 0.7440″ | 0.7476″ |
| β − `erfa.moon98` | −0.0000″ | **0.0052″** | 0.0220″ |
| Δ − `erfa.moon98` | +0.0000 km | **0.0004 km** | 0.0024 km |

- [x] Δ < 0.1 km ✓（实际 rms 0.0004 km）
- [x] β < 0.01″ ✓（实际 rms 0.0052″）
- [x] λ 系统差落在 −0.70″±0.05″ ✓（实测 **−0.7438″**，等效 **1.35 s** 光行时）
- [x] Meeus 例 47.a 三个量全对：λ 133.162655°、β −3.229126°、Δ 368409.7 km

λ 那一行值得单独看：**均值、rms、最大值三者几乎相等**（0.7440 / 0.7440 / 0.7476），
说明黄经差**几乎纯粹是那条常数偏移**，围绕它的散布只有约 0.01″。
即：本实现 ≡ `erfa.moon98` + SOFA 省略的那一项光行时改正。§0.2 的结论现在有了数字。

**第 2 步自审发现的一个真 bug（已修，值得记住）**

`angle::newton` 的收敛判据原是 `1e-10` 度，**低于 f64 在该问题上的残差地板**，
于是循环永远达不到判据，每次都空转到 40 次上限。实测：

| | 修前 | 修后 |
|---|---|---|
| `new_moon` | 590 µs/次 | **50.6 µs/次** |
| `nongli` 测试套件 | ~28 s | **3.35 s** |

地板来自表示精度：`jde ≈ 2.45e6` 处相邻可表示数相距 4.7e-10 天，而月球每天走
12.19°，故残差低于约 6e-9° 后，任何步长都无法再改变 `t`。三个修法：
判据放宽到 `1e-8`°、加一条 **`t - step == t` 就退出**（"已到 f64 极限"的唯一可靠
信号）、上限 40 → 12。

**结果的校验和逐位不变**（`4962131845.005`，CLI 12 场景仍逐字节相同）——
说明修掉的只是白算的迭代，没有改任何数值。

教训值得写进注释：**容差不可能小于表示精度**，定得过小等于没有判据。

**第 3 步开工前的诊断：我上一步设的验收靶子是错的**

上一步我写下"把 mooncheck.py 在 2200–6025 的 β 散布从 3.3″ 压到 < 0.1″"，
并推断那是参考系之差、第 3 步能修。**实测证明这个推断错了。**

诊断（`tools/mooncheck.py` + 两个独立实现）：

| 检查 | 结果 |
|---|---|
| Rust β vs Python **独立实现**的 Meeus β | rms **0.000000″** → 转录无误 |
| `erfa.ecm06` vs `astropy` MeanEcliptic | 差 **0.0000″** → 参考没错 |
| 两者方向之间的**角分离** | 现代 ~0.05″，±4000 年 **12″** |

而且角分离 ≈ √(Δλ² + Δβ²)，说明 λ 和 β **都在**变，不是"λ 对、β 不对"。
几何上更关键：**任何刚性转动只要把 β 改了 θ，就必然把 λ 改 θ·tanβ**——
所以"λ 只差 0.01″、β 差 3″"在几何上不可能，这一开始就说明有问题的是口径。

根因在 SOFA 自己的文档里：

> 7) The Meeus algorithm generates position and velocity in **mean ecliptic
>    coordinates of date**, which the present function then **rotates into GCRS**.

即 `erfa.moon98` 的输出是"用 SOFA 挑的那套岁差把 Meeus 黄道坐标转成 GCRS"的结果。
我们再拿 **IAU2006**（`ecm06`/astropy）转回当日黄道，剩下的就是**两套岁差模型之差**：
J2000 处为零、随 T² 增长，与实测的对称增长和量级一致。

**修正**

- **删掉那条靶子**。它要求我们用 Meeus 的黄道去匹配 IAU2006 的黄道，
  是拿参考的实现选择当我们的误差。
- 第 3 步的验收改用**纯参考系量的对拍**（`pfw06` / `obl06` / `nut06a` / `pnm06a`），
  这些是同一口径的直接比较，与月球理论无关。
- 把上述口径差异写进 `moon_geocentric` 的契约注记与 `mooncheck.py` 的输出说明，
  免得下次又有人把它读成误差。

**新发现：远年的退化来自参考系，不来自月球理论**

| 区间 | λ 散布 | β | Δ |
|---|---|---|---|
| 1500–2200 | ~0.01″ | 0.005″ | 0.0004 km |
| −1000..500 | ~0.5″ | 1.87″ | 0.82 km |
| 2200–6025 | ~2.1″ | 3.28″ | 1.91 km |

两侧用的是**同一套** Meeus 级数，所以级数截断误差在作差时抵消；剩下的只能是
**参考系模型之差** —— 本实现沿用 VSOP87D 的"当日平黄道"，astropy 用 IAU2006。
这正是第 3 步要做的，也给了它一个可量化的验收目标。

### 第 3 步 — 岁差 / 章动 / 黄赤交角 / 赤道坐标（部分完成）

- [x] **IAU2006 岁差**（Fukushima–Williams 角）：**对 `erfa.pfw06` 拟合** 5 次多项式，
      不是凭记忆抄。实测 ±200 世纪内与 ERFA 差 **2.3e-8″**（`tools/gen_frames.py`）
- [x] **平黄赤交角 ε_A**：同一组系数的 `EPSA` 角。对 `erfa.obl06` 差 **5e-10 度 ≈ 1.8e-6″**
- [x] **真黄赤交角** ε = ε_A + Δε
- [x] **章动 Δε**：补上了（此前只有 Δψ）
- [x] **黄道 → 赤道**：`Ecliptic::equatorial(t) -> Equatorial`
- [x] 公开面：`nutation` / `mean_obliquity` / `true_obliquity` / `Equatorial` / `Nutation`
- [x] `tools/framecheck.py` + `examples/frame_check.rs`：帧量的逐点对拍

**验收**

| 判据 | 目标 | 实测 |
|---|---|---|
| 章动 Δψ/Δε 对 `nut06a` | < 0.01″ | **未达标**，见下 |
| 黄赤交角对 `obl06` | < 0.001″ | **1.8e-6″** ✅（超 500 倍） |
| 岁差-章动矩阵对 `pnm06a` | < 0.01″ | 未做（`pnm06a` 需要把岁差与章动组合成矩阵，留到第 4 步） |
| 黄道→赤道代数自洽 | — | ✅ 春分/夏至/秋分三点精确 |

#### 章动：目标未达，而且**拟合的方向是错的**

先纠正一处我上一节的误报：我写"缩写式误差 2.1″"，那是**脚本里把年当世纪用了**
（`linspace(-100,100)` 又乘 `36525`，实际是 ±10000 年）。用公平口径重测：

**±1 世纪、0.2 天步长、36.5 万点**（对 `erfa.nut06a`）：

| | Δψ 最大 | Δψ rms | Δε 最大 | Δε rms |
|---|---|---|---|---|
| **缩写式（4+4 项，现在保留的）** | **0.342″** | **0.118″** | **0.087″** | **0.028″** |
| 拟合 200 项 | 0.745″ | 0.219″ | 0.395″ | 0.097″ |

**拟合出来的级数比现有的缩写式差一倍。** 所以没有采用，`frames.rs` 保留缩写式。
这次尝试本身值得记下来，因为它排除了一个看似显然的方向：

- 对 `erfa.nut06a` 做正交匹配追踪，基为五个基本引数 (l, l′, F, D, Ω) 的整数组含，
  跨 ±2000 年用**抖动窗**采样（36000 点、1 天步长）
- 期间踩了四个坑，都记在 `tools/gen_frames.py` 里：
  1. **采样混叠**——基里有 27 天的月球周期，0.2 世纪步长会把它们折到低频，
     拟合出的系数外推时爆到 1e8 角秒；
  2. **基重复**——`c` 与 `−c` 张成同一平面（`sin`/`cos` 的奇偶性），
     两个都放进基让正规方程奇异；
  3. **可分辨性**——短跨度里频率挨得近的两项系数无法分开，会互相抵消，
     一外推就爆（±10 世纪差到 3.6e3″）；
  4. **窗等距混叠**——窗中心等距排列会以 1/间距 的频率再次混叠，必须抖动。
- 修完四个坑后数值稳定了，但精度**在约 0.3″ 处饱和**：项数从 40 加到 320 几乎无改善。
  **瓶颈是引数多项式本身与 IAU2000A 的差异，不是项数**——这一条是收敛不动的。

**结论**：要把章动做到 0.01″，得用 IAU2000A 的**原始表**（约 1300 项）或 IAU1980 的
106 项表，那是**数据**问题，与月球理论本体同一个阻塞项（见第 1 节）。
0.34″ 对现状已经够用：月球理论本体是 2.9″，章动只占 1/8。

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
