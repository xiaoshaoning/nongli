//! 参考系归算：把理论给出的坐标化到所需的参考系。
//!
//! 这是**一份**知识：从"平黄道平分点"出发，经章动、岁差、黄赤交角，最后到
//! 赤道坐标与站心地平坐标。太阳理论与月球理论都只是往这条链的入口喂数据，
//! 它们不该各自知道这条链。
//!
//! 已实现：黄赤交角、章动、[`Ecliptic::apparent`]、[`Ecliptic::equatorial`]。
//! 尚未实现：岁差（把"当日"化到 J2000/ICRS）与地平坐标，见 `docs/plan.md` 第 4 步。
//!
//! 章动用的是 IAU 1980 缩写式；一次"拟合更好级数"的尝试失败，原因记在
//! `nutation` 的文档与 `docs/plan.md` 第 3 步。

use crate::angle::{norm360, D2R, R2D};
use crate::ecliptic_frame_tables::DYNAMICAL_TO_IAU2006_LON_ARCSEC;
use crate::frames_tables::{EPSA, GAMB, PHIB, PSIB};
use crate::nutation_tables::{NUT_LS, NUT_PL};
use crate::time::{Instant, J2000};

/// 地球与月球**地心黄道**坐标。
///
/// 这是两个理论的公共输出类型。约定（**改动会破坏所有下游**）：
///
/// | 项 | 取值 |
/// |---|---|
/// | 参考点 | 地心 |
/// | 参考面/点 | **平**黄道、**平**分点（*不是*真分点）|
/// | 单位 | 经度、纬度：度；距离：km。经度归化到 [0,360) |
///
/// 注意参考面有两个**不同的**"当日平黄道"，别搞混：
///
/// * VSOP87D 与 ELP2000-82B **自己给的**是 1980 年代那套岁差的当日黄道；
/// * 本类型是 **IAU 2006/P03** 的当日黄道，与章动、恒星时、赤道坐标同一口径。
///
/// 两者差 ~0.003″/年（在 J2000 处为零，±1000 年 ~0.3″，±4000 年 ~25″）。
/// 两个理论输出处的换算都由 `frames.rs` 里的 `dynamical_to_iau2006_lon_offset_deg`
/// 完成，所以从本类型出去的东西已经统一在 IAU 2006 里；**新接一个理论时别忘了加**。
///
/// `erfa.moon98` / VSOP87D 的输出是前者，跨口径直接比会看到那个 0.003″/年的
/// 斜率——见 `docs/plan.md` 里“消掉那 0.003″/年的斜率”一节。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Ecliptic {
    /// 黄经，度，[0,360)。
    pub lon_deg: f64,
    /// 黄纬，度。
    pub lat_deg: f64,
    /// 地心距离，km。
    pub distance_km: f64,
}

/// 同一个地心位置，但参考当日**真**分点（已加章动 Δψ）。
///
/// 与 [`Ecliptic`] 的字段一模一样，存在的唯一理由是**消除一个陷阱**：
/// 两者转赤道时一个要加章动、一个不要，光看数值分不出来。第 4 步就因为这个
/// 把平黄道坐标喂进了 `equatorial()`，与 `atco13` 差 12–20″。
/// 用两个类型，编译器替调用方记住这一步。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ApparentEcliptic {
    /// 黄经，度，[0,360)，真分点起算。
    pub lon_deg: f64,
    /// 黄纬，度。
    pub lat_deg: f64,
    /// 地心距离，km。
    pub distance_km: f64,
}

/// **视**赤道坐标（当日真分点、真赤道）。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Equatorial {
    /// 赤经，度，[0,360)。
    pub ra_deg: f64,
    /// 赤纬，度。
    pub dec_deg: f64,
    /// 与 [`Ecliptic`] 同一个距离，原样带过来。
    pub distance_km: f64,
}

impl Ecliptic {
    /// 平黄道分点 → **真**分点：加章动黄经 Δψ。
    ///
    /// 两个天体加的是**同一个** Δψ，所以在"朔"（日月黄经相等）的判定中它会自动抵消，
    /// 但求太阳黄经达到某个固定值的时刻（节气）时必须计入。
    pub fn apparent(self, t: Instant) -> ApparentEcliptic {
        self.apparent_with(nutation(t))
    }

    /// 章动已知时用这个，免得再算一遍（一次 ~11 µs）。
    fn apparent_with(self, n: Nutation) -> ApparentEcliptic {
        ApparentEcliptic {
            lon_deg: norm360(self.lon_deg + n.dpsi_deg),
            lat_deg: self.lat_deg,
            distance_km: self.distance_km,
        }
    }

    /// **平**黄道 → **视**赤道（当日真分点、真赤道）。
    ///
    /// 内部先加章动（[`Ecliptic::apparent`]）再用真黄赤交角 ε = ε_A + Δε，
    /// 所以**传进来的必须是平黄道坐标**，调用方不必先自己调 `apparent()`。
    ///
    /// 这一点曾经是个陷阱：早先 `equatorial()` 直接用真交角，隐含要求调用方先
    /// 加章动，而类型上看不出来——把平黄道坐标喂进去会差 12–20″（用 `atco13`
    /// 对拍时抓到的）。改成"一个函数做完整件事"，预条件就没有了。
    ///
    /// 标准球面三角：
    /// `sin δ = sin β cos ε + cos β sin ε sin λ`，
    /// `α = atan2(sin λ cos ε − tan β sin ε, cos λ)`。
    /// 章动只算一次，同时喂给黄经改正与真黄赤交角（此前两者各算一遍）。
    pub fn equatorial(self, t: Instant) -> Equatorial {
        let n = nutation(t);
        self.apparent_with(n).rotate(true_obliquity_with(t, n))
    }
}

impl ApparentEcliptic {
    /// **视**黄道 → **视**赤道。用真黄赤交角，**不再**加章动。
    ///
    /// 标准球面三角：
    /// `sin δ = sin β cos ε + cos β sin ε sin λ`，
    /// `α = atan2(sin λ cos ε − tan β sin ε, cos λ)`。
    pub fn equatorial(self, t: Instant) -> Equatorial {
        self.rotate(true_obliquity(t))
    }

    /// 球面三角本体。真黄赤交角 ε 由调用方给定，这样在已知章动时不必再算一遍。
    fn rotate(self, eps_deg: f64) -> Equatorial {
        let eps = eps_deg * D2R;
        let (lam, bet) = (self.lon_deg * D2R, self.lat_deg * D2R);
        let (slam, clam) = lam.sin_cos();
        let (sbet, cbet) = bet.sin_cos();
        let (seps, ceps) = eps.sin_cos();
        let dec = (sbet * ceps + cbet * seps * slam).clamp(-1.0, 1.0).asin();
        let y = slam * ceps - (sbet / cbet) * seps;
        Equatorial {
            ra_deg: norm360(y.atan2(clam) * R2D),
            dec_deg: dec * R2D,
            distance_km: self.distance_km,
        }
    }
}

/// 章动黄经 Δψ 与交角章动 Δε，单位**度**。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Nutation {
    pub dpsi_deg: f64,
    pub deps_deg: f64,
}

/// 章动 Δψ / Δε，单位**度**。
///
/// **IAU 2000A 完整模型**（MHB2000 日月章动 678 项 + 行星章动 687 项），
/// 系数与引数取自 ERFA 的 `src/nut00a.c`（见 `third_party/erfa/`），
/// 转录由 `tools/gen_nutation.py` 完成。
///
/// 实测（对 `erfa.nut00a`）：±4000 年 Δψ 最大 **0.0029″**、Δε 最大 **0.0021″**，
/// 见 `tools/framecheck.py`。
///
/// 此前用的是 4 + 4 项缩写式（Δψ 误差 0.34″）。**当初换掉它的理由后来被推翻了**：
/// 本以为那 0.34″ 是太阳视黄经的主要误差源，实测对比 `astropy` 才发现太阳那边是
/// **VSOP87 截断**主导（−1.2″…+1.7″），换章动模型只动 0.1″。换的正当理由另有两条：
/// ① ERA 与恒星时要用到同一套 Δψ，缩写式在 ±40 世纪会差 2327″；
/// ② 专业级定位本来就不该停在一个 4 项缩写式上。
///
/// 代价：每调用一次 284 组三角函数（系数小于 1e-5″ 的项已剪掉，见
/// `tools/gen_nutation.py`），约 11 µs。`new_moon` **没有**因此变慢——那里章动
/// 在朔的判定中作差抵消，压根不算（见 [`crate::moon::new_moon`]）。
pub fn nutation(t: Instant) -> Nutation {
    // 系数单位是 0.1 微角秒：U2R = 角秒→弧度 / 1e7
    const U2R: f64 = core::f64::consts::PI / (180.0 * 3600.0) / 1e7;
    const TURNAS: f64 = 1296000.0;
    const TAU: f64 = core::f64::consts::TAU;

    let tc = (t.tt_jd() - J2000) / 36525.0;
    // 下面这些常数**逐字抄自** ERFA 的 nut00a.c / fa*.c，故意不加千分位下划线——
    // 这样可以直接和源码对diff。（clippy 的 digit grouping 提示在此不适用。）
    //
    // 角秒多项式 → 归化到一圈 → 弧度（对应 C 里的 fmod(...) * DAS2R）
    let red = |as_poly: f64| (as_poly % TURNAS) * (TAU / TURNAS);

    // ---- 日月章动的引数（IERS 2003；其中 l' 与 D 用 MHB2000 自身的形式）----
    let el = red(485868.249036
        + tc * (1717915923.2178
            + tc * (31.8792 + tc * (0.051635 + tc * (-0.00024470)))));
    let elp = red(1287104.79305
        + tc * (129596581.0481
            + tc * (-0.5532 + tc * (0.000136 + tc * (-0.00001149)))));
    let f = red(335779.526232
        + tc * (1739527262.8478
            + tc * (-12.7512 + tc * (-0.001037 + tc * (0.00000417)))));
    let d = red(1072260.70369
        + tc * (1602961601.2090
            + tc * (-6.3706 + tc * (0.006593 + tc * (-0.00003169)))));
    let om = red(450160.398036
        + tc * (-6962890.5431
            + tc * (7.4722 + tc * (0.007702 + tc * (-0.00005939)))));

    let mut dp = 0.0f64;
    let mut de = 0.0f64;
    for &(nl, nlp, nf, nd, nom, sp, spt, cp, ce, cet, se) in NUT_LS {
        let arg = nl as f64 * el
            + nlp as f64 * elp
            + nf as f64 * f
            + nd as f64 * d
            + nom as f64 * om;
        let (s, c) = arg.sin_cos();
        dp += (sp + spt * tc) * s + cp * c;
        de += (ce + cet * tc) * c + se * s;
    }

    // ---- 行星章动：引数用 MHB2000 自己的形式（与上面故意不同，照抄 ERFA）----
    let al = (2.35555598 + 8_328.6914269554 * tc) % TAU;
    let af = (1.627905234 + 8_433.466158131 * tc) % TAU;
    let ad = (5.198466741 + 7_771.3771468121 * tc) % TAU;
    let aom = (2.18243920 - 33.757045 * tc) % TAU;
    let apa = (0.024381750 + 0.00000538691 * tc) * tc;
    let m_lon = |c0: f64, c1: f64| (c0 + c1 * tc) % TAU;
    let alme = m_lon(4.402608842, 2608.7903141574);
    let alve = m_lon(3.176146697, 1021.3285546211);
    let alea = m_lon(1.753470314, 628.3075849991);
    let alma = m_lon(6.203480913, 334.0612426700);
    let alju = m_lon(0.599546497, 52.9690962641);
    let alsa = m_lon(0.874016757, 21.3299104960);
    let alur = m_lon(5.481293872, 7.4781598567);
    let alne = m_lon(5.321159000, 3.8127774000);

    for &(nl, nf, nd, nom, nme, nve, nea, nma, nju, nsa, nur, nne, npa, sp, cp, ce, se) in NUT_PL {
        let arg = nl as f64 * al
            + nf as f64 * af
            + nd as f64 * ad
            + nom as f64 * aom
            + nme as f64 * alme
            + nve as f64 * alve
            + nea as f64 * alea
            + nma as f64 * alma
            + nju as f64 * alju
            + nsa as f64 * alsa
            + nur as f64 * alur
            + nne as f64 * alne
            + npa as f64 * apa;
        let (s, c) = arg.sin_cos();
        dp += sp * s + cp * c;
        de += se * s + ce * c;
    }

    Nutation {
        dpsi_deg: dp * U2R * R2D,
        deps_deg: de * U2R * R2D,
    }
}

/// 从**动力学平黄道分点**到 **IAU 2006 平黄道分点**的黄经改正，度。
///
/// VSOP87D 与 ELP2000-82B 都把黄经给在**它们自带**的"当日动力学平黄道分点"里
/// （1980 年代那套岁差），而本 crate 其余环节用 IAU 2006/P03。两者差一个随时间
/// 平滑增长的微小角量，**在 J2000 处为零、向两侧线性增长**——对 JPL DE421 实测，
/// 太阳视黄经的误差正是它（1900–1925 中位 −0.25″、2000–2025 +0.02″、2025–2050 +0.09″）。
///
/// 系数是拟合出来的，但**不是拿 DE421 拟合**：用的是 VSOP87 的另一个变体
/// （VSOP87B，黄道 J2000）经 IAU 2006 岁差转回当日黄道，与 VSOP87D 作差——
/// 两边同一套理论、同一个地球位置，所以只剩历元约定。交叉验证：一次项
/// −0.30038″/世纪 与 IAU1976−IAU2006 的周年岁差速率差 0.300405″/世纪 差 0.009%。
/// 见 `tools/gen_ecliptic_frame.py`。
///
/// 只改正黄经；两个黄极间的微小倾斜（全程 |Δβ| ≤ 2.1″）未建模，对月球
/// （黄纬可达 ±5°）留 ~0.17″ 残差——比它自己的截断误差小一个量级。
pub(crate) fn dynamical_to_iau2006_lon_offset_deg(t: Instant) -> f64 {
    let tc = (t.tt_jd() - J2000) / 36525.0;
    let (mut sum, mut tp) = (0.0, 1.0);
    for c in DYNAMICAL_TO_IAU2006_LON_ARCSEC {
        sum += c * tp;
        tp *= tc;
    }
    sum / 3600.0
}

/// 平黄赤交角 ε_A，度。即 Fukushima–Williams 的 `EPSA` 角。
pub fn mean_obliquity(t: Instant) -> f64 {
    poly_arcsec(&EPSA, t) / 3600.0
}

/// **真**黄赤交角 ε = ε_A + Δε，度。
pub fn true_obliquity(t: Instant) -> f64 {
    true_obliquity_with(t, nutation(t))
}

/// 章动已知时用这个，免得再算一遍。
///
/// 私有："真交角 = 平交角 + Δε" 这个式子只在这里写一遍，
/// [`Ecliptic::equatorial`] 需要它时也走这里——否则那个式子会有两份，
/// 日后若像 ERFA `ee00a` 那样再加个岁差速率修正就会只改到一份。
fn true_obliquity_with(t: Instant, n: Nutation) -> f64 {
    mean_obliquity(t) + n.deps_deg
}

// 岁差的三个 Fukushima–Williams 角（γ̄、φ̄、ψ̄）**没有暴露**。
//
// `tools/gen_frames.py` 会拟合并验证全部四个角（±200 世纪内对 erfa.pfw06 差 2.3e-8″），
// 但只有 ε_A（EPSA）有人用——`mean_obliquity` 需要它。
//
// 另外三个角是把"当日"化到 J2000/ICRS 时才需要的。目前没有任何消费者：
// 地平坐标、恒星时都在"当日"系里。",
// 第 7 步换 ELP2000-82B 时发现它给的是**黄道 J2000**，于是这四个角终于有了
// 消费者：拼出 fw2m 的旋转矩阵才能把它化到当日。下面那几行就是那时的"一行"。
//
// 之前这里放过一个 `fukushima_williams()`，靠"它已经验证过了"留着——那是拿
// 沉没成本给死代码找理由，已删。

/// 按 T 的降幂求多项式。系数以**角秒**给出，返回值也是角秒——单位换算由调用方做，
/// 免得像先前那样在两层里各换一次。
fn poly_arcsec(c: &[f64; 6], t: Instant) -> f64 {
    let x = (t.tt_jd() - J2000) / 36525.0;
    c.iter().fold(0.0, |acc, &k| acc * x + k)
}

// ===== 把黄道 J2000 化到当日黄道 =====

/// 3×3 矩阵，行主序。只在本模块内部用来拼岁差旋转。
type Mat3 = [[f64; 3]; 3];

/// ERFA 的 `eraRz`：绕 z 轴转 `a`。
///
/// **符号约定与教科书常见的差一个转置**：ERFA 的 `Rx(a)` 是 `[1][2] = +sin`。
/// 这里照抄 ERFA，好与 `erfa.ecm06` 逐个矩阵对拍（踩过一次：自己那个方向的
/// `Rx` 把黄道 J2000 算成了差 30800″ 的东西）。
fn rz(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// ERFA 的 `eraRx`：绕 x 轴转 `a`。符号约定见 [`rz`]。
fn rx(a: f64) -> Mat3 {
    let (s, c) = a.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}

fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut m = [[0.0; 3]; 3];
    for (i, row) in m.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    m
}

fn mat_vec(m: &Mat3, v: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|k| m[i][k] * v[k]).sum())
}

/// 黄道 J2000 → IAU 2006 当日平黄道的旋转矩阵。
///
/// 复合式：`Rx(ε_A) · fw2m(γ̄,φ̄,ψ̄,ε_A) · Rx(−ε₀)`。
/// 中间那个就是 SOFA 的 `pmat06`（= `fw2m(pfw06)`；frame bias 已含在 φ̄ 里，
/// 不再另乘），所以本函数与 `erfa.ecm06` 只差最右边那个 `Rx(−ε₀)`。
fn mat_ecliptic_j2000_to_date(t: Instant) -> Mat3 {
    let as2r = D2R / 3600.0;
    let (gamb, phib, psib, epsa) = (
        poly_arcsec(&GAMB, t) * as2r,
        poly_arcsec(&PHIB, t) * as2r,
        poly_arcsec(&PSIB, t) * as2r,
        poly_arcsec(&EPSA, t) * as2r,
    );
    // fw2m：R = Rx(−ε)·Rz(−ψ)·Rx(φ)·Rz(γ)，每一步都**左乘**
    let mut fw = rz(gamb);
    fw = mat_mul(&rx(phib), &fw);
    fw = mat_mul(&rz(-psib), &fw);
    fw = mat_mul(&rx(-epsa), &fw);
    // 输入是黄道 J2000，先用 J2000 的平黄赤交角转到赤道 J2000
    let eps0 = poly_arcsec(&EPSA, Instant::from_tt(J2000)) * as2r;
    mat_mul(&rx(epsa), &mat_mul(&fw, &rx(-eps0)))
}

/// 累积黄经岁差 p(T)，弧度：J2000 平春分点在**当日**平黄道里的黄经。
///
/// 定义上就是"从 J2000 到当日，春分点在黄道上走了多少"。不另立多项式，
/// 直接拿上面那个**已验证等于 `erfa.ecm06` 复合式**的矩阵作用在 J2000 春分点
/// 方向 (1,0,0) 上——这样它只可能和矩阵一起对或一起错。
pub(crate) fn precession_in_longitude(t: Instant) -> f64 {
    let v = mat_vec(&mat_ecliptic_j2000_to_date(t), [1.0, 0.0, 0.0]);
    v[1].atan2(v[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 黄道 J2000 → IAU 2006 当日平黄道的归算矩阵，对 `erfa.ecm06 ∘ Rx(−ε₀)` 逐元素比。
    ///
    /// 这是整条参考系链的**枢纽**：日月两个理论都靠它落到"当日"。参考值由 erfa 算出
    /// 后写死（测试不依赖 Python/ERFA）。容差 1e-13，即 2e-8″——比拟合 FW 角时
    /// 那个 2.3e-8″ / ±200 世纪的界还紧，因为下面四个历元都在 1900–2100 内。
    ///
    /// 曾经只临时比过一次就把代码删了，而 `docs/compliance.md` 里写着"已对拍"——
    /// 没有可复现的检查就不该那么写。
    #[test]
    fn ecliptic_frame_matrix_matches_erfa() {
        // (JD_TT, erfa.ecm06(jd) @ Rx(−ε₀) 的三行)
        let cases: [(f64, [[f64; 3]; 3]); 4] = [
            (
                2451545.0,
                [
                    [0.999_999_999_999_994_12, -3.289_700_745_181_088e-8, 1.020_704_461_672_871_4e-7],
                    [3.289_700_407_741_964_6e-8, 0.999_999_999_999_998_78, 3.305_943_631_861_770_2e-8],
                    [-1.020_704_472_548_435_5e-7, -3.305_943_294_933_994e-8, 0.999_999_999_999_994_34],
                ],
            ),
            (
                2448724.5,
                [
                    [0.999_998_227_893_309_28, 0.001_882_606_529_697_559_5, 1.701_567_062_546_497_3e-6],
                    [-0.001_882_606_559_283_476_1, 0.999_998_227_740_587_66, 1.755_639_259_978_608_8e-5],
                    [-1.668_512_267_582_219_9e-6, -1.755_956_486_922_894_6e-5, 0.999_999_999_844_438_99],
                ],
            ),
            (
                2415020.5,
                [
                    [0.999_702_954_712_536_01, 0.024_372_150_326_079_516, 2.505_011_197_008_190_1e-5],
                    [-0.024_372_155_384_189_104, 0.999_702_929_204_021_26, 2.266_778_073_128_435e-4],
                    [-1.951_804_471_796_506_1e-5, -2.272_209_989_596_960_3e-4, 0.999_999_973_994_831_43],
                ],
            ),
            (
                2469807.5,
                [
                    [0.999_925_684_309_800_33, -0.012_191_217_242_490_643, -8.928_826_240_743_986_6e-6],
                    [0.012_191_216_151_154_84, 0.999_925_677_921_916, -1.134_951_653_792_149_4e-4],
                    [1.031_180_684_893_211_2e-5, 1.133_778_776_568_902_3e-4, 0.999_999_993_519_561_7],
                ],
            ),
        ];
        for (jd, want) in cases {
            let got = mat_ecliptic_j2000_to_date(Instant::from_tt(jd));
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (got[i][j] - want[i][j]).abs() < 1e-13,
                        "jd={jd} m[{i}][{j}] 得到 {} 参考 {}",
                        got[i][j],
                        want[i][j]
                    );
                }
            }
        }
    }

    /// 章动与黄赤交角对 ERFA 的参考值。
    ///
    /// 容差 **0.001″**：章动现在是 IAU 2000A 完整模型，与 `nut06a` 的差只剩
    /// 它对 2000–2100 的那个小调整因子（~1e-5″）。换回缩写式（0.34″）会立刻挂。
    #[test]
    fn matches_erfa_reference_values() {
        // (JD_TT, nut06a Δψ, nut06a Δε, obl06 ε_A/度)
        let cases = [
            (2451545.0, -13.9320, -5.7694, 23.4392794), // J2000
            (2448724.5, 16.6065, 1.2218, 23.4402841),   // 1992-04-12
            (2461099.0, 7.2272, 8.8653, 23.4358763),    // 2026-02-17
            (2415020.5, 17.4337, -2.2902, 23.4522889),  // 1900-01-00.5
            (2488069.5, 3.2884, 8.5643, 23.4262699),    // ~2100
        ];
        for (jd, dpsi_as, deps_as, eps_deg) in cases {
            let t = Instant::from_tt(jd);
            let n = nutation(t);
            let got_psi = n.dpsi_deg * 3600.0;
            let got_eps = n.deps_deg * 3600.0;
            assert!(
                (got_psi - dpsi_as).abs() < 0.001,
                "jd={jd} Δψ={got_psi} 参考={dpsi_as}"
            );
            assert!(
                (got_eps - deps_as).abs() < 0.001,
                "jd={jd} Δε={got_eps} 参考={deps_as}"
            );
            assert!(
                (mean_obliquity(t) - eps_deg).abs() < 1e-6,
                "jd={jd} ε_A={} 参考={eps_deg}",
                mean_obliquity(t)
            );
        }
    }

    #[test]
    fn apparent_only_shifts_longitude() {
        let e = Ecliptic {
            lon_deg: 10.0,
            lat_deg: -3.0,
            distance_km: 384400.0,
        };
        let t = Instant::from_tt(2451545.0);
        let a = e.apparent(t);
        assert_eq!(a.lat_deg, e.lat_deg);
        assert_eq!(a.distance_km, e.distance_km);
        assert!((a.lon_deg - 10.0 - nutation(t).dpsi_deg).abs() < 1e-12);
    }

    /// 黄道→赤道：`equatorial` 必须**自己**加章动。
    ///
    /// 用不变量来测，而不是把公式再抄一遍：黄道上的点（β=0）转到赤道后，
    /// δ 只由**真**黄经决定——真黄经为 0/180 时 δ=0，为 90/270 时 |δ|=ε。
    /// 若 `equatorial` 忘了加章动，这两条会差约 Δψ·sin ε ≈ 4″。
    #[test]
    fn ecliptic_to_equatorial_applies_nutation() {
        let t = Instant::from_tt(J2000);
        let dpsi = nutation(t).dpsi_deg;
        let eps = true_obliquity(t);

        // 真黄经 = 0 → δ = 0
        let e = Ecliptic { lon_deg: -dpsi, lat_deg: 0.0, distance_km: 1.0 };
        assert!(e.equatorial(t).dec_deg.abs() < 1e-9, "δ={}", e.equatorial(t).dec_deg);

        // 真黄经 = 90 → δ = +ε
        let e = Ecliptic { lon_deg: 90.0 - dpsi, lat_deg: 0.0, distance_km: 1.0 };
        let q = e.equatorial(t);
        assert!((q.dec_deg - eps).abs() < 1e-9, "δ={} ε={eps}", q.dec_deg);

        // 距离原样带过
        let e = Ecliptic { lon_deg: 33.0, lat_deg: -7.0, distance_km: 384400.0 };
        assert_eq!(e.equatorial(t).distance_km, 384400.0);
    }

}
