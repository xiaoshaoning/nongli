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
use crate::frames_tables::EPSA;
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
/// 注意"平黄道"是 **VSOP87D/Meeus 的口径**。`erfa.moon98` 的输出是 SOFA 用
/// 它自己那套岁差转入 GCRS 的，再用 IAU2006 转回当日黄道会有 T² 增长的口径差
/// （J2000 处为零，±1000 年 ~0.05″，±4000 年 ~12″）。跨口径比较前请先读
/// `docs/plan.md` 第 2 步末尾的诊断。
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
        self.apparent_with(n).rotate(mean_obliquity(t) + n.deps_deg)
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

/// 平黄赤交角 ε_A，度。即 Fukushima–Williams 的 `EPSA` 角。
pub fn mean_obliquity(t: Instant) -> f64 {
    poly_arcsec(&EPSA, t) / 3600.0
}

/// **真**黄赤交角 ε = ε_A + Δε，度。
pub fn true_obliquity(t: Instant) -> f64 {
    mean_obliquity(t) + nutation(t).deps_deg
}

// 岁差的三个 Fukushima–Williams 角（γ̄、φ̄、ψ̄）**没有暴露**。
//
// `tools/gen_frames.py` 会拟合并验证全部四个角（±200 世纪内对 erfa.pfw06 差 2.3e-8″），
// 但只有 ε_A（EPSA）有人用——`mean_obliquity` 需要它。
//
// 另外三个角是把"当日"化到 J2000/ICRS 时才需要的。目前没有任何消费者：
// 地平坐标、恒星时都在"当日"系里，用不上岁差。等真要做 ICRS 输出
// （例如把月球画到星图上）再让生成器把它们一并输出即可，成本是一行。
//
// 之前这里放过一个 `fukushima_williams()`，靠"它已经验证过了"留着——那是拿
// 沉没成本给死代码找理由，已删。

/// 按 T 的降幂求多项式。系数以**角秒**给出，返回值也是角秒——单位换算由调用方做，
/// 免得像先前那样在两层里各换一次。
fn poly_arcsec(c: &[f64; 6], t: Instant) -> f64 {
    let x = (t.tt_jd() - J2000) / 36525.0;
    c.iter().fold(0.0, |acc, &k| acc * x + k)
}

#[cfg(test)]
mod tests {
    use super::*;

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
