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
    pub fn apparent(self, t: Instant) -> Self {
        Self {
            lon_deg: norm360(self.lon_deg + nutation(t).dpsi_deg),
            ..self
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
    pub fn equatorial(self, t: Instant) -> Equatorial {
        let e = self.apparent(t);
        let eps = true_obliquity(t) * D2R;
        let (lam, bet) = (e.lon_deg * D2R, e.lat_deg * D2R);
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

/// 章动 Δψ / Δε。
///
/// 用 **IAU 1980 的缩写式**（4 + 4 项，Meeus 22.3）。实测精度（对 `erfa.nut06a`，
/// ±1 世纪、0.2 天步长、36.5 万点）：
///
/// | 量 | 最大 | rms |
/// |---|---|---|
/// | Δψ | 0.34″ | 0.12″ |
/// | Δε | 0.09″ | 0.03″ |
///
/// **曾经尝试用拟合代替它，失败了。** 见 `tools/gen_frames.py` 与
/// `docs/plan.md` 第 3 步：对 `nut06a` 拟合 200 项（五个基本引数的整数组合作基，
/// 跨 ±2000 年抖动窗采样）得到的级数在 ±1 世纪反而更差（Δψ 0.75″）。
/// 精度在 ~0.3″ 处饱和，加项数无用。所以**保留缩写式**，不引入回归。
pub fn nutation(t: Instant) -> Nutation {
    let jd = t.tt_jd();
    let tc = (jd - J2000) / 36525.0;
    let om = (125.04452 - 1934.136261 * tc + 0.0020708 * tc * tc + tc * tc * tc / 450000.0)
        * D2R;
    let ls = (280.4665 + 36000.7698 * tc) * D2R;
    let lm = (218.3165 + 481267.8813 * tc) * D2R;
    let (s_om, c_om) = om.sin_cos();
    let dpsi = -17.20 * s_om - 1.32 * (2.0 * ls).sin() - 0.23 * (2.0 * lm).sin()
        + 0.21 * (2.0 * om).sin();
    let deps =
        9.20 * c_om + 0.57 * (2.0 * ls).cos() + 0.10 * (2.0 * lm).cos() - 0.09 * (2.0 * om).cos();
    Nutation {
        dpsi_deg: dpsi / 3600.0,
        deps_deg: deps / 3600.0,
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

    /// 与 `erfa.nut06a`、`erfa.obl06` 的参考值比对。
    ///
    /// 参考值由 ERFA 直接给出，单位角秒（Δψ、Δε）与度（ε_A）。
    /// 容忍度取实测上限的约两倍，容不下回归。
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
                (got_psi - dpsi_as).abs() < 1.0,
                "jd={jd} Δψ={got_psi} 参考={dpsi_as}"
            );
            assert!(
                (got_eps - deps_as).abs() < 1.0,
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
