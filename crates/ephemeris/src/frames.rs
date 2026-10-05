//! 参考系归算：把理论给出的坐标化到所需的参考系。
//!
//! 这是**一份**知识：从"平黄道平分点"出发，经章动、岁差、黄赤交角，最后到
//! 赤道坐标与站心地平坐标。太阳理论与月球理论都只是往这条链的入口喂数据，
//! 它们不该各自知道这条链。
//!
//! 目前实现了章动与 [`Ecliptic::apparent`]；岁差 IAU2006、黄赤交角与
//! 黄道↔赤道将在后续步骤加入，见 `docs/plan.md` 第 3 步。

use crate::angle::D2R;
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
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Ecliptic {
    /// 黄经，度，[0,360)。
    pub lon_deg: f64,
    /// 黄纬，度。
    pub lat_deg: f64,
    /// 地心距离，km。
    pub distance_km: f64,
}

impl Ecliptic {
    /// 平黄道分点 → **真**分点：加章动黄经 Δψ。
    ///
    /// 两个天体加的是**同一个** Δψ，所以在"朔"（日月黄经相等）的判定中它会自动抵消，
    /// 但求太阳黄经达到某个固定值的时刻（节气）时必须计入。
    pub fn apparent(self, t: Instant) -> Self {
        Self {
            lon_deg: crate::angle::norm360(self.lon_deg + nutation_longitude(t.tt_jd())),
            ..self
        }
    }
}

/// 章动黄经 Δψ (度)。IAU 1980 主项 (Meeus 22.3)。
///
/// 缩写式，与 IAU2006/2000A 的最大偏差 0.9″，对太阳视黄经即约 20 秒的时间误差。
/// 后续会换成从 `erfa.nut06a` 拟合的紧凑级数 (见 docs/plan.md 第 3 步)。
///
/// 接口说明：返回**度**，与坐标类型一致；内部式子是角秒。
pub fn nutation_longitude(jde_tt: f64) -> f64 {
    let t = (jde_tt - J2000) / 36525.0;
    let om = (125.04452 - 1934.136261 * t + 0.0020708 * t * t + t * t * t / 450000.0) * D2R;
    let ls = (280.4665 + 36000.7698 * t) * D2R;
    let lm = (218.3165 + 481267.8813 * t) * D2R;
    (-17.20 * om.sin() - 1.32 * (2.0 * ls).sin() - 0.23 * (2.0 * lm).sin()
        + 0.21 * (2.0 * om).sin())
        / 3600.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 钉住缩写式在几个历元的取值（转录检查）。
    ///
    /// 括号内是同期 `erfa.nut06a`（IAU2006/2000A）的值，用于说明缩写式的
    /// 偏差量级：这几个历元上是 0.1″–0.2″，全时段最大约 0.9″（见 `docs/plan.md`）。
    #[test]
    fn nutation_matches_abridged_formula() {
        for (jde, want) in [
            (2448724.5, 0.004578583), // 1992-04-12.0   erfa +16.6065"
            (2451545.0, -0.003897599), // J2000        erfa −13.9320"
            (2461099.0, 0.001957487), // 2026-02-17   erfa  +7.2272"
        ] {
            let got = nutation_longitude(jde);
            assert!((got - want).abs() < 1e-9, "jde={jde} got={got} want={want}");
        }
    }

    #[test]
    fn apparent_only_shifts_longitude() {
        let e = Ecliptic {
            lon_deg: 10.0,
            lat_deg: -3.0,
            distance_km: 384400.0,
        };
        let a = e.apparent(Instant::from_tt(2451545.0));
        assert_eq!(a.lat_deg, e.lat_deg);
        assert_eq!(a.distance_km, e.distance_km);
        assert!((a.lon_deg - 10.0 - nutation_longitude(2451545.0)).abs() < 1e-12);
    }
}
