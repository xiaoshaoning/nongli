//! 太阳理论：VSOP87D 地球级数与太阳地心视黄经。
//!
//! 输出参考**真**分点（视黄经），可直接用于比较与求根。

use crate::frames::nutation;
use crate::angle::{norm180, norm360, newton, D2R, R2D};
use crate::tables::*;
use crate::time::{Instant, J2000};

/// 太阳平黄经的线性模型：λ ≈ [`MEAN_LONGITUDE_AT_J2000_DEG`]
/// + [`MEAN_MOTION_DEG_PER_DAY`]·(t − J2000) 度。
///
/// 只用于给牛顿迭代提供初值，故忽略偏心率与行星摄动。它是太阳本身的量，
/// 与任何历法无关，所以放在这里；节气那边只是它的使用者之一。
pub const MEAN_LONGITUDE_AT_J2000_DEG: f64 = 280.46646;
/// 见 [`MEAN_LONGITUDE_AT_J2000_DEG`]。
pub const MEAN_MOTION_DEG_PER_DAY: f64 = 0.98564736;

fn vsop_sum(terms: &[(f64, f64, f64)], tau: f64) -> f64 {
    let mut s = 0.0;
    for &(a, b, c) in terms {
        s += a * (b + c * tau).cos();
    }
    s
}

/// 地球日心黄经、黄纬 (弧度) 与向径 (AU)，VSOP87D，参考当日平黄道分点。
///
/// crate 内部用；公开面只保留调用方真正要的量。
///
/// 返回 `(L, B, R)`。L、B 单位弧度，R 单位 AU。
pub(crate) fn earth_heliocentric(jde_tt: f64) -> (f64, f64, f64) {
    let tau = (jde_tt - J2000) / 365250.0;
    let t2 = tau * tau;
    let t3 = t2 * tau;
    let t4 = t3 * tau;
    let t5 = t4 * tau;
    let l = (vsop_sum(EARTH_L0, tau)
        + vsop_sum(EARTH_L1, tau) * tau
        + vsop_sum(EARTH_L2, tau) * t2
        + vsop_sum(EARTH_L3, tau) * t3
        + vsop_sum(EARTH_L4, tau) * t4
        + vsop_sum(EARTH_L5, tau) * t5)
        / 1e8;
    let b = (vsop_sum(EARTH_B0, tau) + vsop_sum(EARTH_B1, tau) * tau) / 1e8;
    let r = (vsop_sum(EARTH_R0, tau)
        + vsop_sum(EARTH_R1, tau) * tau
        + vsop_sum(EARTH_R2, tau) * t2
        + vsop_sum(EARTH_R3, tau) * t3
        + vsop_sum(EARTH_R4, tau) * t4)
        / 1e8;
    (l, b, r)
}

/// 太阳地心**几何**黄经 (度，[0,360)，当日**平**分点起算)。
///
/// 只含 VSOP87 → FK5 改正：**不含**周年光行差、**不含**章动。
///
/// 什么时候要它：算**月相**。相位角是日–月–地的纯几何三角形，而光行差是
/// 观测方向的效应，并不改变被照亮的那部分比例。用视黄经会让相位角差约 20.5″
/// （见 `tools/moonphasecheck.py`）。Meeus 第 48 章用的也正是几何黄经。
pub fn sun_geometric_longitude(t: Instant) -> f64 {
    let jde_tt = t.tt_jd();
    let (l, b, _r) = earth_heliocentric(jde_tt);
    let tc = (jde_tt - J2000) / 36525.0;
    let theta = l + core::f64::consts::PI;
    let lp = theta - (1.397 * tc + 0.00031 * tc * tc) * D2R;
    let dlam = (-0.09033 + 0.03916 * (lp.cos() + lp.sin()) * b.tan()) / 3600.0;
    norm360(theta * R2D + dlam)
}

/// 太阳的周年光行差，度（负值，约 −20.5″/R）。
///
/// 单独拆出来是因为"朔"的判定要用它：章动 Δψ 对日月是同一个量、作差时**精确抵消**，
/// 但光行差只作用于太阳，不抵消。见 [`new_moon`](crate::new_moon)。
pub(crate) fn sun_aberration_deg(t: Instant) -> f64 {
    -20.4898 / earth_heliocentric(t.tt_jd()).2 / 3600.0
}

/// 太阳地心视黄经 (度，[0,360)，真分点起算)。
///
/// 即 [`sun_geometric_longitude`] 再加周年光行差与章动。
pub fn sun_apparent_longitude(t: Instant) -> f64 {
    norm360(sun_geometric_longitude(t) + sun_aberration_deg(t) + nutation(t).dpsi_deg)
}

/// 日地距离，AU。
///
/// 月相角要用到它（相位角与距离有关，不只是黄经差）。
pub fn sun_distance_au(t: Instant) -> f64 {
    earth_heliocentric(t.tt_jd()).2
}

/// 太阳视黄经等于 `target_deg` 的时刻。
///
/// `near` 需落在该解前后约 ±7 天内，否则会收敛到相邻的那一次。
pub fn sun_longitude_at(target_deg: f64, near: Instant) -> Instant {
    let mut t = near.tt_jd();
    // 先用平黄经的变化率把初值拉近，再交给牛顿迭代
    t -= norm180(sun_apparent_longitude(Instant::from_tt(t)) - target_deg) / MEAN_MOTION_DEG_PER_DAY;
    let f = |x: f64| norm180(sun_apparent_longitude(Instant::from_tt(x)) - target_deg);
    Instant::from_tt(newton(t, f, 0.5, 5.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meeus_sun_example() {
        // Meeus 例 25.b: 1992-10-13.0 TD, 视黄经 199°54'21.56"
        let want = 199.0 + 54.0 / 60.0 + 21.56 / 3600.0;
        let got = sun_apparent_longitude(Instant::from_tt(2448908.5));
        assert!((got - want).abs() < 0.002, "got {got} want {want}");
    }

    #[test]
    fn sun_longitude_at_is_exact() {
        for (target, near) in [
            (0.0, 2451545.0),
            (270.0, 2451534.0),
            (285.0, 2451545.0),
            (123.0, 2448724.0),
        ] {
            let t = sun_longitude_at(target, Instant::from_tt(near));
            let d = norm180(sun_apparent_longitude(t) - target);
            assert!(d.abs() < 1e-6, "target={target} d={d}");
        }
    }
}
