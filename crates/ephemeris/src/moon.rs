//! 月球理论：ELP2000-82B 截断级数 (Meeus 第 47 章)。
//!
//! **这里是整个 crate 的理论之"缝"。** 将来换成 ELP/MPP02 或由 DE440 拟合的
//! 级数时，只需要替换本模块，参考系归算与调用方不动。
//!
//! 精度：对 ELP/MPP02，1950–2100 为 **RMS 2.9″、最坏 18.3″**
//! （位置 6.1 km / 31.7 km）。见 `docs/plan.md` 第 0 节。

use crate::frames::nutation_longitude;
use crate::angle::{norm180, norm360, newton, D2R};
use crate::sun::sun_apparent_longitude;
use crate::tables::MOON_LON;
use crate::time::J2000;

/// 平朔望月，日。
const SYNODIC_MONTH_DAYS: f64 = 29.530588861;
/// J2000 附近的一次朔，作为 [`new_moon_jde`] 序号 `k = 0` 的参考。
const EPOCH_NEW_MOON_JDE: f64 = 2451550.09766;

/// 月球地心黄经 (度)，参考当日**平**分点 (Meeus 第 47 章)。
///
/// 这是 Meeus 级数的原样输出：**含**对月球平黄经的光行时改正，**不含**章动。
/// 注意 `erfa.moon98` 与 SOFA 文档明说省略了该光行时项，因此两者之间有一条
/// **−0.70″ 的系统偏差**（见 `docs/plan.md` §0.2）。
pub fn moon_longitude_mean_equinox(jde_tt: f64) -> f64 {
    let t = (jde_tt - J2000) / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let lp = 218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841.0
        - t4 / 65194000.0;
    let d = 297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868.0
        - t4 / 113065000.0;
    let m = 357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000.0;
    let mp = 134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699.0
        - t4 / 14712000.0;
    let f = 93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000.0
        + t4 / 863310000.0;
    let a1 = 119.75 + 131.849 * t;
    let a2 = 53.09 + 479264.290 * t;
    let e = 1.0 - 0.002516 * t - 0.0000074 * t2;

    let mut sl = 0.0;
    for &(cd, cm, cmp, cf, cl, _cr) in MOON_LON {
        let w = match cm {
            1 | -1 => e,
            2 | -2 => e * e,
            _ => 1.0,
        };
        let arg = (cd as f64 * d + cm as f64 * m + cmp as f64 * mp + cf as f64 * f) * D2R;
        sl += cl as f64 * w * arg.sin();
    }
    sl += 3958.0 * (a1 * D2R).sin()
        + 1962.0 * ((lp - f) * D2R).sin()
        + 318.0 * (a2 * D2R).sin();
    norm360(lp + sl / 1e6)
}

/// 月球地心视黄经 (度，[0,360)，真分点起算)。
pub fn moon_apparent_longitude(jde_tt: f64) -> f64 {
    norm360(moon_longitude_mean_equinox(jde_tt) + nutation_longitude(jde_tt))
}

/// 第 `k` 次朔 (日月地心视黄经相等的时刻) 的 TT 儒略日。
///
/// `k = 0` 对应 2000-01-06 附近；`k` 每 +1 前进一个朔望月。
pub fn new_moon_jde(k: i64) -> f64 {
    let kf = k as f64;
    let t = kf / 1236.85;
    let jde = EPOCH_NEW_MOON_JDE
        + SYNODIC_MONTH_DAYS * kf
        + 0.00015437 * t * t
        - 0.000000150 * t * t * t
        + 0.00000000073 * t * t * t * t;
    let f = |x: f64| norm180(moon_apparent_longitude(x) - sun_apparent_longitude(x));
    newton(jde, f, 0.05, 2.0)
}

/// 离 `jde` 最近的一次朔的序号（允许 ±1 的误差）。
///
/// 与 [`new_moon_jde`] 互逆，用平朔望月估算。
pub fn new_moon_index_near(jde: f64) -> i64 {
    ((jde - EPOCH_NEW_MOON_JDE) / SYNODIC_MONTH_DAYS).round() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meeus_moon_example() {
        // Meeus 例 47.a: 1992-04-12.0 TD, 黄经 133.162655°
        let jde = 2448724.5;
        let got = moon_longitude_mean_equinox(jde);
        assert!((got - 133.162655).abs() < 1e-5, "got {got}");
    }

    #[test]
    fn new_moon_is_syzygy() {
        for k in [-2000i64, -100, 0, 1, 337, 12000] {
            let t = new_moon_jde(k);
            let d = norm180(moon_apparent_longitude(t) - sun_apparent_longitude(t));
            assert!(d.abs() < 1e-6, "k={k} d={d}");
        }
    }

    #[test]
    fn new_moon_index_roundtrip() {
        for k in [-2000i64, -1, 0, 1, 337, 12000] {
            assert_eq!(new_moon_index_near(new_moon_jde(k)), k);
        }
    }
}
