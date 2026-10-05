//! 月球理论：ELP2000-82B 截断级数 (Meeus 第 47 章)。
//!
//! **这里是整个 crate 的理论之"缝"。** 将来换成 ELP/MPP02 或由 DE440 拟合的
//! 级数时，只需替换 [`moon_geocentric`]，参考系归算与所有调用方不动。
//!
//! 精度：对 ELP/MPP02，1950–2100 为 **RMS 2.9″、最坏 18.3″**
//! （位置 6.1 km / 31.7 km）。见 `docs/plan.md` 第 0 节。

use crate::angle::{norm180, norm360, newton, D2R};
use crate::frames::Ecliptic;
use crate::sun::sun_apparent_longitude;
use crate::tables::{MOON_LAT, MOON_LON};
use crate::time::{Instant, J2000};

/// 平朔望月，日。
const SYNODIC_MONTH_DAYS: f64 = 29.530588861;
/// J2000 附近的一次朔，作为 [`new_moon`] 序号 `k = 0` 的参考。
const EPOCH_NEW_MOON_JDE: f64 = 2451550.09766;
/// 地心平均距离，km；Meeus 47 的 Σr 以此为基准。
const MEAN_DISTANCE_KM: f64 = 385000.56;

/// 五个基本引数 (D, M, M', F) 与三个附加引数 A1, A2, A3、偏心率因子 E。
///
/// Meeus 47.1–47.6。单独抽出来是因为经度、纬度、附加项都要用它们，
/// 且都只应在**同一次**求值里算一遍。
struct Arguments {
    lprime: f64,
    d: f64,
    m: f64,
    mprime: f64,
    f: f64,
    a1: f64,
    a2: f64,
    a3: f64,
    e: f64,
}

fn arguments(tt_jd: f64) -> Arguments {
    let t = (tt_jd - J2000) / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    Arguments {
        lprime: 218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841.0
            - t4 / 65194000.0,
        d: 297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868.0
            - t4 / 113065000.0,
        m: 357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000.0,
        mprime: 134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699.0
            - t4 / 14712000.0,
        f: 93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000.0
            + t4 / 863310000.0,
        a1: 119.75 + 131.849 * t,
        a2: 53.09 + 479264.290 * t,
        a3: 313.45 + 481266.484 * t,
        e: 1.0 - 0.002516 * t - 0.0000074 * t2,
    }
}

/// `M` 的幂次对应的偏心率因子 E^|M|（Meeus 47.6）。
fn eccentricity_factor(m: i8, e: f64) -> f64 {
    match m {
        1 | -1 => e,
        2 | -2 => e * e,
        _ => 1.0,
    }
}

/// 月球地心**黄道**位置。
///
/// # 约定（改动会破坏所有下游）
///
/// | 项 | 取值 |
/// |---|---|
/// | 参考点 | **地心** |
/// | 参考面/点 | **平**黄道、**平**分点（*不是*真分点）|
/// | 章动 | **不含**——要真分点用 [`Ecliptic::apparent`] |
/// | 光行时 | **已含**（Meeus 第 47 章所给即如此）|
/// | 单位 | 度、度、km；经度归化到 [0,360) |
///
/// 关于光行时：SOFA 的 `eraMoon98` 文档明说它实现 Meeus 时**省略了对月球平黄经的
/// 光行时改正**，因此本函数的黄经与该例程有约 **−0.70″** 的系统差。这不是误差，
/// 用 `erfa.moon98` 对拍时必须先扣掉（`cargo run -p ephemeris --example moon_check`
/// 会把它量出来）。
pub fn moon_geocentric(t: Instant) -> Ecliptic {
    let a = arguments(t.tt_jd());
    let (d, m, mp, f) = (a.d, a.m, a.mprime, a.f);

    let mut sl = 0.0; // 1e-6 度
    let mut sr = 0.0; // 1e-3 km
    for &(cd, cm, cmp, cf, cl, cr) in MOON_LON {
        let w = eccentricity_factor(cm, a.e);
        let arg = (cd as f64 * d + cm as f64 * m + cmp as f64 * mp + cf as f64 * f) * D2R;
        sl += cl as f64 * w * arg.sin();
        sr += cr as f64 * w * arg.cos();
    }
    // 附加项 (Meeus 47.7)
    sl += 3958.0 * (a.a1 * D2R).sin()
        + 1962.0 * ((a.lprime - f) * D2R).sin()
        + 318.0 * (a.a2 * D2R).sin();

    let mut sb = 0.0; // 1e-6 度
    for &(cd, cm, cmp, cf, cb) in MOON_LAT {
        let w = eccentricity_factor(cm, a.e);
        let arg = (cd as f64 * d + cm as f64 * m + cmp as f64 * mp + cf as f64 * f) * D2R;
        sb += cb as f64 * w * arg.sin();
    }
    // 附加项 (Meeus 47.8)
    sb += -2235.0 * (a.lprime * D2R).sin()
        + 382.0 * (a.a3 * D2R).sin()
        + 175.0 * ((a.a1 - f) * D2R).sin()
        + 175.0 * ((a.a1 + f) * D2R).sin()
        + 127.0 * ((a.lprime - mp) * D2R).sin()
        - 115.0 * ((a.lprime + mp) * D2R).sin();

    Ecliptic {
        lon_deg: norm360(a.lprime + sl / 1e6),
        lat_deg: sb / 1e6,
        distance_km: MEAN_DISTANCE_KM + sr / 1000.0,
    }
}

/// 月球地心视黄经 (度，[0,360)，真分点起算)。
///
/// 等价于 `moon_geocentric(t).apparent(t).lon_deg`，留作便捷入口。
pub fn moon_apparent_longitude(t: Instant) -> f64 {
    moon_geocentric(t).apparent(t).lon_deg
}

/// 第 `k` 次朔 (日月地心视黄经相等的时刻)。
///
/// `k = 0` 对应 2000-01-06 附近；`k` 每 +1 前进一个朔望月。
pub fn new_moon(k: i64) -> Instant {
    let kf = k as f64;
    let t = kf / 1236.85;
    let jde = EPOCH_NEW_MOON_JDE
        + SYNODIC_MONTH_DAYS * kf
        + 0.00015437 * t * t
        - 0.000000150 * t * t * t
        + 0.00000000073 * t * t * t * t;
    let f = |x: f64| {
        let i = Instant::from_tt(x);
        norm180(moon_apparent_longitude(i) - sun_apparent_longitude(i))
    };
    Instant::from_tt(newton(jde, f, 0.05, 2.0))
}

/// 离 `t` 最近的一次朔的序号（允许 ±1 的误差）。
///
/// 与 [`new_moon`] 互逆，用平朔望月估算；`t` 只当作一个粗略的儒略日使用，
/// 因此可以安全地拿民用日编号近似（见 `nongli::new_moon_on_or_before`）。
pub fn new_moon_index_near(t: Instant) -> i64 {
    ((t.tt_jd() - EPOCH_NEW_MOON_JDE) / SYNODIC_MONTH_DAYS).round() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Meeus 例 47.a：1992-04-12.0 TD。三个量都要对上。
    #[test]
    fn meeus_example_47a() {
        let e = moon_geocentric(Instant::from_tt(2448724.5));
        assert!((e.lon_deg - 133.162655).abs() < 1e-5, "λ={}", e.lon_deg);
        assert!((e.lat_deg - -3.229126).abs() < 1e-5, "β={}", e.lat_deg);
        assert!(
            (e.distance_km - 368409.7).abs() < 0.5,
            "Δ={}",
            e.distance_km
        );
    }

    #[test]
    fn new_moon_is_syzygy() {
        for k in [-2000i64, -100, 0, 1, 337, 12000] {
            let t = new_moon(k);
            let d = norm180(moon_apparent_longitude(t) - sun_apparent_longitude(t));
            assert!(d.abs() < 1e-6, "k={k} d={d}");
        }
    }

    #[test]
    fn new_moon_index_roundtrip() {
        for k in [-2000i64, -1, 0, 1, 337, 12000] {
            assert_eq!(new_moon_index_near(new_moon(k)), k);
        }
    }

    /// 距离的物理量级：近地点 ~356 400 km，远地点 ~406 700 km。
    #[test]
    fn distance_stays_in_physical_range() {
        for k in 0..30 {
            let e = moon_geocentric(Instant::from_tt(2451545.0 + k as f64 * 1.0));
            assert!(
                (356_000.0..407_000.0).contains(&e.distance_km),
                "Δ={}",
                e.distance_km
            );
            assert!(e.lat_deg.abs() < 5.5, "β={}", e.lat_deg);
        }
    }
}
