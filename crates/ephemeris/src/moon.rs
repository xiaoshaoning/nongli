//! 月球理论：ELP2000-82B 截断级数 (Meeus 第 47 章)。
//!
//! **这里是整个 crate 的理论之"缝"。** 将来换成 ELP/MPP02 或由 DE440 拟合的
//! 级数时，只需替换 [`moon_geocentric`]，参考系归算与所有调用方不动。
//!
//! 精度：对 ELP/MPP02，1950–2100 为 **RMS 2.9″、最坏 18.3″**
//! （位置 6.1 km / 31.7 km）。见 `docs/plan.md` 第 0 节。

use crate::angle::{norm180, norm360, newton, D2R};
use crate::frames::Ecliptic;
use crate::sun::sun_geometric_longitude;
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

impl Arguments {
    /// 一项的引数 `i·D + j·M + k·M' + l·F`，弧度。
    ///
    /// 经度表与纬度表共用同一个组合方式，所以只写一次。
    fn term_angle(&self, cd: i8, cm: i8, cmp: i8, cf: i8) -> f64 {
        (cd as f64 * self.d
            + cm as f64 * self.m
            + cmp as f64 * self.mprime
            + cf as f64 * self.f)
            * D2R
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
    // 附加项里反复出现这两个，给个短名字免得到处写 a.
    let (mp, f) = (a.mprime, a.f);

    let mut sl = 0.0; // 1e-6 度
    let mut sr = 0.0; // 1e-3 km
    for &(cd, cm, cmp, cf, cl, cr) in MOON_LON {
        let (w, arg) = (eccentricity_factor(cm, a.e), a.term_angle(cd, cm, cmp, cf));
        sl += cl as f64 * w * arg.sin();
        sr += cr as f64 * w * arg.cos();
    }
    // 附加项 (Meeus 47.7)
    sl += 3958.0 * (a.a1 * D2R).sin()
        + 1962.0 * ((a.lprime - f) * D2R).sin()
        + 318.0 * (a.a2 * D2R).sin();

    let mut sb = 0.0; // 1e-6 度
    for &(cd, cm, cmp, cf, cb) in MOON_LAT {
        let (w, arg) = (eccentricity_factor(cm, a.e), a.term_angle(cd, cm, cmp, cf));
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
    // 这里**故意不用**视黄经：朔的条件是日月视黄经相等，而章动 Δψ 对两者是同一个量、
    // 作差时精确抵消。少算两遍 IAU 2000A 章动（约 22 µs），根一模一样。
    // 太阳的周年光行差**不能**省——它只作用于太阳。
    let f = |x: f64| {
        let i = Instant::from_tt(x);
        norm180(
            moon_geocentric(i).lon_deg
                - sun_geometric_longitude(i)
                - crate::sun::sun_aberration_deg(i),
        )
    };
    Instant::from_tt(newton(jde, f, 0.05, 2.0))
}

/// 月球**光学天平动**：从地球看到的月面中心（sub-Earth point）的月面经纬度。
///
/// 月球自转基本均匀，但公转不均匀（轨道偏心），且自转轴与轨道面有夹角，于是我们
/// 能多看一点东、西、南、北——这就是天平动。
///
/// | 分量 | 幅度 | 本实现 |
/// |---|---|---|
/// | 经度方向的光学天平动 | ±7.9° | ✅ |
/// | 纬度方向的光学天平动 | ±6.7° | ✅ |
/// | 物理天平动（月球本体真实摆动） | ±0.04° | ❌ |
/// | 周日天平动（观测者不在圆心） | ±1° | 由 `moon` crate 另行处理 |
///
/// **物理天平动未实现**：幅度 ±0.04°，在月面上约合直径的 0.07%，肉眼看不出；
/// 而 Meeus 第 53 章那套级数需要另一张约 20 项的表。等真要拿月面纹理定位再说。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Libration {
    /// 月面经度 l，度，东正。范围约 ±8°。
    pub lon_deg: f64,
    /// 月面纬度 b，度，北正。范围约 ±7°。
    pub lat_deg: f64,
}

/// 由"观测者到月球"的**黄道方向**算月面中心经纬度。
///
/// 私有：目前只有[`libration`]一个消费者。将来做**周日天平动**（观测者不在地心）
/// 时会需要它——那时再公开，成本是一行。第 3 步就是在这个模式上栽过
/// （`fukushima_williams`，靠"以后要用"留着，已删）。
///
/// `lon_deg` 用当日**平**分点起算的几何黄经——**不要**加章动。实测：Meeus 例 53.a
/// 的 l' 是 −1.206°，加了 Δψ（该历元 16.6″）会得到 −1.201°，不加则 −1.2056°。
/// 原因也说得通：Ω = L' − F 本身就在平分点体系里，两边必须同口径。
fn libration_from_direction(lon_deg: f64, lat_deg: f64, t: Instant) -> Libration {
    use crate::angle::{norm180, D2R, R2D};
    /// 月球赤道对黄道的倾角（Cassini 定律），度。
    const I_DEG: f64 = 1.54242;

    let a = arguments(t.tt_jd());
    // 月球轨道升交点黄经 = L' − F
    let om = (a.lprime - a.f).to_radians();
    let (lam, bet) = (lon_deg * D2R, lat_deg * D2R);
    let inc = I_DEG * D2R;
    let w = lam - om; // 月球黄经与升交点黄经之差

    // Meeus 53.1–53.4
    let big_a =
        (w.sin() * bet.cos() * inc.cos() - bet.sin() * inc.sin()).atan2(w.cos() * bet.cos());
    let l = big_a - a.f.to_radians();
    let b = (-w.sin() * bet.cos() * inc.sin() - bet.sin() * inc.cos())
        .clamp(-1.0, 1.0)
        .asin();

    Libration {
        lon_deg: norm180(l * R2D),
        lat_deg: b * R2D,
    }
}

/// 月球的光学天平动（地心）。
pub fn libration(t: Instant) -> Libration {
    let e = moon_geocentric(t);
    libration_from_direction(e.lon_deg, e.lat_deg, t)
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
        use crate::sun::sun_apparent_longitude;
        for k in [-2000i64, -100, 0, 1, 337, 12000] {
            let t = new_moon(k);
            let d = norm180(moon_apparent_longitude(t) - sun_apparent_longitude(t));
            assert!(d.abs() < 1e-6, "k={k} d={d}");
        }
    }

    /// Meeus 例 53.a：1992-04-12.0 TD。天平动两个分量。
    ///
    /// 注意这里用的是**平**黄经（不加章动）——见 `libration_from_direction` 的说明。
    #[test]
    fn meeus_example_53a() {
        let l = libration(Instant::from_tt(2448724.5));
        assert!((l.lon_deg - -1.206).abs() < 0.002, "l'={}", l.lon_deg);
        assert!((l.lat_deg - 4.194).abs() < 0.002, "b'={}", l.lat_deg);
    }

    /// 天平动的**自洽**检查，不依赖外部参考值：
    ///
    /// * 光经天平动在月球过近地点/远地点附近应接近 0，在两者之间达到极值；
    /// * 光纬天平动在月球过升/降交点附近应接近 0，在两交点之间达到极值；
    /// * 两者幅度都应落在教科书值附近（经度 ±7.9°、纬度 ±6.7°）。
    #[test]
    fn libration_amplitudes_and_phases() {
        // 采样一个交点月与一个近点月
        let mut lo_lon: f64 = 0.0;
        let mut lo_lat: f64 = 0.0;
        let mut n = 0;
        for i in 0..(60 * 24) {
            let t = Instant::from_tt(2451545.0 + i as f64 / 24.0);
            let l = libration(t);
            assert!(l.lon_deg.abs() < 9.0, "l={}", l.lon_deg);
            assert!(l.lat_deg.abs() < 8.0, "b={}", l.lat_deg);
            lo_lon = lo_lon.max(l.lon_deg.abs());
            lo_lat = lo_lat.max(l.lat_deg.abs());
            n += 1;
        }
        assert!(n > 1000);
        // 60 天足以覆盖完整的极值范围
        assert!((7.0..8.5).contains(&lo_lon), "经度天平动幅度 {lo_lon}");
        assert!((5.8..7.2).contains(&lo_lat), "纬度天平动幅度 {lo_lat}");
    }

    /// 天平动的两个分量应当以不同周期变化（一个是近点月、一个是交点月），
    /// 所以 60 天里不应出现"两者同时恒为 0"——那是公式退化的征兆。
    #[test]
    fn libration_is_not_degenerate() {
        let mut both_small = 0;
        for i in 0..(60 * 24) {
            let t = Instant::from_tt(2451545.0 + i as f64 / 24.0);
            let l = libration(t);
            if l.lon_deg.abs() < 0.05 && l.lat_deg.abs() < 0.05 {
                both_small += 1;
            }
        }
        assert!(both_small < 24, "两者同时近零的小时数过多：{both_small}");
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
