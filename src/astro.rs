//! 太阳与月球的地心视黄经 (GB/T 33661-2017 第 5 章的计算模型)。
//!
//! * 太阳：VSOP87D 地球级数 (Meeus 附录 III)，含光行差与岁差改正，得到黄道**真**分点
//!   起算的地心视黄经。
//! * 月球：ELP2000-82B 截断级数 (Meeus 第 47 章)，加 IAU 1980 章动，得到同样的
//!   真分点黄经。
//!
//! 两者参考同一分点，因此"朔"的判定中章动自动抵消。
//!
//! 精度 (与 ERFA/IAU 模型对比，TT 儒略日)：
//! 太阳 ≈ 0.4″ (近现代) ~ 50″ (∓4000 年)；月球 ≈ 0.5″ ~ 10″。章动截断 ≤ 0.9″。

use crate::jd::J2000;
use crate::tables::*;

const D2R: f64 = core::f64::consts::PI / 180.0;
const R2D: f64 = 180.0 / core::f64::consts::PI;

/// 把角度归化到 [0, 360)。
#[inline]
pub fn norm360(x: f64) -> f64 {
    x.rem_euclid(360.0)
}

/// 把角度归化到 (−180, 180]。
#[inline]
pub fn wrap180(x: f64) -> f64 {
    let y = x.rem_euclid(360.0);
    if y > 180.0 {
        y - 360.0
    } else {
        y
    }
}

// ------------------------------------------------------------------ 太阳

fn vsop_sum(terms: &[(f64, f64, f64)], tau: f64) -> f64 {
    let mut s = 0.0;
    for &(a, b, c) in terms {
        s += a * (b + c * tau).cos();
    }
    s
}

/// 地球日心黄经、黄纬 (弧度) 与向径 (AU)，VSOP87D，参考当日平黄道分点。
pub fn earth_heliocentric(jde: f64) -> (f64, f64, f64) {
    let tau = (jde - J2000) / 365250.0;
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

/// 章动黄经 Δψ (度)。IAU 1980 主项 (Meeus 22.3)。
pub fn nutation_longitude(jde: f64) -> f64 {
    let t = (jde - J2000) / 36525.0;
    let om = (125.04452 - 1934.136261 * t + 0.0020708 * t * t + t * t * t / 450000.0) * D2R;
    let ls = (280.4665 + 36000.7698 * t) * D2R;
    let lm = (218.3165 + 481267.8813 * t) * D2R;
    (-17.20 * om.sin() - 1.32 * (2.0 * ls).sin() - 0.23 * (2.0 * lm).sin()
        + 0.21 * (2.0 * om).sin())
        / 3600.0
}

/// 太阳地心视黄经 (度, 真分点起算)。
pub fn sun_apparent_longitude(jde_tt: f64) -> f64 {
    let (l, b, r) = earth_heliocentric(jde_tt);
    let t = (jde_tt - J2000) / 36525.0;
    let theta = l + core::f64::consts::PI;
    // VSOP87 → FK5 改正
    let lp = theta - (1.397 * t + 0.00031 * t * t) * D2R;
    let dlam = (-0.09033 + 0.03916 * (lp.cos() + lp.sin()) * b.tan()) / 3600.0;
    let theta = theta * R2D + dlam;
    // 周年光行差
    let aber = -20.4898 / r / 3600.0;
    norm360(theta + aber + nutation_longitude(jde_tt))
}

// ------------------------------------------------------------------ 月球

/// 月球地心黄经 (度)，参考当日**平**分点 (Meeus 第 47 章)。
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
        let mut w = 1.0;
        match cm {
            1 | -1 => w = e,
            2 | -2 => w = e * e,
            _ => {}
        }
        let arg = (cd as f64 * d + cm as f64 * m + cmp as f64 * mp + cf as f64 * f) * D2R;
        sl += cl as f64 * w * arg.sin();
    }
    sl += 3958.0 * (a1 * D2R).sin()
        + 1962.0 * ((lp - f) * D2R).sin()
        + 318.0 * (a2 * D2R).sin();
    norm360(lp + sl / 1e6)
}

/// 月球地心视黄经 (度, 真分点起算)。
pub fn moon_apparent_longitude(jde_tt: f64) -> f64 {
    norm360(moon_longitude_mean_equinox(jde_tt) + nutation_longitude(jde_tt))
}

// ------------------------------------------------------------------ 求根

/// 求使 `f(t)=0` 的时刻，`f` 返回 (−180,180] 的角度差。
/// 步长限幅以免发散。
fn newton<F: Fn(f64) -> f64>(mut t: f64, f: F, h: f64, max_step: f64) -> f64 {
    for _ in 0..40 {
        let y = f(t);
        if y.abs() < 1e-10 {
            break;
        }
        let dy = (f(t + h) - f(t - h)) / (2.0 * h);
        if dy.abs() < 1e-12 {
            break;
        }
        let mut step = y / dy;
        if step > max_step {
            step = max_step;
        } else if step < -max_step {
            step = -max_step;
        }
        t -= step;
    }
    t
}

// ------------------------------------------------------------------ 朔与节气

/// 太阳平黄经的线性模型：λ ≈ `MEAN_LON_AT_J2000` + `MEAN_MOTION`·(jde − J2000) 度。
/// 只用于给牛顿迭代提供初值，故可忽略偏心率与行星摄动。
const MEAN_LON_AT_J2000: f64 = 280.46646;
const MEAN_MOTION_DEG_PER_DAY: f64 = 0.98564736;
/// 一个节气间隔 15°。
const DEG_PER_TERM: f64 = 15.0;

/// 冬至的太阳地心视黄经，见 GB/T 33661-2017 附录 A 表 A.1 第 1 行。
///
/// [`term_index`] 与 [`crate::calendar::winter_solstice_jde`] 都由它导出，
/// 因此“冬至是第 0 个节气”这个约定只存在一处。
pub const WINTER_SOLSTICE_LONGITUDE: f64 = 270.0;

/// 第 `k` 次朔 (日月黄经相同的时刻) 的 TT 儒略日。
///
/// `k = 0` 对应 2000-01-06 附近；`k` 每 +1 前进一个朔望月。
pub fn new_moon_jde(k: i64) -> f64 {
    let kf = k as f64;
    let t = kf / 1236.85;
    let jde = 2451550.09766 + 29.530588861 * kf + 0.00015437 * t * t
        - 0.000000150 * t * t * t
        + 0.00000000073 * t * t * t * t;
    let f = |x: f64| wrap180(moon_apparent_longitude(x) - sun_apparent_longitude(x));
    newton(jde, f, 0.05, 2.0)
}

/// 离 `jde` 最近的一次朔的序号（允许 ±1 的误差）。
///
/// 与 [`new_moon_jde`] 互逆，用平朔望月估算。
pub fn new_moon_index_near(jde: f64) -> i64 {
    ((jde - 2451550.09766) / 29.530588861).round() as i64
}

/// 离 `jde` 最近的节气的序号 `j`（即太阳黄经 ≈ 15j 度），允许 ±1 的误差。
///
/// 与 [`solar_term_jde`] 互逆。取样、扫描节气序列时用它给定范围。
pub fn solar_term_index_near(jde: f64) -> i64 {
    ((MEAN_MOTION_DEG_PER_DAY * (jde - J2000) + MEAN_LON_AT_J2000) / DEG_PER_TERM).round() as i64
}

/// 太阳地心视黄经等于 `15·j` 度 (mod 360) 的时刻 (TT 儒略日)。
///
/// `j = 18` 附近为 1999 年冬至。
pub fn solar_term_jde(j: i64) -> f64 {
    let target = (DEG_PER_TERM * j as f64).rem_euclid(360.0);
    let jde0 = J2000 + (DEG_PER_TERM * j as f64 - MEAN_LON_AT_J2000) / MEAN_MOTION_DEG_PER_DAY;
    solve_sun_longitude(jde0, target)
}

/// 求解太阳视黄经 = `target` (度)，`jde_approx` 需落在该节气前后约 ±7 天内。
pub fn solve_sun_longitude(jde_approx: f64, target: f64) -> f64 {
    let mut t = jde_approx;
    t -= wrap180(sun_apparent_longitude(t) - target) / MEAN_MOTION_DEG_PER_DAY;
    let f = |x: f64| wrap180(sun_apparent_longitude(x) - target);
    newton(t, f, 0.5, 5.0)
}

/// 节气序号 `j` (太阳黄经 = 15j 度) 在 [`TERM_NAMES`](crate::names::TERM_NAMES)
/// 中的下标：0 = 冬至。
#[inline]
pub fn term_index(j: i64) -> usize {
    (j - (WINTER_SOLSTICE_LONGITUDE / DEG_PER_TERM) as i64).rem_euclid(24) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meeus_sun_example() {
        // Meeus 例 25.b: 1992-10-13.0 TD, 视黄经 199°54'21.56"
        let jde = 2448908.5;
        let want = 199.0 + 54.0 / 60.0 + 21.56 / 3600.0;
        let got = sun_apparent_longitude(jde);
        assert!((got - want).abs() < 0.002, "got {got} want {want}");
    }

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
            let d = wrap180(moon_apparent_longitude(t) - sun_apparent_longitude(t));
            assert!(d.abs() < 1e-6, "k={k} d={d}");
        }
    }

    #[test]
    fn solar_term_is_exact() {
        for j in [18i64, 19, 42, 1000, -500, 48000] {
            let t = solar_term_jde(j);
            let want = (15.0 * j as f64).rem_euclid(360.0);
            let d = wrap180(sun_apparent_longitude(t) - want);
            assert!(d.abs() < 1e-6, "j={j} d={d}");
        }
    }

    #[test]
    fn terms_are_monotone_and_spaced() {
        let mut prev = solar_term_jde(18);
        for j in 19..18 + 24 * 20 {
            let t = solar_term_jde(j);
            let gap = t - prev;
            assert!(gap > 13.5 && gap < 17.0, "j={j} gap={gap}");
            prev = t;
        }
    }
}
