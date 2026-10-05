//! 二十四节气 (GB/T 33661-2017 3.13、附录 A)。
//!
//! 节气的定义是通用的天文事实：太阳地心视黄经为 15° 的整数倍。但本模块保留的是
//! **农历特有的约定** —— "冬至是第 0 个节气"、序号与 [`TERM_NAMES`] 下标的对应。
//! 太阳黄经本身由 `ephemeris::sun` 提供。
//!
//! 序号 `j`：太阳视黄经 = 15·`j` 度 (mod 360)。`j` 每 +1 前进一个节气 (~15.2 天)，
//! 故 `j` 是无界的单调整数，**不是** 1..24。

use crate::names::TERM_NAMES;
use crate::beijing_jdn;
use ephemeris::sun::{sun_longitude_at, MEAN_LONGITUDE_AT_J2000_DEG, MEAN_MOTION_DEG_PER_DAY};
use ephemeris::{Instant, J2000};

/// 相邻节气相差 15°。
const DEG_PER_TERM: f64 = 15.0;

/// 冬至的太阳地心视黄经，附录 A 表 A.1 第 1 行。
///
/// [`term_index`] 与 [`crate::calendar::winter_solstice`] 都由它导出，
/// 因此"冬至是第 0 个节气"这个约定只存在一处。
pub const WINTER_SOLSTICE_LONGITUDE: f64 = 270.0;

/// 太阳视黄经等于 `15·j` 度 (mod 360) 的时刻。
///
/// `j = 18` 附近为 1999 年冬至。
pub fn solar_term(j: i64) -> Instant {
    let target = (DEG_PER_TERM * j as f64).rem_euclid(360.0);
    let near = Instant::from_tt(
        J2000 + (DEG_PER_TERM * j as f64 - MEAN_LONGITUDE_AT_J2000_DEG) / MEAN_MOTION_DEG_PER_DAY,
    );
    sun_longitude_at(target, near)
}

/// 离 `t` 最近的节气的序号 `j`（允许 ±1 的误差）。
///
/// 与 [`solar_term`] 互逆，用平黄经估算。取样、扫描节气序列时用它给定范围。
pub fn solar_term_index_near(t: Instant) -> i64 {
    ((MEAN_MOTION_DEG_PER_DAY * (t.tt_jd() - J2000) + MEAN_LONGITUDE_AT_J2000_DEG) / DEG_PER_TERM)
        .round() as i64
}

/// 节气序号 `j` 在 [`TERM_NAMES`] 中的下标：0 = 冬至。
#[inline]
pub fn term_index(j: i64) -> usize {
    (j - (WINTER_SOLSTICE_LONGITUDE / DEG_PER_TERM) as i64).rem_euclid(24) as usize
}

/// 序号 `j` 是否为中气。
///
/// 中气的黄经为 30° 的整数倍，对应偶数 `j` (3.14)。置闰规则 4.4 只关心中气。
#[inline]
pub fn is_mid_term(j: i64) -> bool {
    j.rem_euclid(2) == 0
}

/// 覆盖北京日 `[jdn_lo, jdn_hi]` 的节气序号范围（宽松 ±2）。
pub(crate) fn term_range(jdn_lo: i64, jdn_hi: i64) -> core::ops::RangeInclusive<i64> {
    let idx = |jdn: i64| solar_term_index_near(Instant::from_tt(ephemeris::jd_from_jdn(jdn)));
    (idx(jdn_lo) - 2)..=(idx(jdn_hi) + 2)
}

/// 北京时间日 `jdn` 上若有节气，返回其名称。
pub fn solar_term_on(jdn: i64) -> Option<&'static str> {
    term_range(jdn, jdn)
        .find(|&j| beijing_jdn(solar_term(j)) == jdn)
        .map(|j| TERM_NAMES[term_index(j)])
}

/// `[jdn_lo, jdn_hi]` 区间内的全部节气：`(北京日, 名称)`。
pub fn solar_terms_between(jdn_lo: i64, jdn_hi: i64) -> Vec<(i64, &'static str)> {
    term_range(jdn_lo, jdn_hi)
        .filter_map(|j| {
            let d = beijing_jdn(solar_term(j));
            (d >= jdn_lo && d <= jdn_hi).then(|| (d, TERM_NAMES[term_index(j)]))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ephemeris::sun::sun_apparent_longitude;

    /// 角度差归化到 (−180, 180]。只测试用。
    fn ang_diff(a: f64, b: f64) -> f64 {
        (a - b + 180.0).rem_euclid(360.0) - 180.0
    }

    #[test]
    fn solar_term_is_exact() {
        for j in [18i64, 19, 42, 1000, -500, 48000] {
            let t = solar_term(j);
            let want = (DEG_PER_TERM * j as f64).rem_euclid(360.0);
            let d = ang_diff(sun_apparent_longitude(t), want);
            assert!(d.abs() < 1e-6, "j={j} d={d}");
        }
    }

    #[test]
    fn terms_are_monotone_and_spaced() {
        let mut prev = solar_term(18).tt_jd();
        for j in 19..18 + 24 * 20 {
            let t = solar_term(j).tt_jd();
            let gap = t - prev;
            assert!(gap > 13.5 && gap < 17.0, "j={j} gap={gap}");
            prev = t;
        }
    }

    #[test]
    fn mid_terms_are_even() {
        // 冬至 (270°) 与春分 (0°) 都应是中气
        assert!(is_mid_term(18));
        assert!(is_mid_term(24));
        assert!(!is_mid_term(19));
    }
}
