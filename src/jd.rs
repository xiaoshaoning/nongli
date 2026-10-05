//! 儒略日、公历/儒略历换算、ΔT (TT−UT) 模型。
//!
//! 农历以**北京时间**为准 (GB/T 33661-2017 4.1)，而天文级数以 TT 为时间变量，
//! 因此需要一个 ΔT 模型把 TT 换算到平太阳时。

/// 历法：普通日历（公历）或儒略历。
///
/// 本库默认使用**外推公历**(proleptic Gregorian)——即 1582-10-15 之前也按公历规则
/// 向前外推，这是绝大多数软件的习惯。若需要历史学家使用的儒略历，选 [`Calendar::Julian`]。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Calendar {
    /// 外推公历 (proleptic Gregorian)。
    Gregorian,
    /// 儒略历 (Julian)。
    Julian,
}

/// 北京时间相对 UT 的偏移 (天)。
pub const BEIJING_OFFSET_DAYS: f64 = 8.0 / 24.0;

/// J2000.0 历元 (2000-01-01 12:00 TT) 的儒略日。
pub const J2000: f64 = 2451545.0;

// ---------------------------------------------------------------- 历日 <-> JDN

/// 由公历/儒略历日期取该日 0 时 (UT) 的儒略日数 JDN (整数)。
///
/// 返回的是"日数"，即 JD = JDN − 0.5 对应的那一天的编号。
pub fn jdn_from_ymd(year: i64, month: u32, day: u32, cal: Calendar) -> i64 {
    let a = (14 - month as i64).div_euclid(12);
    let y = year + 4800 - a;
    let m = month as i64 + 12 * a - 3;
    let jdn = day as i64 + (153 * m + 2).div_euclid(5) + 365 * y + y.div_euclid(4)
        - 32083;
    match cal {
        Calendar::Julian => jdn,
        Calendar::Gregorian => jdn - y.div_euclid(100) + y.div_euclid(400) + 38,
    }
}

/// [`jdn_from_ymd`] 的逆运算。
pub fn ymd_from_jdn(jdn: i64, cal: Calendar) -> (i64, u32, u32) {
    let (c, century) = match cal {
        Calendar::Julian => (jdn + 32082, 0),
        Calendar::Gregorian => {
            let a = jdn + 32044;
            let b = (4 * a + 3).div_euclid(146097);
            (a - 146097 * b / 4, 100 * b)
        }
    };
    let d = (4 * c + 3).div_euclid(1461);
    let e = c - 1461 * d / 4;
    let m = (5 * e + 2).div_euclid(153);
    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = century + d - 4800 + m / 10;
    (year, month as u32, day as u32)
}

/// 日编号 → 该日 0 时的儒略日。
///
/// 不含时区信息：时区只是在“儒略日 → 日编号”这一步加上偏移，
/// 见 [`tt_to_beijing_jdn`]。
#[inline]
pub fn jd_from_jdn(jdn: i64) -> f64 {
    jdn as f64 - 0.5
}

/// 儒略日 → 该瞬间所在的日编号。
#[inline]
pub fn jdn_from_jd(jd: f64) -> i64 {
    (jd + 0.5).floor() as i64
}

/// 一个普通日历上的日期时刻。
///
/// 本身不含时区信息——时区由使用它的方法名决定（例如
/// [`LunarDate::from_datetime`](crate::LunarDate::from_datetime) 按北京时间解释，
/// [`LunarDate::from_utc`](crate::LunarDate::from_utc) 按 UTC 解释）。
///
/// 农历日以北京时间的 0 时为界（标准 3.17），所以求农历时只有“落在哪一天”
/// 才要紧，时分秒仅用于跨日判断。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct DateTime {
    pub year: i64,
    pub month: u32,
    /// 日期，1 起算。
    pub day: u32,
    pub hour: u32,
    pub minute: u32,
    /// 秒（允许小数）。
    pub second: f64,
}

impl DateTime {
    pub fn new(year: i64, month: u32, day: u32, hour: u32, minute: u32, second: f64) -> Self {
        Self {
            year,
            month,
            day,
            hour,
            minute,
            second,
        }
    }

    /// 当日 0 时。
    pub fn date(year: i64, month: u32, day: u32) -> Self {
        Self::new(year, month, day, 0, 0, 0.0)
    }

    /// 该时刻当天 0 时之后的分数（天）。
    fn day_fraction(&self) -> f64 {
        (self.hour as f64 * 3600.0 + self.minute as f64 * 60.0 + self.second) / 86400.0
    }

    /// 该时刻在该历法下的日编号。
    pub fn jdn(&self, cal: Calendar) -> i64 {
        jdn_from_ymd(self.year, self.month, self.day, cal)
    }

    /// 该时刻的儒略日（该历法下日编号的 0 时起算）。
    ///
    /// 不含时区：要得到北京时间的日编号，需再加上
    /// [`BEIJING_OFFSET_DAYS`] 后用 [`jdn_from_jd`] 取整。
    pub fn jd(&self, cal: Calendar) -> f64 {
        jd_from_jdn(self.jdn(cal)) + self.day_fraction()
    }
}

// ---------------------------------------------------------------- ΔT

/// 小数年 (由儒略日估算)。
#[inline]
fn decimal_year_from_jd(jd: f64) -> f64 {
    2000.0 + (jd - J2000) / 365.25
}

/// 1973–2025 年**实测** ΔT (TT−UT1)，每年 1 月 1 日取值，单位秒。
///
/// 取自 IERS 的 EOP 序列（通过随包发布的 IERS-B 表）。Espenak–Meeus 多项式
/// 在这段区间最大偏差 5.3 s（rms 1.4 s），而 5 s 足以让一个贴近午夜的朔/节气
/// 落到前一天，故近几十年直接使用实测值。
static DELTA_T_OBSERVED: &[(f64, f64)] = &[
    (1973.0, 43.38), (1974.0, 44.48), (1975.0, 45.48), (1976.0, 46.46), (1977.0, 47.52), (1978.0, 48.53),
    (1979.0, 49.59), (1980.0, 50.54), (1981.0, 51.38), (1982.0, 52.17), (1983.0, 52.96), (1984.0, 53.79),
    (1985.0, 54.34), (1986.0, 54.87), (1987.0, 55.32), (1988.0, 55.82), (1989.0, 56.30), (1990.0, 56.86),
    (1991.0, 57.57), (1992.0, 58.31), (1993.0, 59.12), (1994.0, 59.98), (1995.0, 60.79), (1996.0, 61.63),
    (1997.0, 62.30), (1998.0, 62.97), (1999.0, 63.47), (2000.0, 63.83), (2001.0, 64.09), (2002.0, 64.30),
    (2003.0, 64.47), (2004.0, 64.57), (2005.0, 64.69), (2006.0, 64.85), (2007.0, 65.15), (2008.0, 65.46),
    (2009.0, 65.78), (2010.0, 66.07), (2011.0, 66.32), (2012.0, 66.60), (2013.0, 66.91), (2014.0, 67.28),
    (2015.0, 67.64), (2016.0, 68.10), (2017.0, 68.59), (2018.0, 68.97), (2019.0, 69.22), (2020.0, 69.36),
    (2021.0, 69.36), (2022.0, 69.29), (2023.0, 69.20), (2024.0, 69.18), (2025.0, 69.14),
];

fn observed_delta_t(y: f64) -> f64 {
    let (first, last) = (DELTA_T_OBSERVED[0].0, DELTA_T_OBSERVED[DELTA_T_OBSERVED.len() - 1].0);
    let y = y.clamp(first, last);
    let i = ((y - first) as usize).min(DELTA_T_OBSERVED.len() - 2);
    let (y0, v0) = DELTA_T_OBSERVED[i];
    let (y1, v1) = DELTA_T_OBSERVED[i + 1];
    v0 + (v1 - v0) * (y - y0) / (y1 - y0)
}

/// ΔT = TT − UT，单位秒。
///
/// 1973–2025 年用 IERS 实测值；其余年份用 Espenak & Meeus (NASA, 2006) 的分段
/// 多项式，并在表的两端做**增量衔接**以保证连续（即只取多项式的变化量）。
///
/// **精度说明**：近几十年优于 0.1 s；1900–2005 年优于 1 s；再往外 ΔT 本身
/// 不确定度迅速增大——±4000 年可达小时量级，这是超长期农历日界的主要限制。
pub fn delta_t_seconds(decimal_year: f64) -> f64 {
    let (first, last) = (DELTA_T_OBSERVED[0].0, DELTA_T_OBSERVED[DELTA_T_OBSERVED.len() - 1].0);
    if (first..=last).contains(&decimal_year) {
        return observed_delta_t(decimal_year);
    }
    // 表外：用多项式，但把多项式在表端点的值与实测值对齐（保证连续、不改形状）。
    let anchor = if decimal_year < first { first } else { last };
    observed_delta_t(anchor) + (delta_t_series(decimal_year) - delta_t_series(anchor))
}

/// Espenak & Meeus (2006) 分段多项式，直接返回 ΔT 的绝对值。
fn delta_t_series(decimal_year: f64) -> f64 {
    let y = decimal_year;
    if y < -500.0 {
        let u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0 * u * u;
    }
    if y < 500.0 {
        let u = y / 100.0;
        return 10583.6 - 1014.41 * u + 33.78311 * u.powi(2) - 5.952053 * u.powi(3)
            - 0.1798452 * u.powi(4)
            + 0.022174192 * u.powi(5)
            + 0.0090316521 * u.powi(6);
    }
    if y < 1600.0 {
        let u = (y - 1000.0) / 100.0;
        return 1574.2 - 556.01 * u + 71.23472 * u.powi(2) + 0.319781 * u.powi(3)
            - 0.8503463 * u.powi(4)
            - 0.005050998 * u.powi(5)
            + 0.0083572073 * u.powi(6);
    }
    if y < 1700.0 {
        let t = y - 1600.0;
        return 120.0 - 0.9808 * t - 0.01532 * t * t + t * t * t / 7129.0;
    }
    if y < 1800.0 {
        let t = y - 1700.0;
        return 8.83 + 0.1603 * t - 0.0059285 * t * t + 0.00013336 * t * t * t
            - t.powi(4) / 1174000.0;
    }
    if y < 1860.0 {
        let t = y - 1800.0;
        return 13.72 - 0.332447 * t + 0.0068612 * t * t + 0.0041116 * t.powi(3)
            - 0.00037436 * t.powi(4)
            + 0.0000121272 * t.powi(5)
            - 0.0000001699 * t.powi(6)
            + 0.000000000875 * t.powi(7);
    }
    if y < 1900.0 {
        let t = y - 1860.0;
        return 7.62 + 0.5737 * t - 0.251754 * t * t + 0.01680668 * t.powi(3)
            - 0.0004473624 * t.powi(4)
            + t.powi(5) / 233174.0;
    }
    if y < 1920.0 {
        let t = y - 1900.0;
        return -2.79 + 1.494119 * t - 0.0598939 * t * t + 0.0061966 * t.powi(3)
            - 0.000197 * t.powi(4);
    }
    if y < 1941.0 {
        let t = y - 1920.0;
        return 21.20 + 0.84493 * t - 0.076100 * t * t + 0.0020936 * t.powi(3);
    }
    if y < 1961.0 {
        let t = y - 1950.0;
        return 29.07 + 0.407 * t - t * t / 233.0 + t * t * t / 2547.0;
    }
    if y < 1986.0 {
        let t = y - 1975.0;
        return 45.45 + 1.067 * t - t * t / 260.0 - t * t * t / 718.0;
    }
    if y < 2005.0 {
        let t = y - 2000.0;
        return 63.86 + 0.3345 * t - 0.060374 * t * t + 0.0017275 * t.powi(3)
            + 0.000651814 * t.powi(4)
            + 0.00002373599 * t.powi(5);
    }
    if y < 2050.0 {
        let t = y - 2000.0;
        return 62.92 + 0.32217 * t + 0.005589 * t * t;
    }
    if y < 2150.0 {
        return -20.0 + 32.0 * ((y - 1820.0) / 100.0).powi(2) - 0.5628 * (2150.0 - y);
    }
    let u = (y - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

/// TT 儒略日 → 北京时间日编号 (农历日)。
///
/// 步骤：TT → UT1 (减 ΔT) → 北京时间 (+8 h) → 取整日。
pub fn tt_to_beijing_jdn(jd_tt: f64) -> i64 {
    let dt = delta_t_seconds(decimal_year_from_jd(jd_tt));
    let jd_ut = jd_tt - dt / 86400.0;
    jdn_from_jd(jd_ut + BEIJING_OFFSET_DAYS)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jdn_roundtrip() {
        for (y, m, d) in [
            (2000, 1, 1),
            (1970, 1, 1),
            (1582, 10, 15),
            (1, 1, 1),
            (-1974, 3, 5),
            (6026, 12, 31),
            (2024, 2, 29),
        ] {
            let jdn = jdn_from_ymd(y, m, d, Calendar::Gregorian);
            assert_eq!(ymd_from_jdn(jdn, Calendar::Gregorian), (y, m, d), "y={y}");
            let jj = jdn_from_ymd(y, m, d, Calendar::Julian);
            assert_eq!(ymd_from_jdn(jj, Calendar::Julian), (y, m, d), "julian y={y}");
        }
    }

    #[test]
    fn known_jdn() {
        // 1949-10-01 (GB/T 33661-2017 6.3.2 的干支纪日参考日)
        assert_eq!(jdn_from_ymd(1949, 10, 1, Calendar::Gregorian), 2433191);
        assert_eq!(jdn_from_ymd(2000, 1, 1, Calendar::Gregorian), 2451545);
    }

    #[test]
    fn known_delta_t() {
        // 1973 年之后使用 IERS 实测值
        assert!((delta_t_seconds(2000.0) - 63.83).abs() < 0.1, "{}", delta_t_seconds(2000.0));
        assert!((delta_t_seconds(2020.0) - 69.36).abs() < 0.1, "{}", delta_t_seconds(2020.0));
        assert!((delta_t_seconds(1975.5) - 45.98).abs() < 0.2, "{}", delta_t_seconds(1975.5));
        // 1900 年前后用多项式分支
        assert!((delta_t_seconds(1900.0) + 2.8).abs() < 0.5);
        // 表端须连续
        for y in [1972.9, 1973.0, 1973.1, 2024.9, 2025.0, 2025.1] {
            let a = delta_t_seconds(y - 1e-6);
            let b = delta_t_seconds(y + 1e-6);
            assert!((a - b).abs() < 0.01, "y={y} a={a} b={b}");
        }
    }
}
