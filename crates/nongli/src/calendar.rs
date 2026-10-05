//! 农历编排规则 (GB/T 33661-2017 第 4 章)。
//!
//! * 4.1 以北京时间为标准时间
//! * 4.2 朔日为农历月的第一个农历日
//! * 4.3 包含冬至的农历月为十一月
//! * 4.4 某十一月到下一个十一月(不含)之间有 13 个月时置闰，取最先出现且不含中气的月为闰月
//! * 4.5 十一月之后第 2 个(不计闰月)农历月为年的起始月

use crate::names::{day_name, ganzhi, MONTH_NAMES, ZODIAC};
use crate::terms::{is_mid_term, solar_term_jde, solar_term_on, term_range,
                   WINTER_SOLSTICE_LONGITUDE};
use crate::{tt_to_beijing_jdn, BEIJING_OFFSET_HOURS};
use ephemeris::sun::sun_longitude_at;
use ephemeris::{jd_from_jdn, jdn_from_ymd, jdn_in_offset, new_moon_index_near, new_moon_jde,
                ymd_from_jdn, Calendar, DateTime};

/// 一个农历日期。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct LunarDate {
    /// 农历年 = 该年正月初一所在的公历年份 (6.4 的"年首所在公历年份")。
    pub year: i64,
    /// 农历月，1..=12。
    pub month: u8,
    /// 是否闰月。
    pub leap: bool,
    /// 农历日，1..=30。
    pub day: u8,
    /// 该日对应的北京时间日编号 (儒略日数)。
    pub jdn: i64,
}

/// 干支纪年循环参考：农历甲子年开始于 1984-02-02 (6.1.1)。
const GANZHI_YEAR_EPOCH: i64 = 1984;
/// 干支纪日循环参考：1949-10-01 的农历日为甲子日 (6.3.2)。
const GANZHI_DAY_EPOCH_JDN: i64 = 2433191;

/// 一次朔：序号及它所在的北京时间日。
///
/// 农历月以朔日为首日 (4.2)，故这也是“某月从哪一天开始”的答案。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct NewMoon {
    /// [`new_moon_jde`] 的序号。
    pub index: i64,
    /// 该朔所在的北京时间日编号 (日编号，非儒略日)。
    pub jdn: i64,
}

/// 不晚于 `jdn` 的最后一次朔。
pub fn new_moon_on_or_before(jdn: i64) -> NewMoon {
    // 初值只需落在真值的 ±1 个朔望月内，故直接用日编号当儒略日估算。
    let mut k = new_moon_index_near(jd_from_jdn(jdn));
    for _ in 0..8 {
        let s = tt_to_beijing_jdn(new_moon_jde(k));
        if s > jdn {
            k -= 1;
        } else if tt_to_beijing_jdn(new_moon_jde(k + 1)) <= jdn {
            k += 1;
        } else {
            return NewMoon { index: k, jdn: s };
        }
    }
    NewMoon {
        index: k,
        jdn: tt_to_beijing_jdn(new_moon_jde(k)),
    }
}

/// 冬至的 TT 儒略日，取给定公历年年末的那一次。
pub fn winter_solstice_jde(greg_year: i64) -> f64 {
    let approx = jd_from_jdn(jdn_from_ymd(greg_year, 12, 22, Calendar::Gregorian));
    sun_longitude_at(WINTER_SOLSTICE_LONGITUDE, approx)
}

/// 冬至所在的北京时间日。
pub fn winter_solstice_jdn(greg_year: i64) -> i64 {
    tt_to_beijing_jdn(winter_solstice_jde(greg_year))
}

/// 由"岁"内的月序号求月名序号。
///
/// 序号 0 恒为十一月；跳过闰月时数序递增 1 到 12 后回到 1。闰月沿用前一月序号。
fn month_number(i: usize, leap_idx: Option<usize>) -> (u8, bool) {
    let mut num = 11u8;
    let mut leap = false;
    for idx in 1..=i {
        if Some(idx) == leap_idx {
            leap = true;
        } else {
            num = num % 12 + 1;
            leap = false;
        }
    }
    (num, leap)
}

/// 由北京时间日编号 JDN 求农历日期。
pub fn lunar_from_jdn(jdn: i64) -> LunarDate {
    let (gy, _, _) = ymd_from_jdn(jdn, Calendar::Gregorian);

    // 定位所属的"岁"：从含冬至的十一月到下一个含冬至的十一月(不含)。
    let mut anchor = gy;
    let mut month11 = new_moon_on_or_before(winter_solstice_jdn(anchor));
    if jdn < month11.jdn {
        anchor -= 1;
        month11 = new_moon_on_or_before(winter_solstice_jdn(anchor));
    }
    let k0 = month11.index;
    let m0 = month11.jdn;
    let next_month11 = new_moon_on_or_before(winter_solstice_jdn(anchor + 1));
    let nm = (next_month11.index - k0) as usize;
    let m_end = next_month11.jdn;

    let starts: Vec<i64> = (0..=nm)
        .map(|i| tt_to_beijing_jdn(new_moon_jde(k0 + i as i64)))
        .collect();
    debug_assert_eq!(starts[0], m0);
    debug_assert_eq!(starts[nm], m_end);

    // 4.4 置闰：13 个月的岁中，最先不含中气的月为闰月。
    let leap_idx = if nm == 13 {
        let mut has_zhongqi = vec![false; nm];
        for j in term_range(m0, m_end - 1) {
            if !is_mid_term(j) {
                continue;
            }
            let day = tt_to_beijing_jdn(solar_term_jde(j));
            if day >= m0 && day < m_end {
                let i = starts.partition_point(|&s| s <= day) - 1;
                has_zhongqi[i] = true;
            }
        }
        // 4.4：取最先出现的不含中气之月。13 个月的岁必至少有一个这样的月；
        // 也可能有两个（例：公元 -124 年），此时取最先者。or() 仅为防御性兼容。
        (1..nm).find(|&i| !has_zhongqi[i]).or(Some(nm - 1))
    } else {
        None
    };

    let i = starts.partition_point(|&s| s <= jdn) - 1;
    let (month, leap) = month_number(i, leap_idx);
    let day = (jdn - starts[i] + 1) as u8;

    // 4.5 十一月之后第 2 个(不计闰月)月为年首；农历年以正月初一所在公历年命名。
    let zhengyue = (0..nm)
        .find(|&x| month_number(x, leap_idx).0 == 1)
        .expect("岁内必有正月");
    let zy_start = starts[zhengyue];
    let zy_year = ymd_from_jdn(zy_start, Calendar::Gregorian).0;
    let year = if jdn < zy_start { zy_year - 1 } else { zy_year };

    LunarDate {
        year,
        month,
        leap,
        day,
        jdn,
    }
}

// ------------------------------------------------------------- 公历 -> 农历

impl LunarDate {
    /// 由北京时间日编号求农历日期。
    pub fn from_jdn(jdn: i64) -> Self {
        lunar_from_jdn(jdn)
    }

    /// 由普通日历 (公历) 的北京日期/时刻求农历日期。
    ///
    /// 农历日以北京时间为界 (3.17)，故时刻只影响 [0,24) 内落在哪一天。
    /// 由北京时间的日期时刻求农历日期。
    ///
    /// `dt` 按北京时间解释；农历日以北京时间的 0 时为界 (3.17)，
    /// 因此同一个北京日内的任何时刻结果相同。
    pub fn from_datetime(dt: DateTime, cal: Calendar) -> Self {
        lunar_from_jdn(dt.jdn(cal))
    }

    /// 由普通日历的日期 (0 时) 求农历日期。
    pub fn from_date(year: i64, month: u32, day: u32, cal: Calendar) -> Self {
        Self::from_datetime(DateTime::date(year, month, day), cal)
    }

    /// 由 UTC 日期时刻求农历日期 (换算到北京时间)。
    ///
    /// `dt` 按 UTC 解释，日期部分按公历。北京时间 = UTC + 8 h，
    /// 所以 UTC 时间早于 16:00 时落在同一个北京日，否则落到次日。
    pub fn from_utc(dt: DateTime) -> Self {
        let jd_ut = dt.jd(Calendar::Gregorian);
        lunar_from_jdn(jdn_in_offset(jd_ut, BEIJING_OFFSET_HOURS))
    }

    /// 该日对应的公历日期。
    pub fn gregorian(&self) -> (i64, u32, u32) {
        ymd_from_jdn(self.jdn, Calendar::Gregorian)
    }

    /// 干支纪年 (6.1.1)，如 "丙午"。
    pub fn ganzhi_year(&self) -> String {
        ganzhi(self.year - GANZHI_YEAR_EPOCH)
    }

    /// 生肖纪年 (6.1.2)，如 "马"。
    pub fn zodiac(&self) -> &'static str {
        ZODIAC[(self.year - GANZHI_YEAR_EPOCH).rem_euclid(12) as usize]
    }

    /// 数序纪月 (6.2)，如 "正月" / "闰四月"。
    pub fn month_name(&self) -> String {
        let n = MONTH_NAMES[(self.month - 1) as usize];
        if self.leap {
            format!("闰{n}")
        } else {
            n.to_string()
        }
    }

    /// 数序纪日 (6.3.1)，如 "初一"。
    pub fn day_name(&self) -> String {
        day_name(self.day as u32)
    }

    /// 干支纪日 (6.3.2)，如 "甲子"。
    pub fn ganzhi_day(&self) -> String {
        ganzhi(self.jdn - GANZHI_DAY_EPOCH_JDN)
    }

    /// 这一农历日上若有节气，返回节气名。
    pub fn solar_term(&self) -> Option<&'static str> {
        solar_term_on(self.jdn)
    }
}

impl core::fmt::Display for LunarDate {
    /// 6.4 农历日期的表示方法：农历 + 年名 + 月名 + 日名。
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "农历{}年{}{}",
            self.ganzhi_year(),
            self.month_name(),
            self.day_name()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lunar(y: i64, m: u32, d: u32) -> LunarDate {
        LunarDate::from_date(y, m, d, Calendar::Gregorian)
    }

    /// 春节 (正月初一) 的公历日期。
    #[test]
    fn spring_festivals() {
        let cases = [
            (1949, 1, 29),
            (1976, 1, 31),
            (1984, 2, 2), // 甲子年，标准 6.1.1 的参考
            (1990, 1, 27),
            (1999, 2, 16),
            (2000, 2, 5),
            (2017, 1, 28),
            (2018, 2, 16),
            (2019, 2, 5),
            (2020, 1, 25),
            (2021, 2, 12),
            (2022, 2, 1),
            (2023, 1, 22),
            (2024, 2, 10),
            (2025, 1, 29),
            (2026, 2, 17),
            (2027, 2, 6),
            (2028, 1, 26),
            (2029, 2, 13),
            (2030, 2, 3),
            (2031, 1, 23),
            (2032, 2, 11),
            (2033, 1, 31),
            (2034, 2, 19),
            (2035, 2, 8),
        ];
        for (y, m, d) in cases {
            let l = lunar(y, m, d);
            assert_eq!(
                (l.year, l.month, l.leap, l.day),
                (y, 1, false, 1),
                "{y}-{m}-{d} => {l}"
            );
        }
    }

    #[test]
    fn epoch_references() {
        // 6.1.1 / 6.1.2：1984 农历年为甲子年、鼠年
        let l = lunar(1984, 2, 2);
        assert_eq!(l.ganzhi_year(), "甲子");
        assert_eq!(l.zodiac(), "鼠");
        // 6.3.2：1949-10-01 的农历日为甲子日
        let l = lunar(1949, 10, 1);
        assert_eq!(l.ganzhi_day(), "甲子");
        assert_eq!((l.year, l.month, l.day), (1949, 8, 10));
    }

    #[test]
    fn leap_months() {
        // 已知闰月：起始日 (公历) 与名称
        let cases = [
            (2017, 7, 23, "闰六月"),
            (2020, 5, 23, "闰四月"),
            (2023, 3, 22, "闰二月"),
            (2025, 7, 25, "闰六月"),
            (2033, 12, 22, "闰十一月"), // 著名的 2033 年闰十一月
        ];
        for (y, m, d, want) in cases {
            let l = lunar(y, m, d);
            assert!(l.leap, "{y}-{m}-{d} 应为闰月，得到 {l}");
            assert_eq!(l.month_name(), want, "{y}-{m}-{d}");
            assert_eq!(l.day, 1);
        }
    }

    #[test]
    fn previous_day_is_end_of_previous_month() {
        // 春节前一天 = 上一农历年的十二月(或闰十二月)末日
        for (y, m, d) in [(2026, 2, 17), (2025, 1, 29), (2033, 1, 31)] {
            let eve = lunar(y, m, d - 1);
            assert_eq!(eve.year, y - 1);
            assert!(eve.month == 12 || eve.leap);
        }
    }

    /// 传统节日的农历日期 (附录 B) 对应的公历日期。
    #[test]
    fn festival_dates() {
        let cases = [
            (2023, 6, 22, 5, 5),   // 端午
            (2024, 6, 10, 5, 5),
            (2025, 5, 31, 5, 5),
            (2026, 6, 19, 5, 5),
            (2024, 9, 17, 8, 15),  // 中秋
            (2025, 10, 6, 8, 15),
            (2024, 2, 24, 1, 15),  // 元宵
            (2025, 2, 12, 1, 15),
            (2025, 1, 7, 12, 8),   // 腊八 (甲辰年十二月初八)
        ];
        for (y, m, d, lm, ld) in cases {
            let l = lunar(y, m, d);
            assert_eq!((l.month, l.day), (lm, ld), "{y}-{m}-{d} => {l}");
        }
    }

    /// 清明节 = 包含节气清明在内的农历日。
    #[test]
    fn qingming_is_the_day_containing_the_term() {
        for (y, m, d) in [(2023, 4, 5), (2024, 4, 4), (2025, 4, 4), (2030, 4, 5)] {
            let l = lunar(y, m, d);
            assert_eq!(l.solar_term(), Some("清明"), "{y}-{m}-{d}");
        }
    }

    /// 岁内的月序应连续递增，闰月重复前一月名。
    #[test]
    fn sui_month_numbering_2033() {
        // 2033 年闰十一月：岁内月首与月名
        let expect = [
            ((2033, 11, 22), "十一月"),
            ((2033, 12, 22), "闰十一月"),
            ((2034, 1, 20), "十二月"),
            ((2034, 2, 19), "正月"),
            ((2034, 3, 20), "二月"),
            ((2034, 4, 19), "三月"),
            ((2034, 5, 18), "四月"),
            ((2034, 6, 16), "五月"),
            ((2034, 7, 16), "六月"),
            ((2034, 8, 14), "七月"),
            ((2034, 9, 13), "八月"),
            ((2034, 10, 12), "九月"),
            ((2034, 11, 11), "十月"),
        ];
        for ((y, m, d), name) in expect {
            let l = lunar(y, m, d);
            assert_eq!(l.day, 1, "{y}-{m}-{d}");
            assert_eq!(l.month_name(), name, "{y}-{m}-{d}");
        }
        assert_eq!(lunar(2033, 11, 22).year, 2033);
        assert_eq!(lunar(2034, 2, 19).year, 2034);
    }

    /// 公元 -124 年的岁中有**两个**无中气之月（索引 9 与 12），4.4 规定取最先者，
    /// 故闰月应为索引 9（闰七月），而不是索引 12。
    #[test]
    fn two_zhongqi_free_months_takes_the_first() {
        let ws = winter_solstice_jdn(-124);
        let m11 = new_moon_on_or_before(ws);
        let month_start = |i: i64| tt_to_beijing_jdn(new_moon_jde(m11.index + i));

        let l = lunar_from_jdn(month_start(9));
        assert!(l.leap, "{l}");
        assert_eq!(l.month_name(), "闰七月");

        // 索引 0..8 都不是闰月
        let mut j = m11.jdn;
        while j < month_start(9) {
            assert!(!lunar_from_jdn(j).leap);
            j += 1;
        }
        // 索引 12 的月首（另一个无中气之月）不是闰月
        assert!(!lunar_from_jdn(month_start(12)).leap);
    }

    /// 在整个 2026±4000 年范围内抽样校验编排规则的内部自洽性。
    #[test]
    fn sweep_whole_range_invariants() {
        for year in (-1974..=6026).step_by(37) {
            for &(y, m, d) in &[(year, 1, 1), (year, 6, 15), (year, 12, 31)] {
                let jdn = jdn_from_ymd(y, m, d, Calendar::Gregorian);
                let l = lunar_from_jdn(jdn);
                assert!((1..=12).contains(&l.month), "{y}-{m}-{d} => {l:?}");
                assert!((1..=30).contains(&l.day), "{y}-{m}-{d} => {l:?}");
                assert!(l.year <= y && l.year >= y - 1, "{y}-{m}-{d} => {l:?}");
                // 干支/生肖应在各自周期内
                assert!(matches!(l.ganzhi_year().chars().count(), 2));
                // 月首日应是该月第 1 日；次日必须递增
                let next = lunar_from_jdn(jdn + 1);
                assert!(next.day == 1 || next.day == l.day + 1, "{y}-{m}-{d} => {l:?} {next:?}");
                assert!(next.day != 1 || next.month != l.month || next.leap != l.leap || l.day == 30);
            }
        }
    }

    /// 岁内不变量：12 或 13 个月；13 个月时至少一个月无中气，
    /// 实现取其中最先出现的一个（4.4）。
    #[test]
    fn sweep_sui_invariants() {
        for year in (-1974..=6025).step_by(37) {
            let m11 = new_moon_on_or_before(winter_solstice_jdn(year));
            let next11 = new_moon_on_or_before(winter_solstice_jdn(year + 1));
            let (m0, m1) = (m11.jdn, next11.jdn);
            let k0 = m11.index;
            let nm = (next11.index - k0) as usize;
            assert!(nm == 12 || nm == 13, "year={year} nm={nm}");
            let ws = winter_solstice_jdn(year);
            assert!(m0 <= ws && ws < m1, "year={year}: 冬至不在本岁第一个月内");
            if nm == 13 {
                let starts: Vec<i64> = (0..=nm)
                    .map(|i| tt_to_beijing_jdn(new_moon_jde(k0 + i as i64)))
                    .collect();
                let mut has = vec![false; nm];
                for j in term_range(m0, m1 - 1) {
                    if j.rem_euclid(2) != 0 {
                        continue;
                    }
                    let day = tt_to_beijing_jdn(solar_term_jde(j));
                    if day >= m0 && day < m1 {
                        has[starts.partition_point(|&s| s <= day) - 1] = true;
                    }
                }
                // 十一月含冬至，故索引 0 不为空。
                assert!(has[0], "year={year}: 十一月必须含冬至");
                // 13 个月的岁至少有一个月无中气，否则 4.4 的置闰规则无从适用。
                // 注意可能有**两个**月无中气（例：公元 -124 年），此时取最先者。
                let empty: Vec<usize> = (1..nm).filter(|&i| !has[i]).collect();
                assert!(!empty.is_empty(), "year={year}: 13 个月的岁应有月无中气");
                // 实现取的是最先出现的那个，且只有一个闰月。
                for (idx, &s) in starts[..nm].iter().enumerate() {
                    assert_eq!(
                        lunar_from_jdn(s).leap,
                        idx == empty[0],
                        "year={year} idx={idx}"
                    );
                }
            } else {
                assert_eq!(nm, 12);
            }
        }
    }

    #[test]
    fn consecutive_days_are_consecutive() {
        let start = jdn_from_ymd(2023, 1, 1, Calendar::Gregorian);
        let mut prev = lunar_from_jdn(start);
        for j in start + 1..start + 1500 {
            let cur = lunar_from_jdn(j);
            assert!(
                cur.day == 1 || cur.day == prev.day + 1,
                "j={j} prev={prev:?} cur={cur:?}"
            );
            assert!(
                cur.day == 1 || (cur.month == prev.month && cur.leap == prev.leap),
                "j={j} prev={prev:?} cur={cur:?}"
            );
            assert!(cur.day <= 30, "j={j} cur={cur:?}");
            prev = cur;
        }
    }
}
