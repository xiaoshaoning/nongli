//! # nongli —— 中国农历编算
//!
//! 按 GB/T 33661-2017《农历的编算和颁行》实现：
//!
//! * 以**北京时间**为标准时间 (4.1)；
//! * 朔日为农历月首日 (4.2)；含冬至之月为十一月 (4.3)；
//! * 十三月之岁置闰，取最先不含中气之月为闰月 (4.4)；
//! * 十一月后第 2 个(不计闰月)月为年首 (4.5)；
//! * 干支纪年、生肖纪年、数序纪月、数序/干支纪日 (第 6 章)。
//!
//! 天文计算全部来自 [`ephemeris`]：本 crate 只负责**编排规则**与命名。
//! 朔、节气、太阳黄经的精度与限制见 `ephemeris` 的文档。
//!
//! ## 例子
//!
//! ```
//! use nongli::{Calendar, DateTime, LunarDate};
//!
//! // 只要日期
//! let d = LunarDate::from_date(2026, 2, 17, Calendar::Gregorian);
//! assert_eq!(d.to_string(), "农历丙午年正月初一");
//! assert_eq!(d.zodiac(), "马");
//!
//! // 日期 + 时刻（北京时间）：农历日以北京时间的 0 时为界
//! let d = LunarDate::from_datetime(DateTime::new(2026, 2, 17, 23, 59, 59.0), Calendar::Gregorian);
//! assert_eq!(d.day, 1);
//! let d = LunarDate::from_datetime(DateTime::new(2026, 2, 18, 0, 0, 0.0), Calendar::Gregorian);
//! assert_eq!(d.day, 2);
//!
//! // UTC 输入
//! let d = LunarDate::from_utc(DateTime::new(2026, 2, 17, 7, 30, 0.0));
//! assert_eq!(d.day, 1); // 北京时间 15:30
//! ```
//!
//! ## 适用范围
//!
//! 2026 ± 4000 年都可给出结果，但可信度随年代递减：太阳/月球位置的级数截断误差
//! 在近现代约 1″ (数秒)，到 ±4000 年约 1′ (数分钟)；ΔT 本身在 ±4000 年有小时量级的
//! 不确定度，会直接平移北京时间的日界。近现代 (±500 年) 与紫金山天文台颁行的
//! 农历一致。
//!
//! 唯一性（"算到多远就没有第二种答案"）的完整实测分析见仓库中的 `docs/ambiguity.md`。

pub mod calendar;
pub mod names;
pub mod terms;

/// 标准 4.1：农历以北京时间为准（东经 120° 标准时）。
///
/// 这是**农历的政策**，因此留在本 crate；`ephemeris` 不预设任何时区。
pub const BEIJING_OFFSET_HOURS: f64 = 8.0;

/// 某个瞬间落在哪一个**农历日**（北京时间的日编号）。
///
/// 农历日以北京时间 0 时为界（标准 3.17），所以朔、节气、冬至的“哪一天”全走这里。
/// 把 +8 收进函数，避免每个调用点都重复一遍这个政策。
#[inline]
pub fn beijing_jdn(t: Instant) -> i64 {
    t.jdn_in_offset(BEIJING_OFFSET_HOURS)
}

pub use calendar::{lunar_from_jdn, new_moon_on_or_before, winter_solstice, winter_solstice_jdn,
                   LunarDate, NewMoon};
pub use names::{day_name, ganzhi, GAN, MONTH_NAMES, TERM_NAMES, ZHI, ZODIAC};
// 这两个出现在本 crate 公开函数的签名里。调用方不该为了能写出参数类型而额外
// 依赖 `ephemeris`；其余星历接口请直接用 `ephemeris`。
pub use ephemeris::{Calendar, DateTime, Instant};
pub use terms::{
    is_mid_term, solar_term, solar_term_index_near, solar_term_on, solar_terms_between, term_index,
    WINTER_SOLSTICE_LONGITUDE,
};
