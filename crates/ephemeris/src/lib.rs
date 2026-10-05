//! # ephemeris —— 通用星历
//!
//! 时间尺度、太阳与月球理论、参考系归算。**这个 crate 不知道任何农历概念**：
//! "冬至是第几个节气"、"正月初一"、"闰月" 都留在 `nongli` 里。
//!
//! 划分的依据是**知识**而非执行顺序：
//!
//! * [`time`] —— 儒略日、历法换算、ΔT。不预设时区。
//! * [`sun`] —— VSOP87D 地球级数、太阳地心视黄经。
//! * [`moon`] —— ELP2000-82B 截断级数。**理论之"缝"**，换更高精度理论只动这里。
//! * [`frames`] —— 参考系归算（章动，岁差等后续加入）。
//!
//! ## 时间尺度约定
//!
//! 所有天文函数的自变量一律是 **TT 儒略日**（形参名 `jde_tt` 或 `*_tt_jd`），
//! 所有"哪一天"的返回值一律是**民用日编号**，且时区由调用方以
//! `utc_offset_hours` 指定。两者之间的换算见 [`tt_to_jdn`]。
//!
//! 混用这两种量是本领域最常见的错误，所以参数名一律带 `_tt` 或 `_jdn` 后缀。
//!
//! ## 精度（必读）
//!
//! 月球理论是 Meeus 第 47 章的 60 项截断级数，**绝对精度为 RMS 2.9″、最坏 18.3″**
//! （位置 6.1 / 31.7 km，对 ELP/MPP02，1950–2100）。参考系归算可以做到比这高得多，
//! 但在理论被卡住时提高归算位数不会改善最终结果。详见 `docs/plan.md`。

pub mod frames;
pub mod moon;
pub mod sun;
pub mod time;

pub mod angle;

// 级数中的相位常数形如 3.14 / 3.142，并非圆周率，clippy 的近似常数检查在此无意义。
#[allow(clippy::approx_constant)]
pub(crate) mod tables;

pub use moon::{
    moon_apparent_longitude, moon_longitude_mean_equinox, new_moon_index_near, new_moon_jde,
};
pub use sun::{
    earth_heliocentric, sun_apparent_longitude, sun_longitude_at, MEAN_LONGITUDE_AT_J2000_DEG,
    MEAN_MOTION_DEG_PER_DAY,
};
pub use time::{
    delta_t_seconds, jd_from_jdn, jdn_from_jd, jdn_from_ymd, jdn_in_offset, tt_to_jdn, ymd_from_jdn,
    Calendar, DateTime, J2000,
};
