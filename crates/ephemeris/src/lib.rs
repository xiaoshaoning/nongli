//! # ephemeris —— 通用星历
//!
//! 时间尺度、太阳与月球理论、参考系归算。**这个 crate 不知道任何农历概念**：
//! "冬至是第几个节气"、"正月初一"、"闰月" 都留在 `nongli` 里。
//!
//! 划分的依据是**知识**而非执行顺序：
//!
//! * [`time`] —— 儒略日、历法换算、ΔT。不预设时区。
//! * [`sun`] —— VSOP87D 地球级数、太阳地心视黄经。
//! * [`moon`] —— ELP2000-82B 月球理论。**理论之"缝"**，换更高精度理论只动这里。
//! * [`frames`] —— 参考系归算（岁差、黄赤交角、章动、黄道↔赤道）。
//! * [`observer`] —— 恒星时、站心视差、大气折射、地平坐标。
//! * `angle` —— 角度归化，crate 内部工具。
//!
//! ## 时间尺度约定
//!
//! * **自变量**：天文函数一律要 TT 儒略日。
//! * **“哪一天”的返回**：一律是**民用日编号**（`jdn`），时区由调用方以
//!   `utc_offset_hours` 指定，本 crate 不预设立场。
//! * 两者的换算：[`Instant::jdn_in_offset`]（由 [`Instant`] 的 TT 出发）、
//!   [`jdn_in_offset`]（由 UT1 出发）。
//!
//! 混用 TT 与 UT1 是本领域最常见的错误（月球 0.55″/秒），所以时间量一律走
//! [`Instant`]，它把 TT 与 UT1 绑在一起，"把 TT 当 UT1 传"这类错误因此不存在。
//!
//! ## 精度（必读）
//!
//! 对 **JPL DE421 真值**（1900–2050）实测，中位\|差\| / 最大\|差\|：
//!
//! | | 中位\|差\| | 最大\|差\| |
//! |---|---|---|
//! | 太阳视黄经 | 0.0095″ | 0.025″ |
//! | 月球视黄经 | 0.044″ | 0.11″ |
//! | 月球视黄纬 | 0.011″ | 0.11″ |
//! | 月地距离 | 0.019 km | 0.091 km |
//!
//! **不要拿 `erfa.moon98` 或 astropy 的 builtin 当月球的真值**：它们对 JPL 的误差是
//! 中位 19.8″、最大 43.7″。验证用 `tools/truthcheck.py`（对 JPL 内核）。
//!
//! 超出 1900–2053 是理论外推，未经真值验证；那边的日界由 ΔT 的不确定度主导
//! （±4000 年可达小时量级），见 `docs/plan.md`。

pub mod frames;
pub mod moon;
pub mod observer;
pub mod sun;
pub mod time;

mod angle;

// 级数中的相位常数形如 3.14 / 3.142，并非圆周率，clippy 的近似常数检查在此无意义。
#[allow(clippy::approx_constant)]
pub(crate) mod elp2000_tables;

#[allow(clippy::approx_constant)]
pub(crate) mod frames_tables;

#[allow(clippy::approx_constant)]
pub(crate) mod ecliptic_frame_tables;

#[allow(clippy::approx_constant)]
pub(crate) mod nutation_tables;

#[allow(clippy::approx_constant)]
pub(crate) mod vsop87_tables;

// R2D 有真实消费者（moon crate 用 4 处），故对外导出。D2R 没有——已收成 pub(crate)。
pub use frames::{
    mean_obliquity, nutation, true_obliquity, ApparentEcliptic, Ecliptic, Equatorial, Nutation,
};
pub use angle::R2D;
pub use observer::{gast, gmst, Atmosphere, Horizontal, Observer};
pub use moon::{
    libration, moon_apparent_longitude, moon_geocentric, new_moon, new_moon_index_near, Libration,
};
pub use sun::{
    sun_apparent_longitude, sun_distance_au, sun_geometric_longitude,
    sun_longitude_at,
    MEAN_LONGITUDE_AT_J2000_DEG,
    MEAN_MOTION_DEG_PER_DAY,
};
pub use time::{
    delta_t_seconds, jd_from_jdn, jdn_from_ymd, jdn_in_offset, ymd_from_jdn, Calendar,
    DateTime, Instant, J2000,
};
