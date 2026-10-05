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
//! * [`frames`] —— 参考系归算（黄赤交角、章动、黄道↔赤道）。
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
//! 混用 TT 与 UT1 是本领域最常见的错误（月球 0.55″/秒），所以每个形参名都带
//! 时间尺度后缀。**已知不一致**：`_jde` 与 `_tt_jd` 两种拼法并存，是第 1 步
//! 搬移留下的；第 2 步引入 `Instant` 后会取代这些裸 `f64` 参数，届时统一。
//! 已将就之事记录在 `docs/plan.md` 第 2 步。
//!
//! ## 精度（必读）
//!
//! 月球理论是 Meeus 第 47 章的 60 项截断级数，**绝对精度为 RMS 2.9″、最坏 18.3″**
//! （位置 6.1 / 31.7 km，对 ELP/MPP02，1950–2100）。参考系归算可以做到比这高得多，
//! 但在理论被卡住时提高归算位数不会改善最终结果。详见 `docs/plan.md`。

pub mod frames;
pub mod moon;
pub mod observer;
pub mod sun;
pub mod time;

mod angle;

// 级数中的相位常数形如 3.14 / 3.142，并非圆周率，clippy 的近似常数检查在此无意义。
#[allow(clippy::approx_constant)]
pub(crate) mod tables;

#[allow(clippy::approx_constant)]
pub(crate) mod frames_tables;

#[allow(clippy::approx_constant)]
pub(crate) mod nutation_tables;

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
