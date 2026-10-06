//! 日食与月食的**发生与否**与**类型**。
//!
//! 只做**地心几何**的部分：由太阳与月球的地心位置判断这次朔/望是不是食、
//! 以及属于哪一类。**不做**当地见食情况（食分、时刻、食带）——那要指定观测者，
//! 是另一件事（见 `docs/plan.md`）。
//!
//! # 为什么这里可以做得准
//!
//! 判据全是几个角的比较，而太阳与月球的地心位置本 crate 已经验到 0.01″/0.04″
//! （见 `docs/accuracy.md`）。所以这里的分类误差主要来自**几何模型的简化**，
//! 不是来自星历——两者差着两三个数量级。
//!
//! # 简化（都写在各自的地方）
//!
//! * 地球取**球**（R = 6371.0 km），不用椭球。食的判据只关心地球的投影半径，
//!   椭球带来的差 < 0.3%。
//! * **本影不放大**。真实地球大气会把本影放大 ~2%（Danjon），影响的是
//!   "擦边"那几次的偏/全分类，不是食是否发生。本实现用纯几何值，
//!   这与权威目录的差异在 `tools/eclipsecheck.py` 里逐条列出。
//! * 月球与太阳都取圆盘。

use crate::angle::{D2R, R2D};
use crate::moon::moon_geocentric;
use crate::sun::{sun_distance_au, sun_geometric_longitude};
use crate::time::Instant;

/// 地球半径，km（球近似）。
const R_EARTH_KM: f64 = 6371.0;
/// 月球半径，km。
const R_MOON_KM: f64 = 1737.4;
/// 太阳半径，km。
const R_SUN_KM: f64 = 696_000.0;
/// 天文单位，km。
const AU_KM: f64 = 1.495_978_707e8;

/// 月食类型。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LunarEclipseKind {
    /// 半影月食：月球只进地球的半影。
    Penumbral,
    /// 月偏食：月球部分进入本影。
    Partial,
    /// 月全食。
    Total,
}

/// 一次月食。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct LunarEclipse {
    /// 最大食时刻（地心几何意义上两心角距最小的时刻）。
    pub greatest: Instant,
    pub kind: LunarEclipseKind,
    /// 最大食时月心与地影轴的角距，度。
    pub separation_deg: f64,
}

/// 日食类型（**只到中心食/偏食这一层**）。
///
/// 全食/环食/全环食的区分取决于观测者站在哪里（同一场日食在食带不同段落是
/// 不同的），权威目录给的是"最大食点处"的分类。本实现不追那一层，
/// 只判"影轴是否打中地球"——这是**不依赖观测者**的部分。
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SolarEclipseKind {
    /// 偏食：半影扫到地球，但影轴没打中。
    Partial,
    /// 中心食：影轴打中地球（对应目录里的全食/环食/全环食）。
    Central,
}

/// 一次日食。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SolarEclipse {
    /// 最大食时刻（影轴与地球中心最近时刻的近似）。
    pub greatest: Instant,
    pub kind: SolarEclipseKind,
    /// 影轴到地心的最近距离，km。
    pub axis_distance_km: f64,
}

/// 黄道球坐标 (经度°, 纬度°, 距离 km) → 直角坐标 km，当日平黄道分点系。
fn to_xyz(lon_deg: f64, lat_deg: f64, r_km: f64) -> [f64; 3] {
    let (l, b) = (lon_deg * D2R, lat_deg * D2R);
    [r_km * b.cos() * l.cos(), r_km * b.cos() * l.sin(), r_km * b.sin()]
}

/// 太阳的地心直角坐标，km。
///
/// 用**几何**黄经：食的判据是纯几何，加光行差只会平移结果（对分类无影响）。
fn sun_xyz(t: Instant) -> [f64; 3] {
    to_xyz(sun_geometric_longitude(t), 0.0, sun_distance_au(t) * AU_KM)
}

/// 判断 `t` 附近的望是不是月食。`t` 取 [`crate::moon::full_moon`] 给出的时刻。
///
/// 返回 `None` 表示没有食。
pub fn lunar_eclipse_near(t: Instant) -> Option<LunarEclipse> {
    // 两心角距（从地球看，月心到地影轴）在望附近取极小。做几步三分搜索。
    let sep = |x: f64| {
        let i = Instant::from_tt(x);
        let s = sun_xyz(i);
        let d_sun = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        let m = moon_geocentric(i);
        let u = to_xyz(m.lon_deg, m.lat_deg, 1.0);
        // 地影轴方向 = 背离太阳
        let k = -1.0 / d_sun;
        let anti = [s[0] * k, s[1] * k, s[2] * k];
        // 月心方向与影轴方向的夹角
        let dot = (u[0] * anti[0] + u[1] * anti[1] + u[2] * anti[2]).clamp(-1.0, 1.0);
        (dot.acos(), m.distance_km, d_sun)
    };
    let (mut lo, mut hi) = (t.tt_jd() - 0.5, t.tt_jd() + 0.5);
    for _ in 0..60 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if sep(a).0 < sep(b).0 {
            hi = b;
        } else {
            lo = a;
        }
    }
    let (ang, d_moon, d_sun) = sep(0.5 * (lo + hi));
    let greatest = Instant::from_tt(0.5 * (lo + hi));

    // 地影在本影/半影处的线性半径（标准锥几何，纯几何值，未放大）
    let r_umbra = R_EARTH_KM - d_moon * (R_SUN_KM - R_EARTH_KM) / d_sun;
    let r_penumbra = R_EARTH_KM + d_moon * (R_SUN_KM + R_EARTH_KM) / d_sun;
    // 换算成从地心看的角半径
    let (a_u, a_p, a_m) = (
        r_umbra / d_moon,
        r_penumbra / d_moon,
        (R_MOON_KM / d_moon).asin(),
    );

    let kind = if ang <= a_u - a_m {
        Some(LunarEclipseKind::Total)
    } else if ang <= a_u + a_m {
        Some(LunarEclipseKind::Partial)
    } else if ang <= a_p + a_m {
        Some(LunarEclipseKind::Penumbral)
    } else {
        None
    }?;
    Some(LunarEclipse { greatest, kind, separation_deg: ang * R2D })
}

/// 判断 `t` 附近的朔是不是日食。`t` 取 [`crate::moon::new_moon`] 给出的时刻。
pub fn solar_eclipse_near(t: Instant) -> Option<SolarEclipse> {
    // 影轴到地心的距离在朔附近取极小；同时算出该处的半影半径判"扫没扫到"。
    let axis = |x: f64| {
        let i = Instant::from_tt(x);
        let s = sun_xyz(i);
        let m = moon_geocentric(i);
        let mv = to_xyz(m.lon_deg, m.lat_deg, m.distance_km);
        // 影轴方向：从太阳经月球继续向前
        let d = [mv[0] - s[0], mv[1] - s[1], mv[2] - s[2]];
        let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let d = [d[0] / n, d[1] / n, d[2] / n];
        // 轴上离地心最近的点：轴过 mv、方向 d，故参数 t = −(mv·d)
        let t = -(mv[0] * d[0] + mv[1] * d[1] + mv[2] * d[2]);
        let p = [mv[0] + t * d[0], mv[1] + t * d[1], mv[2] + t * d[2]];
        let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        // 该点处的半影半径。**锥是朝远离太阳的方向张开的**，所以用 |t|：
        // 朔时月球在日地之间，t 是负的（最近点在地球这一侧），直接用 t 会把
        // 半影算成收缩的，判出来的偏食少一个量级。
        let r_pen = R_MOON_KM + t.abs() * (R_SUN_KM + R_MOON_KM) / (sun_distance_au(i) * AU_KM);
        (dist, r_pen)
    };
    let (mut lo, mut hi) = (t.tt_jd() - 0.5, t.tt_jd() + 0.5);
    for _ in 0..60 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if axis(a).0 < axis(b).0 {
            hi = b;
        } else {
            lo = a;
        }
    }
    let jd = 0.5 * (lo + hi);
    let (dist, r_pen) = axis(jd);
    if dist > R_EARTH_KM + r_pen {
        return None;
    }
    Some(SolarEclipse {
        greatest: Instant::from_tt(jd),
        kind: if dist < R_EARTH_KM {
            SolarEclipseKind::Central
        } else {
            SolarEclipseKind::Partial
        },
        axis_distance_km: dist,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::moon::{full_moon, new_moon, new_moon_index_near};

    /// 对 **NASA《五千年日月食目录》** 定死的几个历元。
    ///
    /// 真值由 `tools/eclipsecheck.py` 那套解析给出后写死（测试不该依赖把 HTML
    /// 解析一遍）。全量比对见那个工具：1901–2100 的 **452 次日食日期全中**、
    /// 6 次中心/偏判反；**457 次月食配上 451 次**、13 次类型判反。所有差异都是
    /// 同一类——本实现用的是**纯几何**的本影/半影，没有按地球大气放大 ~2%
    /// （Danjon）。那正是"擦边"那几次偏/全/半影分类翻转的原因。
    #[test]
    fn matches_nasa_catalog_samples() {
        // (朔/望附近的 JD, 目录日期, 期望类型)
        let solar = [
            (2451401.5, "1999-08-11", SolarEclipseKind::Central), // 著名的 1999 全食
            (2455565.5, "2011-01-04", SolarEclipseKind::Partial),
        ];
        for (jd, what, want) in solar {
            let k = new_moon_index_near(Instant::from_tt(jd));
            let e = solar_eclipse_near(new_moon(k)).unwrap_or_else(|| panic!("{what} 应为日食"));
            assert_eq!(e.kind, want, "{what}");
        }

        let lunar = [
            (2451918.5, "2001-01-09", LunarEclipseKind::Total),
            (2458680.5, "2019-07-16", LunarEclipseKind::Partial),
            (2458858.5, "2020-01-10", LunarEclipseKind::Penumbral),
            (2459891.5, "2022-11-08", LunarEclipseKind::Total),
        ];
        for (jd, what, want) in lunar {
            let k = new_moon_index_near(Instant::from_tt(jd));
            let e = lunar_eclipse_near(full_moon(k)).unwrap_or_else(|| panic!("{what} 应为月食"));
            assert_eq!(e.kind, want, "{what}");
        }
    }

    /// 非食的朔望不该被误判。朔望每月都有，食一年才几次。
    #[test]
    fn most_syzygies_are_not_eclipses() {
        let (mut s, mut l, mut n) = (0, 0, 0);
        for k in 0..120 {
            n += 1;
            if solar_eclipse_near(new_moon(k)).is_some() {
                s += 1;
            }
            if lunar_eclipse_near(full_moon(k)).is_some() {
                l += 1;
            }
        }
        // 10 年（120 个朔望月）里：日食约 24 次、月食约 24 次，不该是全部
        assert!(s < n / 3 && l < n / 3, "日食 {s} 月食 {l} / {n}——判据太宽了");
        assert!(s > 10 && l > 10, "日食 {s} 月食 {l}——判据太紧了");
    }
}
