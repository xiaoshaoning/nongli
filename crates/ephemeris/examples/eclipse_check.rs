//! 日食/月食的验收数据源，配 `tools/eclipsecheck.py`。
//!
//! 输出两种行（JD 为 TT）：
//!
//! ```text
//! SOLAR <yyyy-mm-dd> <P|C>
//! LUNAR <yyyy-mm-dd> <N|P|T>
//! ```
//!
//! 日期取**世界时（UT1）**的那一天——NASA 目录给的就是它。用 crate 自己的
//! `ymd_from_jdn`，Python 侧不必再抄一遍历法换算。
//!
//! 遍历 1901–2100 的每一个朔与望，给出判断（没有食就不输出）。
//! 与 NASA 的《五千年日月食目录》逐条对日期与类型。

use ephemeris::eclipse::{lunar_eclipse_near, solar_eclipse_near, LunarEclipseKind, SolarEclipseKind};
use ephemeris::{full_moon, new_moon, new_moon_index_near, ymd_from_jdn, Calendar, Instant};

/// 世界时的日期字符串。
fn ut_date(t: Instant) -> String {
    let (y, m, d) = ymd_from_jdn(t.jdn_in_offset(0.0), Calendar::Gregorian);
    format!("{y:04}-{m:02}-{d:02}")
}

fn main() {
    // 1901-01-01 与 2101-01-01（TT 近似即可，判据是比较日期而非时刻）
    const LO: f64 = 2415020.5;
    const HI: f64 = 2488433.5;

    let mut k = new_moon_index_near(Instant::from_tt(LO)) - 2;
    loop {
        let nm = new_moon(k);
        let fm = full_moon(k);
        if nm.tt_jd() > HI {
            break;
        }
        if nm.tt_jd() >= LO {
            if let Some(e) = solar_eclipse_near(nm) {
                let c = match e.kind {
                    SolarEclipseKind::Partial => 'P',
                    SolarEclipseKind::Central => 'C',
                };
                println!("SOLAR {} {c}", ut_date(e.greatest));
            }
            // 望在同一个 k 的下半月（约 k+0.5），用它判月食
            if let Some(e) = lunar_eclipse_near(fm) {
                let c = match e.kind {
                    LunarEclipseKind::Penumbral => 'N',
                    LunarEclipseKind::Partial => 'P',
                    LunarEclipseKind::Total => 'T',
                };
                println!("LUNAR {} {c}", ut_date(e.greatest));
            }
        }
        k += 1;
    }
}
