//! `moon` —— 打印任意时刻、任意地点的月相与月球位置。
//!
//! ```text
//! moon                                  当前位置(默认北京)、当前时刻
//! moon 2026-02-17T20:00:00              指定北京时间
//! moon --lat 51.5 --lon -0.13 --alt 11  指定地点
//! moon --utc 2026-02-17T12:00:00        按 UTC 解释时间
//! moon --days 30                        连续 30 天的月相表
//! ```

use ephemeris::{Atmosphere, Calendar, DateTime, Instant, Observer};
use moon::{observe, rise_set};

/// 默认地点：北京。
const DEFAULT_LAT: f64 = 39.9042;
const DEFAULT_LON: f64 = 116.4074;
const DEFAULT_ALT: f64 = 44.0;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut lat = DEFAULT_LAT;
    let mut lon = DEFAULT_LON;
    let mut alt = DEFAULT_ALT;
    let mut utc = false;
    let mut refraction = true;
    let mut days: Option<i64> = None;
    let mut when: Option<DateTime> = None;

    let mut i = 0;
    while i < args.len() {
        let a = args[i].as_str();
        let num = |i: &mut usize| -> f64 {
            *i += 1;
            args.get(*i)
                .and_then(|s| s.parse().ok())
                .unwrap_or_else(|| {
                    eprintln!("{} 需要一个数值", args[*i - 1]);
                    std::process::exit(2);
                })
        };
        match a {
            "-h" | "--help" => {
                print_help();
                return;
            }
            "--utc" => utc = true,
            "--no-refraction" => refraction = false,
            "--lat" => lat = num(&mut i),
            "--lon" => lon = num(&mut i),
            "--alt" => alt = num(&mut i),
            "--days" => days = Some(num(&mut i) as i64),
            x if x.starts_with("--") => {
                eprintln!("未知选项: {x}");
                std::process::exit(2);
            }
            x => when = Some(parse_datetime(x)),
        }
        i += 1;
    }

    let observer = Observer {
        lat_deg: lat,
        lon_deg: lon,
        height_m: alt,
    };
    let atm = if refraction {
        Some(Atmosphere::default())
    } else {
        None
    };
    let t0 = when
        .map(|dt| to_instant(dt, utc))
        .unwrap_or_else(now_instant);

    println!(
        "地点 {:.4}°{}, {:.4}°{}, 海拔 {:.0} m",
        lat.abs(),
        if lat >= 0.0 { "N" } else { "S" },
        lon.abs(),
        if lon >= 0.0 { "E" } else { "W" },
        alt
    );

    if let Some(n) = days {
        println!("\n{:<22} {:>7} {:>7} {:>7} {:>9} {:>9}", "北京时间", "月龄", "照亮", "相位角", "地心距", "视直径");
        for d in 0..n {
            let t = Instant::from_tt(t0.tt_jd() + d as f64);
            let m = observe(t, &observer, atm);
            println!(
                "{:<22} {:>7.2} {:>6.1}% {:>7.2}° {:>7.0} km {:>7.2}′",
                fmt_beijing(t),
                m.phase.age_days,
                m.phase.illuminated * 100.0,
                m.phase.phase_angle_deg,
                m.phase.distance_km,
                m.phase.angular_diameter_deg * 60.0
            );
        }
        return;
    }

    print_moment(t0, &observer, atm);
}

fn print_moment(t: Instant, observer: &Observer, atm: Option<Atmosphere>) {
    let m = observe(t, observer, atm);
    let p = &m.phase;
    println!("时刻     {}（北京时间）", fmt_beijing(t));
    println!(
        "地心黄道 (via 赤道)  赤经 {:.4}°  赤纬 {:+.4}°  距离 {:.0} km",
        m.equatorial.ra_deg, m.equatorial.dec_deg, m.equatorial.distance_km
    );
    println!(
        "地平    方位角 {:.3}°  高度 {:+.3}°{}",
        m.horizontal.azimuth_deg,
        m.horizontal.altitude_deg,
        if atm.is_some() { "（含折射）" } else { "（几何）" }
    );
    println!(
        "相位    月龄 {:.2} 天  照亮 {:.2}%  相位角 {:.2}°  {}",
        p.age_days,
        p.illuminated * 100.0,
        p.phase_angle_deg,
        if p.waxing { "盈" } else { "亏" }
    );
    println!(
        "形状    视直径 {:.2}′  亮边方位 {:.1}°（自北向东）",
        p.angular_diameter_deg * 60.0,
        p.bright_limb_deg
    );

    let rs = rise_set(t, observer);
    let f = |x: Option<Instant>| x.map(fmt_beijing).unwrap_or_else(|| "—".to_string());
    println!(
        "出没    升 {}   中天 {}   落 {}",
        f(rs.rise),
        f(rs.transit),
        f(rs.set)
    );
}

/// 把 `Instant` 按北京时间写成 `YYYY-MM-DD HH:MM:SS`。
fn fmt_beijing(t: Instant) -> String {
    let jd_bt = t.ut1_jd() + 8.0 / 24.0;
    let mut jdn = (jd_bt + 0.5).floor() as i64;
    let frac = jd_bt + 0.5 - jdn as f64;
    let mut secs = (frac * 86400.0).round() as i64;
    // 四舍五入可能正好进到 86400，要进位到次日，否则会打印出 "24:00:00"
    if secs >= 86400 {
        secs -= 86400;
        jdn += 1;
    }
    let (y, m, d) = ephemeris::ymd_from_jdn(jdn, Calendar::Gregorian);
    format!(
        "{y:04}-{m:02}-{d:02} {:02}:{:02}:{:02}",
        secs / 3600,
        (secs / 60) % 60,
        secs % 60
    )
}

fn parse_datetime(s: &str) -> DateTime {
    let (date, time) = match s.split_once(['T', ' ']) {
        Some((a, b)) => (a, Some(b)),
        None => (s, None),
    };
    let mut it = date.split('-');
    let y: i64 = if date.starts_with('-') {
        it.next();
        -it.next().and_then(|x| x.parse::<i64>().ok()).unwrap_or(2000)
    } else {
        it.next().and_then(|x| x.parse().ok()).unwrap_or(2000)
    };
    let m: u32 = it.next().and_then(|x| x.parse().ok()).unwrap_or(1);
    let d: u32 = it.next().and_then(|x| x.parse().ok()).unwrap_or(1);
    let (h, mi, se) = match time {
        None => (0, 0, 0.0),
        Some(t) => {
            let p: Vec<&str> = t.split(':').collect();
            (
                p.first().and_then(|x| x.parse().ok()).unwrap_or(0),
                p.get(1).and_then(|x| x.parse().ok()).unwrap_or(0),
                p.get(2).and_then(|x| x.parse().ok()).unwrap_or(0.0),
            )
        }
    };
    DateTime::new(y, m, d, h, mi, se)
}

/// 把日期时刻转成 [`Instant`]。`utc` 为真时按 UTC 解释，否则按北京时间。
fn to_instant(dt: DateTime, utc: bool) -> Instant {
    let jd = dt.jd(Calendar::Gregorian);
    let jd_ut1 = if utc { jd } else { jd - 8.0 / 24.0 };
    let dt_sec = ephemeris::delta_t_seconds(2000.0 + (jd_ut1 - 2451545.0) / 365.25);
    Instant::from_tt(jd_ut1 + dt_sec / 86400.0)
}

/// 当前时刻。取系统时间（UTC），换算到 `Instant`。
fn now_instant() -> Instant {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs_f64())
        .unwrap_or(0.0);
    let jd_utc = 2440587.5 + secs / 86400.0;
    Instant::from_tt(jd_utc + ephemeris::delta_t_seconds(2000.0 + (jd_utc - 2451545.0) / 365.25) / 86400.0)
}

fn print_help() {
    println!(
        "moon —— 月相与月球位置\n\
         \n\
         用法:\n\
         \x20 moon [日期][T时刻]            打印该时刻的月相（默认当前、北京）\n\
         \x20 moon --days <n>               打印连续 n 天的月相表\n\
         \n\
         选项:\n\
         \x20 --lat <度>  --lon <度>  --alt <米>   观测地点（北纬/东经为正）\n\
         \x20 --utc                                输入按 UTC 解释\n\
         \x20 --no-refraction                      地平高度不给大气折射\n\
         \n\
         示例:\n\
         \x20 moon 2026-02-17T20:00:00\n\
         \x20 moon --lat 51.5 --lon -0.13 --alt 11 2026-06-01T22:00\n\
         \x20 moon --days 30"
    );
}
