//! 把太阳/月球位置、朔与节气时刻按机器可读格式输出，供 tools/crosscheck.py
//! 与 ERFA / astropy 参考实现比对。
//!
//! ```text
//! cargo run --release --example ephemeris_dump -- -1975 6025 [nm_step] [term_step]
//! ```

use nongli::astro::{moon_apparent_longitude, sun_apparent_longitude};
use nongli::{new_moon_index_near, new_moon_jde, solar_term_index_near, solar_term_jde};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let y0: i64 = args.first().and_then(|s| s.parse().ok()).unwrap_or(-1975);
    let y1: i64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(6025);
    let nm_step: i64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(97);
    let term_step: i64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(61);

    // 1) 以约 2 年一步采样整整两个世纪以上的经度
    let jd0 = 2451545.0 + (y0 as f64 - 2000.0) * 365.25;
    let jd1 = 2451545.0 + (y1 as f64 - 2000.0) * 365.25;
    let n = 4000;
    for i in 0..=n {
        let jd = jd0 + (jd1 - jd0) * i as f64 / n as f64;
        println!(
            "EPH {jd:.6} {:.10} {:.10}",
            sun_apparent_longitude(jd),
            moon_apparent_longitude(jd)
        );
    }

    // 2) 每 `nm_step` 个朔取一个采样点（覆盖整个范围且不穷举）
    let k0 = new_moon_index_near(jd0) - 1;
    let k1 = new_moon_index_near(jd1) + 1;
    let mut k = k0;
    while k <= k1 {
        println!("NM {k} {:.8}", new_moon_jde(k));
        k += nm_step;
    }

    // 3) 节气采样
    let j0 = solar_term_index_near(jd0) - 1;
    let j1 = solar_term_index_near(jd1) + 1;
    let mut j = j0;
    while j <= j1 {
        println!("ST {j} {:.8}", solar_term_jde(j));
        j += term_step;
    }
}
