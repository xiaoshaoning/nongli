//! 导出月球的地心黄道位置 (λ, β, Δ)，供 `tools/mooncheck.py` 对 ERFA / astropy 比对。
//!
//! ```text
//! cargo run --release -p ephemeris --example moon_check -- -1975 6025
//! ```
//!
//! 与 `ephemeris_dump` 分开：那个导出的是**求根结果**（朔、节气），
//! 这个导出的是**理论输出**（位置）。

use ephemeris::{moon_geocentric, Instant};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let y0: i64 = args.first().and_then(|s| s.parse().ok()).unwrap_or(-1975);
    let y1: i64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(6025);
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4000);

    let jd0 = 2451545.0 + (y0 as f64 - 2000.0) * 365.25;
    let jd1 = 2451545.0 + (y1 as f64 - 2000.0) * 365.25;
    for i in 0..=n {
        let jd = jd0 + (jd1 - jd0) * i as f64 / n as f64;
        let e = moon_geocentric(Instant::from_tt(jd));
        println!("{jd:.6} {:.10} {:.10} {:.6}", e.lon_deg, e.lat_deg, e.distance_km);
    }
}
