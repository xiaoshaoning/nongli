//! 导出参考系量（章动、黄赤交角、黄道→赤道），供 `tools/framecheck.py` 对 ERFA 比对。
//!
//! ```text
//! cargo run --release -p ephemeris --example frame_check -- -4000 4000
//! ```

use ephemeris::{mean_obliquity, nutation, true_obliquity, Ecliptic, Instant};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let y0: i64 = args.first().and_then(|s| s.parse().ok()).unwrap_or(-4000);
    let y1: i64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(4000);
    let n = 4000;
    for i in 0..=n {
        let y = y0 as f64 + (y1 - y0) as f64 * i as f64 / n as f64;
        let t = Instant::from_tt(2451545.0 + (y - 2000.0) * 365.25);
        let nut = nutation(t);
        println!(
            "NUT {y:.6} {:.9} {:.9} {:.9} {:.9}",
            nut.dpsi_deg * 3600.0,
            nut.deps_deg * 3600.0,
            mean_obliquity(t),
            true_obliquity(t)
        );
    }
    // 黄道→赤道：固定一组 (λ, β)，看 α、δ
    let t = Instant::from_tt(2451545.0);
    for (lam, bet) in [
        (0.0, 0.0),
        (90.0, 0.0),
        (45.0, 30.0),
        (200.0, -60.0),
        (359.0, 5.0),
    ] {
        let q = Ecliptic {
            lon_deg: lam,
            lat_deg: bet,
            distance_km: 1.0,
        }
        .equatorial(t);
        println!("EQ {lam} {bet} {:.12} {:.12}", q.ra_deg, q.dec_deg);
    }
}
