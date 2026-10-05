//! 导出恒星时与地平坐标，供 `tools/observercheck.py` 对 ERFA / astropy 比对。
//!
//! ```text
//! cargo run --release -p ephemeris --example observer_check
//! ```

use ephemeris::{gast, gmst, Atmosphere, Ecliptic, Equatorial, Instant, Observer};

/// 覆盖南北半球、赤道、高纬、负海拔。
const SITES: &[(f64, f64, f64, &str)] = &[
    (39.9042, 116.4074, 44.0, "Beijing"),
    (0.0, 0.0, 0.0, "equator"),
    (-33.8688, 151.2093, 58.0, "Sydney"),
    (-77.8463, 166.6682, 30.0, "McMurdo"),
    (64.1355, -21.8954, 20.0, "Reykjavik"),
    (28.0, -16.6, 2400.0, "Teide"),
    (31.5, 35.5, -400.0, "DeadSea"),
];

fn main() {
    // 1) 恒星时：±4000 年内 400 点
    for i in 0..=400 {
        let y = -4000.0 + 8000.0 * i as f64 / 400.0;
        let t = Instant::from_tt(2451545.0 + (y - 2000.0) * 365.25);
        println!("GST {y:.6} {:.12} {:.12}", gmst(t), gast(t));
    }

    // 2) 地平坐标：对若干历元 × 若干站点，取一组固定黄道坐标
    let ecls = [
        (0.0, 0.0, 384400.0),
        (90.0, 20.0, 384400.0),
        (200.0, -60.0, 384400.0),
        (359.0, 5.0, 384400.0),
    ];
    for y in [1900.0, 2000.0, 2026.0, 2100.0, 3000.0] {
        let t = Instant::from_tt(2451545.0 + (y - 2000.0) * 365.25);
        for &(lat, lon, h, _name) in SITES {
            let o = Observer {
                lat_deg: lat,
                lon_deg: lon,
                height_m: h,
            };
            for (lam, bet, dist) in ecls {
                let eq: Equatorial = Ecliptic {
                    lon_deg: lam,
                    lat_deg: bet,
                    distance_km: dist,
                }
                .equatorial(t);
                let vacuum = eq.horizontal(t, &o, None);
                let refr = eq.horizontal(t, &o, Some(Atmosphere::default()));
                println!(
                    "HZ {y} {lat} {lon} {h} {lam} {bet} {:.12} {:.12} {:.12} {:.12}",
                    vacuum.azimuth_deg, vacuum.altitude_deg, refr.azimuth_deg, refr.altitude_deg
                );
            }
        }
    }
}
