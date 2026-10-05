//! 导出月相量，供 `tools/moonphasecheck.py` 对 ERFA 独立比对。
//!
//! ```text
//! cargo run --release -p moon --example phase_check
//! ```

use ephemeris::{Atmosphere, Instant, Observer};
use moon::observe;

fn main() {
    let o = Observer {
        lat_deg: 39.9042,
        lon_deg: 116.4074,
        height_m: 44.0,
    };
    // 覆盖约 20 年，步长 1 天，足以扫过各种相位与距离
    for i in 0..7300 {
        let t = Instant::from_tt(2451545.0 + i as f64);
        let m = observe(t, &o, Some(Atmosphere::default()));
        let p = m.phase;
        println!(
            "{:.6} {:.12} {:.9} {:.9} {:.6} {:.9} {:.9}",
            t.tt_jd(),
            p.illuminated,
            p.phase_angle_deg,
            p.bright_limb_deg,
            p.distance_km,
            p.angular_diameter_deg,
            m.horizontal.altitude_deg
        );
    }
}
