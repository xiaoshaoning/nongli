//! 对 JPL 星历内核的验收数据源，配 `tools/truthcheck.py`。
//!
//! 输出两种行（空白分隔，JD 为 TT）：
//! `SUN  <jd> <视黄经/度> <日地距离/AU>`
//! `MOON <jd> <地心黄经/度> <地心黄纬/度> <距离/km>`
//!
//! 时间跨度取 1900-01-01 起每 5 天到 2050 年——正好是 DE421 覆盖的范围。
//!
//! 为什么要有这个：`*check.py` 那些工具比的是 ERFA / astropy 的**模型**
//! （VSOP2000、moon98 等），它们各自也有历元约定与截断误差，会把我们的
//! 真误差掩掉。JPL 内核是**真值**，`tools/daydiff.py` 那种"看节气落在哪天"
//! 的间接验收也换不来它。

use ephemeris::{moon_geocentric, sun_apparent_longitude, sun_distance_au, Instant};

fn main() {
    let mut jd = 2415020.5;
    while jd < 2469807.5 {
        let t = Instant::from_tt(jd);
        let m = moon_geocentric(t);
        println!("SUN  {jd:.5} {:.9} {:.9}", sun_apparent_longitude(t), sun_distance_au(t));
        println!("MOON {jd:.5} {:.9} {:.9} {:.9}", m.lon_deg, m.lat_deg, m.distance_km);
        jd += 5.0;
    }
}
