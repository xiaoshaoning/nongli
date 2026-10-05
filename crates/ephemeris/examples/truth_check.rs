//! 对 JPL 星历内核的验收数据源，配 `tools/truthcheck.py`。
//!
//! 输出两种行（空白分隔，JD 为 TT）：
//! `SUN  <jd> <视黄经/度> <几何黄经/度> <日地距离/AU>`
//! `MOON <jd> <地心黄经/度> <地心黄纬/度> <距离/km>`
//!
//! 时间跨度取 1900-01-01 起每 5 天到 2050 年——正好是 DE421 覆盖的范围。
//!
//! **几何黄经那一列是关键**：它不含光行差与章动，所以可以直接与
//! "DE421 的地心矢量 → `erfa.ecm06`" 比，验收不依赖我们自己的任何模型。
//! 视黄经那一列要在真值侧补上光行差与章动才能比，等于复刻了我们的口径，
//! 只能当辅助。

use ephemeris::{
    moon_geocentric, sun_apparent_longitude, sun_distance_au, sun_geometric_longitude, Instant,
};

fn main() {
    let mut jd = 2415020.5;
    while jd < 2469807.5 {
        let t = Instant::from_tt(jd);
        let m = moon_geocentric(t);
        println!(
            "SUN  {jd:.5} {:.9} {:.9} {:.9}",
            sun_apparent_longitude(t),
            sun_geometric_longitude(t),
            sun_distance_au(t),
        );
        println!("MOON {jd:.5} {:.9} {:.9} {:.9}", m.lon_deg, m.lat_deg, m.distance_km);
        jd += 5.0;
    }
}
