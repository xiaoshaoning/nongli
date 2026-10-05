//! 太阳位置的验收数据源，配 `tools/suncheck.py`。
//!
//! 输出格式（每行空白分隔）：`SUN <年> <视黄经/度> <日地距离/AU>`，
//! 年是从 −4000 到 4000 的等间隔抽样，对应 JDE = 2451545 + (年−2000)·365.25。
//!
//! 这个例子是 `tools/suncheck.py` 的输入。之所以要它：太阳此前只有 `daydiff.py`
//! 这种**间接**验收（看节气落在哪天），于是"我们的 VSOP87D 与 astropy 用的
//! VSOP2000 差 0.7″"这件事一直没被单独看见——它被误当成截断误差了。

use ephemeris::{sun_apparent_longitude, sun_distance_au, Instant};

fn main() {
    const J2000: f64 = 2451545.0;
    for y in (-4000..=4000).step_by(200) {
        let jd = J2000 + (y as f64 - 2000.0) * 365.25;
        let t = Instant::from_tt(jd);
        println!(
            "SUN {y} {:.9} {:.9}",
            sun_apparent_longitude(t),
            sun_distance_au(t),
        );
    }
}
