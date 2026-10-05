//! 导出 ΔT 采样点，供 `probe/dtreference.py --check` 核对 Python 侧的参考实现
//! 与 `src/jd.rs` 是否一致。
//!
//! 两者是跨语言的同一份知识（IERS 实测表 + E-M 多项式），必须能机器校验，
//! 否则改了一边忘了另一边会静默地让 docs/ambiguity.md 的结论失效。
//!
//! ```text
//! cargo run --release --example dt_dump
//! ```

use nongli::delta_t_seconds;

fn main() {
    // 覆盖 ΔT 的每个分支边界，以及实测表的端点。
    let mut years = Vec::new();
    let mut y = -2000i32;
    while y <= 6100 {
        years.push(y as f64);
        y += 37;
    }
    years.extend([
        -500.0, -499.0, 499.0, 500.0, 1599.0, 1600.0, 1699.0, 1700.0, 1799.0, 1800.0, 1859.0,
        1860.0, 1899.0, 1900.0, 1919.0, 1920.0, 1940.0, 1941.0, 1960.0, 1961.0, 1972.0, 1973.0,
        1973.5, 1985.0, 1986.0, 2004.0, 2005.0, 2025.0, 2025.5, 2049.0, 2050.0, 2149.0, 2150.0,
    ]);
    for y in years {
        println!("DT {y} {:.6}", delta_t_seconds(y));
    }
    // 儒略日入口也要一致（平闰/负年无关，这里只查时间尺度换算）
    for jd in [2451545.0f64, 2440000.0, 990000.0, 3900000.0] {
        println!("JD {jd} {:.6}", delta_t_seconds(2000.0 + (jd - 2451545.0) / 365.25));
    }
}
