//! GB/T 33661-2017 第 **5.2** 条（"朔和节气的北京时间计算精度应达到 1 s"）的验收数据源。
//!
//! 配 `tools/spec_check.py`。输出两种行：
//!
//! ```text
//! NEWMOON <k> <JD_TT>
//! TERM    <j> <JD_TT>
//! ```
//!
//! 给的是事件时刻本身，**不换算北京时间**——两边用同一套 ΔT，所以作差得到的
//! 正是 5.2 要量的"模型精度"（ΔT 自身的不确定度是另一回事，见 README 那张表）。
//!
//! 范围取 1900-01-01 … 2050-01-01，因为 JPL DE421（真值）只覆盖 1900–2053。

use ephemeris::{new_moon, new_moon_index_near, Instant};

fn main() {
    const LO: f64 = 2415020.5;
    const HI: f64 = 2469807.5;

    // 朔
    let mut k = new_moon_index_near(Instant::from_tt(LO)) - 2;
    loop {
        let t = new_moon(k);
        if t.tt_jd() > HI {
            break;
        }
        if t.tt_jd() >= LO {
            println!("NEWMOON {k} {:.9}", t.tt_jd());
        }
        k += 1;
    }

    // 节气：序号 j 定义为太阳视黄经 = 15j (mod 360)
    let mut j = nongli::solar_term_index_near(Instant::from_tt(LO)) - 2;
    loop {
        let t = nongli::solar_term(j);
        if t.tt_jd() > HI {
            break;
        }
        if t.tt_jd() >= LO {
            println!("TERM {j} {:.9}", t.tt_jd());
        }
        j += 1;
    }
}
