//! 把月相画到终端上。
//!
//! 只画**明暗分界（晨昏线）**，不画月面纹理——纹理要额外知道月球自转轴在天上的
//! 方位角，那是另一样东西（见 `docs/plan.md` 第 6 步末尾）。
//!
//! # 怎么才算"画得对"
//!
//! 月球被太阳照亮的那一半，投影到天上是个**椭圆**的大圆（晨昏线是球面上的大圆）。
//! 设 `i` 为相位角，沿"亮边方向"取坐标 `u`、垂直方向取 `v`，则圆面上一点 (u, v)
//! 被照亮当且仅当
//!
//! ```text
//! u ≥ −cos(i) · √(1 − v²)
//! ```
//!
//! 三个极端一验就对：`i = 0`（望）时右端是 −√(1−v²)，即整个圆面都被照亮；
//! `i = 180°`（朔）时是 +√(1−v²)，只剩右边缘一条线，等于全暗；
//! `i = 90°`（弦）时右端为 0，正好一半。中间就是月牙与凸月。
//!
//! 亮边方向由 [`MoonPhase::bright_limb_deg`] 给出（自北起向东为正），
//! 屏幕取"北在上、东在左"的天文习惯，故该方向的单位向量为 `(−sin χ, cos χ)`。

use crate::MoonPhase;

/// 背景（天空）。
const SKY: char = ' ';
/// 月面未被照亮的一侧。比天空亮一点点，好看出圆盘的轮廓（真实里是地照）。
const DARK: char = '░';
/// 被照亮的一侧。
const LIT: char = '█';

/// 把月相画成一个字符网格。
///
/// `cols` 是横向字符数；纵向行数取一半——终端字符高约为宽的两倍，这样月亮才是圆的。
/// 返回的行**不含**换行，行尾也没有填充空格。
pub fn render(phase: &MoonPhase, cols: usize) -> Vec<String> {
    let cols = cols.max(8);
    let rows = (cols / 2).max(4);
    // x 与 y 要**各自**归一到半径 1：网格是 cols × cols/2 个字符，而终端字符
    // 高约为宽的两倍，所以屏幕上才是圆的。早先两边共用同一个 radius，
    // 结果把圆盘压成了半高——是"画出的面积不等于照亮比例"这个测试抓到的。
    let rx = cols as f64 / 2.0 - 1.0;
    let ry = rows as f64 / 2.0 - 1.0;
    let cx = (cols as f64 - 1.0) / 2.0;
    let cy = (rows as f64 - 1.0) / 2.0;

    let cos_i = phase.phase_angle_deg.to_radians().cos();
    let chi = phase.bright_limb_deg.to_radians();
    let (s_chi, c_chi) = chi.sin_cos();

    let mut out = Vec::with_capacity(rows);
    for j in 0..rows {
        let mut line = String::with_capacity(cols);
        for i in 0..cols {
            // 归一到"圆盘半径 = 1"，x 正为屏幕右，y 正为上
            let x = (i as f64 - cx) / rx;
            let y = (cy - j as f64) / ry;
            let rr = x * x + y * y;
            if rr > 1.0 {
                line.push(SKY);
                continue;
            }
            // u 沿亮边方向，v 与之垂直
            let u = -x * s_chi + y * c_chi;
            let v = x * c_chi + y * s_chi;
            let terminator = -cos_i * (1.0 - v * v).max(0.0).sqrt();
            line.push(if u >= terminator { LIT } else { DARK });
        }
        out.push(line.trim_end().to_string());
    }
    out
}

/// 网格里被照亮像素占圆盘像素的比例。**只给测试用**——它是验证手段，不是接口。
#[cfg(test)]
fn lit_fraction(lines: &[String]) -> f64 {
    let (mut lit, mut disk) = (0usize, 0usize);
    for l in lines {
        for c in l.chars() {
            match c {
                LIT => {
                    lit += 1;
                    disk += 1;
                }
                DARK => disk += 1,
                _ => {}
            }
        }
    }
    if disk == 0 {
        0.0
    } else {
        lit as f64 / disk as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MoonPhase;

    fn phase(illuminated: f64, bright_limb_deg: f64) -> MoonPhase {
        // phase_angle 由照亮比例反推：k = (1+cos i)/2
        let i = (2.0 * illuminated - 1.0).clamp(-1.0, 1.0).acos().to_degrees();
        MoonPhase {
            phase_angle_deg: i,
            illuminated,
            age_days: 0.0,
            waxing: true,
            bright_limb_deg,
            distance_km: 384_400.0,
            angular_diameter_deg: 0.0,
            libration: ephemeris::Libration { lon_deg: 0.0, lat_deg: 0.0 },
        }
    }

    /// 渲染出来的被照亮像素比例必须与照亮比例相符。
    ///
    /// 这正是验收里"朔/望/上下弦肉眼与物理预期一致"的可机器检验版本——
    /// 比肉眼看图可靠，而且能进 CI。
    #[test]
    fn rendered_area_matches_illuminated_fraction() {
        for k in [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0] {
            for chi in [0.0, 45.0, 90.0, 180.0, 270.0] {
                let lines = render(&phase(k, chi), 80);
                let got = lit_fraction(&lines);
                // 像素化误差：圆周长 / 面积 量级，80 列时约百分之几
                assert!(
                    (got - k).abs() < 0.04,
                    "k={k} χ={chi}: 画出来 {got:.4}"
                );
            }
        }
    }

    /// 四个几何点的直观形状。
    #[test]
    fn four_shapes() {
        let full = render(&phase(1.0, 0.0), 40);
        assert!(full.iter().all(|l| !l.contains(DARK)), "望应全亮");
        let new = render(&phase(0.0, 0.0), 40);
        assert!(new.iter().all(|l| !l.contains(LIT)), "朔应全暗");
        // 上下弦：恰好一半
        for chi in [90.0, 270.0] {
            let q = lit_fraction(&render(&phase(0.5, chi), 80));
            assert!((q - 0.5).abs() < 0.03, "弦 χ={chi}: {q}");
        }
    }

    /// 亮边方向必须真的跟着 χ 走：亮的一侧应在 χ 所指的那半边。
    #[test]
    fn bright_side_follows_position_angle() {
        // χ=90°（东，屏幕上为左）时，亮的是左半边
        let lines = render(&phase(0.5, 90.0), 80);
        let w = lines.iter().map(|l| l.chars().count()).max().unwrap();
        let (mut left, mut right) = (0, 0);
        for l in &lines {
            for (i, c) in l.chars().enumerate() {
                if c == LIT {
                    if i < w / 2 {
                        left += 1;
                    } else {
                        right += 1;
                    }
                }
            }
        }
        assert!(left > right * 5, "χ=90° 应亮在左：左 {left} 右 {right}");

        // χ=270°（西，屏幕上为右）时反过来
        let lines = render(&phase(0.5, 270.0), 80);
        let (mut left, mut right) = (0, 0);
        for l in &lines {
            for (i, c) in l.chars().enumerate() {
                if c == LIT {
                    if i < w / 2 {
                        left += 1;
                    } else {
                        right += 1;
                    }
                }
            }
        }
        assert!(right > left * 5, "χ=270° 应亮在右：左 {left} 右 {right}");
    }

    /// 圆盘占满网格、而且是圆的。
    ///
    /// 只比"字符宽高比 2:1"是测不出问题的——早先把圆压成半高时该比值照样是 2。
    /// 真正管用的是：最宽行接近 cols、非空行数接近 rows。
    #[test]
    fn disk_fills_the_grid_and_is_round() {
        let cols = 60;
        let lines = render(&phase(1.0, 0.0), cols);
        let nonempty = lines.iter().filter(|l| !l.trim().is_empty()).count();
        let widest = lines.iter().map(|l| l.chars().count()).max().unwrap();
        let rows = cols / 2;
        assert!(
            nonempty >= rows - 2,
            "非空行数 {nonempty} 应接近 {rows}（圆盘被压扁了？）"
        );
        assert!(
            widest >= cols - 3,
            "最宽行 {widest} 应接近 {cols}"
        );
    }
}
