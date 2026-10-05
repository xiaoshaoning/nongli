//! 角度运算与角度-时刻求根。
//!
//! [`norm360`] / [`norm180`] 是对外接口（任何比较角度的调用方都需要）；
//! [`newton`] 是 crate 内部的时间求根器。所有角度均为**度**。

/// 度 → 弧度。
pub const D2R: f64 = core::f64::consts::PI / 180.0;
/// 弧度 → 度。
pub const R2D: f64 = 180.0 / core::f64::consts::PI;

/// 归化到 [0, 360)。
#[inline]
pub fn norm360(x: f64) -> f64 {
    x.rem_euclid(360.0)
}

/// 归化到 (−180, 180]。
#[inline]
pub fn norm180(x: f64) -> f64 {
    let y = x.rem_euclid(360.0);
    if y > 180.0 {
        y - 360.0
    } else {
        y
    }
}

/// 求使 `f(t) = 0` 的时刻。
///
/// `f` 必须返回一个归化到 (−180, 180] 的**角度差**，于是它过零时斜率约为
/// 该天体每日转过的度数。`h` 是数值微分的半步长（日），`max_step` 限幅以免发散。
///
/// 迭代上限 40 次、收敛判据 1e-10 度（约 1e-4 秒），对朔与节气都足够。
pub(crate) fn newton<F: Fn(f64) -> f64>(mut t: f64, f: F, h: f64, max_step: f64) -> f64 {
    for _ in 0..40 {
        let y = f(t);
        if y.abs() < 1e-10 {
            break;
        }
        let dy = (f(t + h) - f(t - h)) / (2.0 * h);
        if dy.abs() < 1e-12 {
            break;
        }
        let mut step = y / dy;
        if step > max_step {
            step = max_step;
        } else if step < -max_step {
            step = -max_step;
        }
        t -= step;
    }
    t
}
