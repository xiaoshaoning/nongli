//! 角度运算与角度-时刻求根。
//!
//! [`norm360`] / [`norm180`] 是对外接口（任何比较角度的调用方都需要）；
//! [`newton`] 是 crate 内部的时间求根器。所有角度均为**度**。

/// 度 → 弧度。
pub(crate) const D2R: f64 = core::f64::consts::PI / 180.0;
/// 弧度 → 度。
pub const R2D: f64 = 180.0 / core::f64::consts::PI;

/// 归化到 [0, 360)。
#[inline]
pub(crate) fn norm360(x: f64) -> f64 {
    x.rem_euclid(360.0)
}

/// 归化到 (−180, 180]。
#[inline]
pub(crate) fn norm180(x: f64) -> f64 {
    let y = x.rem_euclid(360.0);
    if y > 180.0 {
        y - 360.0
    } else {
        y
    }
}

/// 角度求根的收敛判据，单位**度**。
///
/// 它必须**大于** `jde ≈ 2.45e6` 处 f64 能达到的残差地板：该处相邻可表示数相距
/// 4.7e-10 天，而月球每天走 12.19°，所以残差低于约 6e-9° 后，任何步长都无法再
/// 改变 `t`。取 1e-8° 相当于月球 2e-8 秒、太阳 2e-7 秒，远小于任何有意义的量
/// （农历日界以秒计）。
///
/// 定得过小（曾用 1e-10°）会让循环直到上限都退不出，白算 35 次。
const ANGLE_TOLERANCE_DEG: f64 = 1e-8;

/// 牛顿迭代上限。收敛是二次的，实测 3–4 次即到底；留 12 次以防病态初值。
const MAX_ITERATIONS: u32 = 12;

/// 求使 `f(t) = 0` 的时刻。
///
/// `f` 必须返回一个归化到 (−180, 180] 的**角度差**，于是它过零时斜率约为
/// 该天体每日转过的度数。`h` 是数值微分的半步长（日），`max_step` 限幅以免发散。
pub(crate) fn newton<F: Fn(f64) -> f64>(mut t: f64, f: F, h: f64, max_step: f64) -> f64 {
    for _ in 0..MAX_ITERATIONS {
        let y = f(t);
        if y.abs() < ANGLE_TOLERANCE_DEG {
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
        // 步长小于 `t` 自身的表示精度时，再迭代只是空转。这是"已到 f64 极限"的
        // 唯一可靠信号；光靠调小容差解决不了（容差不可能小于表示精度）。
        if t - step == t {
            break;
        }
        t -= step;
    }
    t
}
