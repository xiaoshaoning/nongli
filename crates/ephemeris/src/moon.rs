//! 月球理论：**ELP2000-82B**（Chapront-Touzé, Chapront & Francou 1985；
//! 2002 年按 LLR 重定轨道参数）。
//!
//! **这里是整个 crate 的理论之"缝"。** 系数取自 `third_party/elp2000-82b/`
//! （MIT），按作者标注截断到 0.001″ / 0.001 km，36 张表共 3402 项。
//! 换更细的截断、或换 ELP/MPP02，都只需替换 [`moon_geocentric`]。
//!
//! 此前用的是 Meeus《Astronomical Algorithms》第 47 章那 60 项截断级数
//! （对 ELP/MPP02 的标称精度 RMS 2.9″、最坏 18.3″）。换过来的实测收益见
//! `docs/plan.md`，量它的工具是 `tools/truthcheck.py`（对 JPL 真值）。

use crate::angle::{norm180, norm360, newton, D2R, R2D};
// 通配导入：这里 36 张表 + 类型 + 常数全部都用得上，逐个列反而更难核对。
use crate::elp2000_tables::*;
use crate::frames::{precession_in_longitude, Ecliptic};
use crate::sun::sun_geometric_longitude;
use crate::time::{Instant, J2000};

/// 平朔望月，日。
const SYNODIC_MONTH_DAYS: f64 = 29.530588861;
/// J2000 附近的一次朔，作为 [`new_moon`] 序号 `k = 0` 的参考。
const EPOCH_NEW_MOON_JDE: f64 = 2451550.09766;
/// 角秒 → 弧度。ELP2000 的引数与振幅都以角秒给出（距离以 km）。
const ARCSEC: f64 = core::f64::consts::PI / 648000.0;

/// W1、W2、W3、T、ϖ′ 的多项式（弧秒），取到 `deg` 阶。
///
/// ELP2000 内部两处用**不同阶数**：主问题（ELP01-03）用完整的 5 阶，
/// 其余各表只用常数与一次项。这是理论本身的分工（扰动项的频率不需要那么准），
/// 不是我们做的近似。
fn argument_polynomials(t: f64, deg: usize) -> [f64; 5] {
    let mut out = [0.0; 5];
    for (o, c) in out.iter_mut().zip(ELP_ARG_POLY.iter()) {
        let (mut tn, mut v) = (1.0, 0.0);
        for &ck in c.iter().take(deg) {
            v += ck * tn;
            tn *= t;
        }
        *o = v;
    }
    out
}

/// 四个 Delaunay 引数 `(D, l′, l, F)`，弧秒。
fn delaunay(p: &[f64; 5]) -> [f64; 4] {
    [
        p[0] - p[3] + 648000.0, // D = W1 − T + 648000
        p[3] - p[4],            // l′ = T − ϖ′   太阳平近点角
        p[0] - p[1],            // l  = W1 − W2  月球平近点角
        p[0] - p[2],            // F  = W1 − W3  升交点角距
    ]
}

/// 主问题表里 `A` 的拟合改正量（ELP2000 第 7 节，拟合到 DE200/LE200）。
///
/// `nu`、`np` 是**速率**的相对改正（所以除以 W1 的一次项），其余是角量。
struct Corrections {
    nu: f64,
    np: f64,
    gamma: f64,
    e: f64,
    ep: f64,
}

impl Corrections {
    fn new() -> Self {
        let n = ELP_ARG_POLY[0][1];
        Self {
            nu: ELP_DELTANU_ARCSEC * ARCSEC / n,
            np: ELP_DELTANP_ARCSEC * ARCSEC / n,
            gamma: ELP_DELTAGAMMA_ARCSEC * ARCSEC,
            e: ELP_DELTAE_ARCSEC * ARCSEC,
            ep: ELP_DELTAEP_ARCSEC * ARCSEC,
        }
    }
}

/// 主问题（ELP01-03）：`Σ A′·sin(i1·D + i2·l′ + i3·l + i4·F)`。
///
/// `cosine` 为真时改算 `cos`——ELP03 给的是距离，用的是余弦级数。
/// `A′` 是在 `A` 上加按 DE200/LE200 拟合的改正；距离那条还多一个 `−2/3·A·δν`。
fn main_problem_sum(terms: &[ElpMainTerm], d: &[f64; 4], cosine: bool, c: &Corrections) -> f64 {
    let mut s = 0.0;
    for &(i1, i2, i3, i4, a, b1, b2, b3, b4, b5, _b6) in terms {
        let mut amp = a
            + (b1 + ELP_DTASM * b5) * (c.np - ELP_AM * c.nu)
            + b2 * c.gamma
            + b3 * c.e
            + b4 * c.ep;
        if cosine {
            amp -= 2.0 / 3.0 * a * c.nu;
        }
        let arg =
            (i1 as f64 * d[0] + i2 as f64 * d[1] + i3 as f64 * d[2] + i4 as f64 * d[3]) * ARCSEC;
        s += amp * if cosine { arg.cos() } else { arg.sin() };
    }
    s
}

/// 非行星扰动（ELP04-09、ELP22-36）：
/// `Σ A·sin(i1·ζ + i2·D + i3·l′ + i4·l + i5·F + φ)`，`φ` 为度。
///
/// 潮汐、月体形状、相对论、太阳偏心率几组（ELP22-36）不含 ζ，调用时传 0。
fn non_planetary_sum(terms: &[ElpTerm], zeta: f64, d: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for &(i1, i2, i3, i4, i5, phi, a) in terms {
        let arg = (i1 as f64 * zeta
            + i2 as f64 * d[0]
            + i3 as f64 * d[1]
            + i4 as f64 * d[2]
            + i5 as f64 * d[3])
            * ARCSEC
            + phi * D2R;
        s += a * arg.sin();
    }
    s
}

/// 行星扰动（ELP10-21）：前七个乘子作用于七颗行星的平黄经，后四个作用于 `tail`。
///
/// 两张表的**结构完全相同**，只是后四项的含义不同，所以一个函数就够——
/// 表 1（ELP10-15）的 `tail` 是 `[海王星, D, l, F]`，表 2（ELP16-21）是 `[D, l′, l, F]`。
/// 这样调用处直接看得出"这张表的后四项是什么"，不必写两个近乎一样的求和。
fn planetary_sum(terms: &[ElpPlanetaryTerm], planets: &[f64; 7], tail: [f64; 4]) -> f64 {
    let mut s = 0.0;
    for &(i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, phi, a) in terms {
        let arg = (i1 as f64 * planets[0]
            + i2 as f64 * planets[1]
            + i3 as f64 * planets[2]
            + i4 as f64 * planets[3]
            + i5 as f64 * planets[4]
            + i6 as f64 * planets[5]
            + i7 as f64 * planets[6]
            + i8 as f64 * tail[0]
            + i9 as f64 * tail[1]
            + i10 as f64 * tail[2]
            + i11 as f64 * tail[3])
            * ARCSEC
            + phi * D2R;
        s += a * arg.sin();
    }
    s
}

/// 月球地心**黄道**位置。
///
/// # 约定（改动会破坏所有下游）
///
/// | 项 | 取值 |
/// |---|---|
/// | 参考点 | **地心** |
/// | 参考面/点 | **平**黄道、**平**分点（*不是*真分点）|
/// | 章动 | **不含**——要真分点用 [`Ecliptic::apparent`] |
/// | 光行时 | **不含**（ELP2000-82B 给的就是时刻 `t` 的几何位置）|
/// | 单位 | 度、度、km；经度归化到 [0,360) |
///
/// 参考面是 **IAU 2006** 的当日平黄道。ELP2000 原始输出的口径不统一
/// （黄经以 J2000 平春分点起算、黄纬已按当日黄道），见函数体的说明。
///
/// 关于光行时：Meeus 第 47 章那份截断级数把月球平黄经的光行时包进去了，
/// 所以 SOFA 的 `eraMoon98` 与它有约 −0.70″ 的系统差。ELP2000-82B 不含，
/// 与 JPL 的几何位置同口径——换用它之后**不该**再扣那 −0.70″。
fn elp_position(tc: f64) -> Ecliptic {
    let p_full = argument_polynomials(tc, 5);
    let p_lin = argument_polynomials(tc, 2);
    let d_full = delaunay(&p_full);
    let d_lin = delaunay(&p_lin);
    // ζ：岁差角，只进 ELP04-09（地球形状摄动）
    let zeta = ELP_ARG_POLY[0][0] + tc * (ELP_ARG_POLY[0][1] + ELP_PRECESSION_CONSTANT);
    // 八颗行星的平黄经
    let planets: [f64; 8] =
        core::array::from_fn(|i| ELP_PLANET_POLY[i][0] + ELP_PLANET_POLY[i][1] * tc);
    // 行星表 1 与 2 的后四项（见 planetary_sum）
    let p7: &[f64; 7] = planets[..7].try_into().unwrap(); // 长度由类型保证
    let tail1 = [planets[7], d_lin[0], d_lin[2], d_lin[3]];
    let tail2 = d_lin; // 表 2 的后四项就是四个 Delaunay 引数本身

    let c = Corrections::new();
    let t2 = tc * tc;

    // 下面三块按 ELP2000 的表结构逐项相加（黄经用第 1 类表、黄纬第 2 类、距离第 0 类），
    // 每项后面的因子是 T 的幂。写法与原始表一一对应，不做归并——这张表就是理论本身。
    let lon_as = main_problem_sum(ELP01, &d_full, false, &c)
        + p_full[0]
        + non_planetary_sum(ELP04, zeta, &d_lin)
        + non_planetary_sum(ELP07, zeta, &d_lin) * tc
        + planetary_sum(ELP10, p7, tail1)
        + planetary_sum(ELP13, p7, tail1) * tc
        + planetary_sum(ELP16, p7, tail2)
        + planetary_sum(ELP19, p7, tail2) * tc
        + non_planetary_sum(ELP22, 0.0, &d_lin)
        + non_planetary_sum(ELP25, 0.0, &d_lin) * tc
        + non_planetary_sum(ELP28, 0.0, &d_lin)
        + non_planetary_sum(ELP31, 0.0, &d_lin)
        + non_planetary_sum(ELP34, 0.0, &d_lin) * t2;

    let lat_as = main_problem_sum(ELP02, &d_full, false, &c)
        + non_planetary_sum(ELP05, zeta, &d_lin)
        + non_planetary_sum(ELP08, zeta, &d_lin) * tc
        + planetary_sum(ELP11, p7, tail1)
        + planetary_sum(ELP14, p7, tail1) * tc
        + planetary_sum(ELP17, p7, tail2)
        + planetary_sum(ELP20, p7, tail2) * tc
        + non_planetary_sum(ELP23, 0.0, &d_lin)
        + non_planetary_sum(ELP26, 0.0, &d_lin) * tc
        + non_planetary_sum(ELP29, 0.0, &d_lin)
        + non_planetary_sum(ELP32, 0.0, &d_lin)
        + non_planetary_sum(ELP35, 0.0, &d_lin) * t2;

    let dist_km = main_problem_sum(ELP03, &d_full, true, &c)
        + non_planetary_sum(ELP06, zeta, &d_lin)
        + non_planetary_sum(ELP09, zeta, &d_lin) * tc
        + planetary_sum(ELP12, p7, tail1)
        + planetary_sum(ELP15, p7, tail1) * tc
        + planetary_sum(ELP18, p7, tail2)
        + planetary_sum(ELP21, p7, tail2) * tc
        + non_planetary_sum(ELP24, 0.0, &d_lin)
        + non_planetary_sum(ELP27, 0.0, &d_lin) * tc
        + non_planetary_sum(ELP30, 0.0, &d_lin)
        + non_planetary_sum(ELP33, 0.0, &d_lin)
        + non_planetary_sum(ELP36, 0.0, &d_lin) * t2;

    // ELP2000-82B 的输出是**混合口径**，这一点参考实现（名字叫 `...OfDate`）没说清：
    //   * 黄经以 **J2000 平春分点** 起算（它的 W1 是恒星黄经，不含岁差）；
    //   * 黄纬**已经**是按当日黄道。
    // 实测（对 JPL DE421，1900–2050）：这样理解后 Δλ 中位 0.05″、Δβ 中位 0.06″；
    // 若当成"整块都是黄道 J2000"再整体做岁差，Δβ 会差到 46″。
    // 所以只补黄经岁差，黄纬不动。
    Ecliptic {
        lon_deg: norm360(lon_as * ARCSEC * R2D),
        lat_deg: lat_as * ARCSEC * R2D,
        distance_km: dist_km,
    }
}

/// 光速，km/s。
const C_KM_S: f64 = 299_792.458;

/// 月球地心**视**位置（用于农历）。
///
/// 在 `elp_position`（私有的纯 ELP2000 求值）的**几何**位置之后补两件观测效应：
///
/// 1. **光行时**：看到的是 τ = Δ/c ≈ 1.28 s 之前的月球。ELP2000-82B 给的是几何位置
///    （Meeus 第 47 章那份截断级数倒是把这一项包进去了），所以换用它之后必须自己补，
///    否则"朔"会漂 ~1.4 s。做法就是**在 t−τ 处重算一遍**——τ 依赖 Δ，而 Δ 依赖位置，
///    所以先算一次拿 Δ，再算一次取结果。
/// 2. **周年光行差**：与太阳同一个量级（≈20.5″）。对"朔"（日月视黄经之差）它会抵消，
///    但绝对黄经要用它才对得上 JPL 的 apparent 位置。
pub fn moon_geocentric(t: Instant) -> Ecliptic {
    let tc = (t.tt_jd() - J2000) / 36525.0;
    let first = elp_position(tc);
    // 光行时（天）：τ = Δ/c
    let tau_days = first.distance_km / C_KM_S / 86400.0;
    let mut e = elp_position(tc - tau_days / 36525.0);
    // 把几何距离留成第一次的那个：光行时对 Δ 的影响只有 ~0.06 km，不值得再算
    e.distance_km = first.distance_km;
    // 周年光行差：**与太阳同一个量、同一个符号**——观测者速度对两个天体是同一个
    // 矢量，所以黄经上的位移也相同（约 −20.5″）。这里曾写成减号，结果月球被
    // 推向太阳的反方向，两者的光行差不但没抵消，还叠成两倍，"朔"整差 ~80 s。
    e.lon_deg = norm360(e.lon_deg + crate::sun::sun_aberration_deg(t));
    // 黄经岁差（ELP 的恒星口径 → 当日平春分点）。放这里而不是 elp_position 里，
    // 是因为它只需要 t：用 t 还是 t−τ 差 3e-9″，不必较真。
    e.lon_deg = norm360(e.lon_deg + precession_in_longitude(t) * R2D);
    e
}

/// 月球地心视黄经 (度，[0,360)，真分点起算)。
///
/// 等价于 `moon_geocentric(t).apparent(t).lon_deg`，留作便捷入口。
pub fn moon_apparent_longitude(t: Instant) -> f64 {
    moon_geocentric(t).apparent(t).lon_deg
}

/// 第 `k` 次朔（日月地心视黄经相等）。
///
/// `k = 0` 对应 2000-01-06 附近；`k` 每 +1 前进一个朔望月。
pub fn new_moon(k: i64) -> Instant {
    syzygy(k, 0.0)
}

/// 第 `k` 次望（日月地心视黄经相差 180°）。
///
/// 与 [`new_moon`] 是**同一个方程**，只是目标差 180°——所以共用下面的求根，
/// 不另写一套。日食判在朔、月食判在望，两者都要用到。
pub fn full_moon(k: i64) -> Instant {
    syzygy(k, 180.0)
}

/// 朔望的一般式：解 `λ_moon − λ_sun_apparent = target`（度）。
///
/// 这里**故意不用**视黄经：条件是日月视黄经之差，而章动 Δψ 对两者是同一个量、
/// 作差时精确抵消。少算两遍 IAU 2000A 章动，根一模一样。
/// 太阳的周年光行差**不能**省——它只作用于太阳。
fn syzygy(k: i64, target_deg: f64) -> Instant {
    let kf = k as f64;
    let t = kf / 1236.85;
    let jde = EPOCH_NEW_MOON_JDE
        + SYNODIC_MONTH_DAYS * kf
        + 0.00015437 * t * t
        - 0.000000150 * t * t * t
        + 0.00000000073 * t * t * t * t;
    // 望的初值要从朔挪半个朔望月，否则牛顿会收敛回同一个朔
    let jde = jde + SYNODIC_MONTH_DAYS * target_deg / 360.0;
    let f = |x: f64| {
        let i = Instant::from_tt(x);
        norm180(
            moon_geocentric(i).lon_deg
                - sun_geometric_longitude(i)
                - crate::sun::sun_aberration_deg(i)
                - target_deg,
        )
    };
    Instant::from_tt(newton(jde, f, 0.05, 2.0))
}

/// 月球**光学天平动**：从地球看到的月面中心（sub-Earth point）的月面经纬度。
///
/// 月球自转基本均匀，但公转不均匀（轨道偏心），且自转轴与轨道面有夹角，于是我们
/// 能多看一点东、西、南、北——这就是天平动。
///
/// | 分量 | 幅度 | 本实现 |
/// |---|---|---|
/// | 经度方向的光学天平动 | ±7.9° | ✅ |
/// | 纬度方向的光学天平动 | ±6.7° | ✅ |
/// | 物理天平动（月球本体真实摆动） | ±0.04° | ❌ |
/// | 周日天平动（观测者不在圆心） | ±1° | 由 `moon` crate 另行处理 |
///
/// **物理天平动未实现**：幅度 ±0.04°，在月面上约合直径的 0.07%，肉眼看不出；
/// 而 Meeus 第 53 章那套级数需要另一张约 20 项的表。等真要拿月面纹理定位再说。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Libration {
    /// 月面经度 l，度，东正。范围约 ±8°。
    pub lon_deg: f64,
    /// 月面纬度 b，度，北正。范围约 ±7°。
    pub lat_deg: f64,
}

/// 由"观测者到月球"的**黄道方向**算月面中心经纬度。
///
/// 私有：目前只有 [`libration`] 一个消费者。将来做**周日天平动**（观测者不在地心）
/// 时会需要它——那时再公开，成本是一行。第 3 步就是在这个模式上栽过
/// （`fukushima_williams`，靠"以后要用"留着，已删）。
///
/// `lon_deg` 用当日**平**分点起算的几何黄经——**不要**加章动。实测：Meeus 例 53.a
/// 的 l′ 是 −1.206°，加了 Δψ（该历元 16.6″）会得到 −1.201°，不加则 −1.2056°。
/// 原因也说得通：Ω 本身就在平分点体系里，两边必须同口径。
fn libration_from_direction(lon_deg: f64, lat_deg: f64, t: Instant) -> Libration {
    /// 月球赤道对黄道的倾角（Cassini 定律），度。
    const I_DEG: f64 = 1.54242;

    // 升交点黄经 Ω = W3。用 ELP2000 自己的 W3，不再为这一个量另留一套引数。
    //
    // **要加黄经岁差**：W3 与 ELP 的 W1 一样是恒星口径（J2000 春分点起算），
    // 而传进来的方向已经换算到当日系了，两者不同口径会让 l′ 整差一个 p(T)
    // （1992 年 ~0.108°）。漏掉这一项时 Meeus 例 53.a 会差 0.108°，一眼可见。
    let p = argument_polynomials((t.tt_jd() - J2000) / 36525.0, 5);
    let om = p[2] * ARCSEC + precession_in_longitude(t); // 弧度
    let f = (p[0] - p[2]) * ARCSEC; // F = W1 − W3，弧度
    let (lam, bet) = (lon_deg * D2R, lat_deg * D2R);
    let inc = I_DEG * D2R;
    let w = lam - om; // 月球黄经与升交点黄经之差

    // Meeus 53.1–53.4
    let big_a =
        (w.sin() * bet.cos() * inc.cos() - bet.sin() * inc.sin()).atan2(w.cos() * bet.cos());
    let l = big_a - f;
    let b = (-w.sin() * bet.cos() * inc.sin() - bet.sin() * inc.cos())
        .clamp(-1.0, 1.0)
        .asin();

    Libration {
        lon_deg: norm180(l * R2D),
        lat_deg: b * R2D,
    }
}

/// 月球的光学天平动，**站心**（含周日天平动）。
///
/// 天平动的定义是"从观测者看到的月面中心在哪"，所以**必须给观测者**：观测者不在地心，
/// 看月面的角度与从地心看最多差 ~1°（周日天平动，随观测者经度改变）。地心那个极限
/// 情况由（私有的）`libration_from_direction` 直接表达，测试用它做 Meeus 例 53.a
/// 的回归。
///
/// 站心方向由（crate 内的）`Equatorial::topocentric_ecliptic` 给出——它做的是
/// 站心视差，与地平坐标同一个来源。
pub fn libration(t: Instant, observer: &crate::observer::Observer) -> Libration {
    let eq = moon_geocentric(t).apparent(t).equatorial(t);
    let (lon, lat) = eq.topocentric_ecliptic(t, observer);
    libration_from_direction(lon, lat, t)
}

/// 离 `t` 最近的一次朔的序号（允许 ±1 的误差）。
///
/// 与 [`new_moon`] 互逆，用平朔望月估算；`t` 只当作一个粗略的儒略日使用，
/// 因此可以安全地拿民用日编号近似（见 `nongli::new_moon_on_or_before`）。
pub fn new_moon_index_near(t: Instant) -> i64 {
    ((t.tt_jd() - EPOCH_NEW_MOON_JDE) / SYNODIC_MONTH_DAYS).round() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 对 **JPL DE421** 的真值（1900–2050 都在它的覆盖内）。
    ///
    /// 真值由 `tools/truthcheck.py` 那套算法算出后写死在这里——测试不该依赖
    /// `kernels/` 里那个不进版本库的文件。表里是**视**位置：月球的光行时按
    /// Δ/c 前推、再加周年光行差，与 `moon_geocentric` 的口径逐条对应。
    /// 容差 0.2″ / 0.05 km：实测三个历元的偏差是 0.02″–0.05″、0.009–0.028 km。
    ///
    /// 原先这里是 Meeus 例 47.a。换成 JPL 是因为**那个例子的值来自 Meeus 那 60 项
    /// 截断级数**，与真值本来就差 1.2″/4.2 km——拿它当判据反而比我们的实现还松。
    #[test]
    fn matches_jpl_de421() {
        for (jd, lam, bet, dist) in [
            (2448724.5, 133.156433018, -3.229189660, 368405.5430),
            (2451545.0, 223.312951927, 5.170871895, 402448.6401),
            (2469807.5, 18.654779809, 3.391949445, 378667.6805),
        ] {
            let e = moon_geocentric(Instant::from_tt(jd));
            let arcsec = |a: f64, b: f64| norm180(a - b).abs() * 3600.0;
            assert!(arcsec(e.lon_deg, lam) < 0.2, "jd={jd} dλ={}", arcsec(e.lon_deg, lam));
            assert!(arcsec(e.lat_deg, bet) < 0.2, "jd={jd} dβ={}", arcsec(e.lat_deg, bet));
            assert!(
                (e.distance_km - dist).abs() < 0.05,
                "jd={jd} dΔ={}",
                e.distance_km - dist
            );
        }
    }

    #[test]
    fn new_moon_is_syzygy() {
        use crate::sun::sun_apparent_longitude;
        for k in [-2000i64, -100, 0, 1, 337, 12000] {
            let t = new_moon(k);
            let d = norm180(moon_apparent_longitude(t) - sun_apparent_longitude(t));
            assert!(d.abs() < 1e-6, "k={k} d={d}");
        }
    }

    /// Meeus 例 53.a：1992-04-12.0 TD。天平动两个分量。
    ///
    /// 这里用的是**平**黄经（不加章动）——见 `libration_from_direction` 的说明。
    ///
    /// 容差 0.008° 而不是 0.002°：这里吃的是 `moon_geocentric` 的**视**方向
    /// （含 ~20″ 周年光行差），而 Meeus 53.a 给的是**几何**方向。20″ 折合 0.0056°，
    /// 在月面上约 37 km，显示上无所谓，但正好超出原来那个容差。
    ///
    /// 直接测 `libration_from_direction` 的**地心**输入：那是 Meeus 那个例子的口径，
    /// 也让这个内部函数有了消费者（否则它只是为"周日天平动"留着的死代码）。
    #[test]
    fn meeus_example_53a() {
        let t = Instant::from_tt(2448724.5);
        let e = moon_geocentric(t);
        let l = libration_from_direction(e.lon_deg, e.lat_deg, t);
        assert!((l.lon_deg - -1.206).abs() < 0.008, "l′={}", l.lon_deg);
        assert!((l.lat_deg - 4.194).abs() < 0.008, "b′={}", l.lat_deg);
    }

    /// 地心天平动（极限情况）。
    ///
    /// 教科书的 ±7.9°/±6.7° 说的就是这个；而 `libration` 给的是**站心**的，
    /// 含着周日天平动，幅度更大。
    fn geocentric_libration(t: Instant) -> Libration {
        let e = moon_geocentric(t);
        libration_from_direction(e.lon_deg, e.lat_deg, t)
    }

    /// 周日天平动：站心与地心之差。
    ///
    /// 这是"观测者不在地心"带来的那一部分，量级 ~1°（地球半径 6371 km 对月地距离
    /// 384400 km 张的角 ≈ 0.95°）。它必须**随观测者改变符号**——东边的观测者看到
    /// 月面偏一侧、西边的偏另一侧。地心天平动对所有人都一样，所以这一条正好把
    /// 两者分开。
    #[test]
    fn diurnal_libration_is_about_a_degree_and_flips_with_the_observer() {
        let t = Instant::from_tt(2451545.0 + 40.0);
        let geo = geocentric_libration(t);
        let east = crate::observer::Observer { lat_deg: 0.0, lon_deg: 90.0, height_m: 0.0 };
        let west = crate::observer::Observer { lat_deg: 0.0, lon_deg: -90.0, height_m: 0.0 };
        let de = libration(t, &east).lon_deg - geo.lon_deg;
        let dw = libration(t, &west).lon_deg - geo.lon_deg;
        assert!(
            (0.2..1.3).contains(&de.abs()) && (0.2..1.3).contains(&dw.abs()),
            "周日天平动量级不对：东 {de} 西 {dw}"
        );
        assert!(de * dw < 0.0, "两侧观测者的偏移应当反号：东 {de} 西 {dw}");
    }

    /// 天平动的**自洽**检查，不依赖外部参考值：
    ///
    /// * 光经天平动在月球过近地点/远地点附近应接近 0，在两者之间达到极值；
    /// * 光纬天平动在月球过升/降交点附近应接近 0，在两交点之间达到极值；
    /// * 两者幅度都应落在教科书值附近（经度 ±7.9°、纬度 ±6.7°）。
    #[test]
    fn libration_amplitudes_and_phases() {
        // 采样一个交点月与一个近点月
        let mut lo_lon: f64 = 0.0;
        let mut lo_lat: f64 = 0.0;
        let mut n = 0;
        for i in 0..(60 * 24) {
            let t = Instant::from_tt(2451545.0 + i as f64 / 24.0);
            let l = geocentric_libration(t);
            assert!(l.lon_deg.abs() < 9.0, "l={}", l.lon_deg);
            assert!(l.lat_deg.abs() < 8.0, "b={}", l.lat_deg);
            lo_lon = lo_lon.max(l.lon_deg.abs());
            lo_lat = lo_lat.max(l.lat_deg.abs());
            n += 1;
        }
        assert!(n > 1000);
        // 60 天足以覆盖完整的极值范围
        assert!((7.0..8.5).contains(&lo_lon), "经度天平动幅度 {lo_lon}");
        assert!((5.8..7.2).contains(&lo_lat), "纬度天平动幅度 {lo_lat}");
    }

    /// 天平动的两个分量应当以不同周期变化（一个是近点月、一个是交点月），
    /// 所以 60 天里不应出现"两者同时恒为 0"——那是公式退化的征兆。
    #[test]
    fn libration_is_not_degenerate() {
        let mut both_small = 0;
        for i in 0..(60 * 24) {
            let t = Instant::from_tt(2451545.0 + i as f64 / 24.0);
            let l = geocentric_libration(t);
            if l.lon_deg.abs() < 0.05 && l.lat_deg.abs() < 0.05 {
                both_small += 1;
            }
        }
        assert!(both_small < 24, "两者同时近零的小时数过多：{both_small}");
    }

    #[test]
    fn new_moon_index_roundtrip() {
        for k in [-2000i64, -1, 0, 1, 337, 12000] {
            assert_eq!(new_moon_index_near(new_moon(k)), k);
        }
    }

    /// 距离的物理量级：近地点 ~356 400 km，远地点 ~406 700 km。
    #[test]
    fn distance_stays_in_physical_range() {
        for k in 0..30 {
            let e = moon_geocentric(Instant::from_tt(2451545.0 + k as f64 * 1.0));
            assert!(
                (356_000.0..407_000.0).contains(&e.distance_km),
                "Δ={}",
                e.distance_km
            );
            assert!(e.lat_deg.abs() < 5.5, "β={}", e.lat_deg);
        }
    }
}
