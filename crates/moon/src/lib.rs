//! # moon —— 月相与月球位置
//!
//! 在 [`ephemeris`] 之上算月球"看起来什么样"：相位角、照亮比例、月龄、
//! 视直径、亮边方位角，以及出没与中天时刻。
//!
//! ```
//! use ephemeris::{Atmosphere, Instant, Observer};
//! use moon::{observe, rise_set};
//!
//! let t = Instant::from_tt(2451545.0);
//! let o = Observer { lat_deg: 39.9, lon_deg: 116.4, height_m: 44.0 };
//! let m = observe(t, &o, Some(Atmosphere::default()));
//! assert!((0.0..=1.0).contains(&m.phase.illuminated));
//! let r = rise_set(t, &o);
//! assert!(r.rise.is_some() || r.set.is_some());
//! ```
//!
//! ## 口径
//!
//! * 所有角度为**度**，距离 **km**，时刻 [`Instant`]。
//! * 照亮比例是**地心**所见；从地面看因视差有极小差别，不计。
//! * 出没的定义见 [`rise_set`]。

pub mod render;

use ephemeris::{
    libration, moon_geocentric, Atmosphere, Ecliptic, Equatorial, Horizontal, Instant, Libration,
    Observer, R2D,
};

/// 月球半径，km。用于视直径。
const MOON_RADIUS_KM: f64 = 1737.4;
/// 天文单位，km。
const AU_KM: f64 = 149_597_870.7;
/// 二分求根的容差，日。
const ROOT_TOL_DAYS: f64 = 1e-9;

/// 月球当前的相位与形状。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MoonPhase {
    /// 相位角 i，度：**日–月–地**夹角。0° = 望（全亮），180° = 朔（全暗）。
    pub phase_angle_deg: f64,
    /// 照亮比例 k = (1 + cos i)/2，0..=1。
    pub illuminated: f64,
    /// 月龄：距上一次**朔**的天数，0 ~ 29.5。
    pub age_days: f64,
    /// 是否在盈（朔→望）阶段。由黄经差判定。
    pub waxing: bool,
    /// 亮边方位角 χ，度，自北起向东为正：**月面亮的一侧朝向**。
    ///
    /// 画月牙朝向时用它。朔、望附近亮边退化，此时无意义。
    pub bright_limb_deg: f64,
    /// 地心距离，km。
    pub distance_km: f64,
    /// 视直径，度。
    pub angular_diameter_deg: f64,
    /// 月面中心的月面经纬度（天平动）。
    ///
    /// 与 [`render()`](crate::render::render) 无关——只画晨昏线用不到它；它决定的是**月面纹理**
    /// 该摆在哪。见 `docs/plan.md` 第 6 步末尾的说明。
    pub libration: Libration,
}

/// 月球的一次完整观测：位置 + 相位。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MoonObservation {
    /// 地心视赤道坐标。
    pub equatorial: Equatorial,
    /// 站心地平坐标（`atmosphere` 为 `Some` 时是含折射的**视**高度）。
    pub horizontal: Horizontal,
    /// 相位与形状。
    pub phase: MoonPhase,
}

/// 出没与中天时刻。`None` 表示所查的 24 小时内不发生——高纬地区月球可能整天
/// 不升或不落，这是真实情况而非错误。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct RiseSet {
    /// 升起：中心高度自下而上过 0。
    pub rise: Option<Instant>,
    /// 上中天：时角过 0。
    pub transit: Option<Instant>,
    /// 落下：中心高度自上而下过 0。
    pub set: Option<Instant>,
}

/// 月相的几何部分：给定月球**平**黄道坐标、太阳**几何**黄经与日地距离，算相位。
///
/// **三者必须同口径（都用当日的平分点、都不含光行差）**：相位角是日–月–地的
/// 纯几何三角形。混用"真分点"与"平分点"会让黄经差带上 Δψ（~14″），
/// 混入光行差则再差 ~20″——`tools/moonphasecheck.py` 里两条都实测到了。
///
/// 单独拆出来是为了能直接验证四个几何点（朔/上弦/望/下弦），不必构造真实历元。
pub fn phase(
    moon: Ecliptic,
    sun_lon_deg: f64,
    sun_distance_au: f64,
    t: Instant,
) -> MoonPhase {
    let dlon = (moon.lon_deg - sun_lon_deg).rem_euclid(360.0);
    let d = moon.distance_km;
    let r = sun_distance_au * AU_KM;

    // 相位角 (Meeus 48.3)：tan i = R sin ψ / (Δ − R cos ψ)。
    // 用 atan2 而非 atan，i 接近 90° 时才不失稳；再取绝对值——atan2 在 ψ>180°
    // 时给出 −i（亏相），而相位角按定义落在 [0,180]。
    // 照亮比例 (1+cos i)/2 对符号不敏感，但报出来的角度必须是正的。
    // 黄纬不能丢：Meeus 48.2 的"伸长"是 cos ψ = cos β · cos(λm − λs)。
    // 只取经度差会让相位角差到 0.1°（β 最大 5°，cos β 的影响约 β²/2）。
    let dlam = dlon.to_radians();
    let psi = (moon.lat_deg.to_radians().cos() * dlam.cos()).clamp(-1.0, 1.0).acos();
    // 亏相时 sin ψ 为负（下同），仍需靠 dlam 判断盈/亏
    let psi_signed = if dlam > core::f64::consts::PI {
        -psi
    } else {
        psi
    };
    let i = (r * psi_signed.sin()).atan2(d - r * psi_signed.cos()).abs();

    // 亮边方位角 (Meeus 48.5)：需要太阳与月球的**视赤道**坐标。
    // 太阳的黄纬按 0 处理（<1.2″）。
    // 亮边方位角要赤道坐标；`Ecliptic::equatorial` 内部会加章动
    let sun_eq = Ecliptic {
        lon_deg: sun_lon_deg,
        lat_deg: 0.0,
        distance_km: r,
    }
    .equatorial(t);
    let me = moon.equatorial(t);
    let dra = (sun_eq.ra_deg - me.ra_deg).to_radians();
    let (ssd, _) = sun_eq.dec_deg.to_radians().sin_cos();
    let (smd, cmd) = me.dec_deg.to_radians().sin_cos();
    let chi = (sun_eq.dec_deg.to_radians().cos() * dra.sin())
        .atan2(ssd * cmd - sun_eq.dec_deg.to_radians().cos() * smd * dra.cos());

    MoonPhase {
        phase_angle_deg: i * R2D,
        illuminated: (1.0 + i.cos()) / 2.0,
        age_days: age_days(t),
        waxing: dlon < 180.0,
        bright_limb_deg: chi.rem_euclid(core::f64::consts::TAU) * R2D,
        distance_km: d,
        angular_diameter_deg: 2.0 * (MOON_RADIUS_KM / d).asin() * R2D,
        libration: libration(t),
    }
}

/// 距上一次**朔**的天数。
pub fn age_days(t: Instant) -> f64 {
    let mut k = ephemeris::new_moon_index_near(t);
    let mut last = ephemeris::new_moon(k);
    while last.tt_jd() > t.tt_jd() {
        k -= 1;
        last = ephemeris::new_moon(k);
    }
    while ephemeris::new_moon(k + 1).tt_jd() <= t.tt_jd() {
        k += 1;
        last = ephemeris::new_moon(k);
    }
    t.tt_jd() - last.tt_jd()
}

/// 月球的地心视赤道坐标。
pub fn equatorial(t: Instant) -> Equatorial {
    moon_geocentric(t).apparent(t).equatorial(t)
}

/// 观测月球：位置、地平坐标、相位。
pub fn observe(t: Instant, observer: &Observer, atmosphere: Option<Atmosphere>) -> MoonObservation {
    // 位置只算一次，赤道坐标与相位都用它。
    let geo = moon_geocentric(t);
    let eq = geo.apparent(t).equatorial(t);
    MoonObservation {
        equatorial: eq,
        horizontal: eq.horizontal(t, observer, atmosphere),
        // 相位用**几何**量：月球取平黄道（`geo`），太阳取几何黄经。
        // 两者同为当日的平分点、都不含光行差，三角形才是对的。
        phase: phase(
            geo,
            ephemeris::sun_geometric_longitude(t),
            ephemeris::sun_distance_au(t),
            t,
        ),
    }
}

/// 由一个视赤道坐标出发算 `(几何地平高度, 本地时角)`，弧度。
///
/// 两者都要用同一时刻的月球位置，所以放在一起算——早先在扫描里分别调用，
/// 每个采样点把月球位置算了两遍。
fn altitude_and_hour_angle(eq: Equatorial, t: Instant, observer: &Observer) -> (f64, f64) {
    let alt = eq.horizontal(t, observer, None).altitude_deg.to_radians();
    let last = (ephemeris::gast(t) + observer.lon_deg).to_radians();
    let ha = (last - eq.ra_deg.to_radians() + core::f64::consts::PI)
        .rem_euclid(core::f64::consts::TAU)
        - core::f64::consts::PI;
    (alt, ha)
}

/// 找出 `t0` 起 24 小时内的出没与上中天。
///
/// **定义**：升/落是月球**中心**的几何高度自下/自上穿过 0——既不含折射，也不含
/// 视差带来的半直径修正（真实的出没比这个早/晚几分钟，且随天气变）。中天是
/// 时角穿过 0。
///
/// 高纬地区月球可能一天不升或一天不落，对应字段就是 `None`。
pub fn rise_set(t0: Instant, observer: &Observer) -> RiseSet {
    const STEPS: i64 = 48; // 半小时一格
    let td = |i: i64| Instant::from_tt(t0.tt_jd() + i as f64 / STEPS as f64);

    // 扫描用同一个 Equatorial 同时得到高度与时角
    let sample = |i: i64| {
        let t = td(i);
        let (a, h) = altitude_and_hour_angle(equatorial(t), t, observer);
        (t, a, h)
    };
    let altitude_at = |t: Instant| altitude_and_hour_angle(equatorial(t), t, observer).0;
    let hour_angle_at = |t: Instant| altitude_and_hour_angle(equatorial(t), t, observer).1;

    let mut rise = None;
    let mut set = None;
    let mut transit = None;
    let mut prev = sample(0);
    for i in 1..=STEPS {
        let cur = sample(i);
        let (pt, pa, ph) = prev;
        let (ct, ca, ch) = cur;
        if pa <= 0.0 && ca > 0.0 {
            rise = Some(bisect(pt, ct, altitude_at));
        }
        if pa >= 0.0 && ca < 0.0 {
            set = Some(bisect(pt, ct, altitude_at));
        }
        if ph < 0.0 && ch >= 0.0 {
            transit = Some(bisect(pt, ct, hour_angle_at));
        }
        prev = cur;
    }
    RiseSet {
        rise,
        transit,
        set,
    }
}

/// 二分求根：`f` 在 `a`、`b` 处异号。
fn bisect<F: Fn(Instant) -> f64>(a: Instant, b: Instant, f: F) -> Instant {
    let (mut lo, mut hi) = (a.tt_jd(), b.tt_jd());
    let flo = f(a);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if (f(Instant::from_tt(mid)) <= 0.0) == (flo <= 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < ROOT_TOL_DAYS {
            break;
        }
    }
    Instant::from_tt(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    const T0: f64 = 2451545.0;
    /// 视差可忽略的远距离，用来隔离几何。
    const FAR: f64 = 384_400.0;

    fn phase_at_elongation(elong: f64) -> MoonPhase {
        let t = Instant::from_tt(T0);
        phase(
            Ecliptic {
                lon_deg: elong,
                lat_deg: 0.0,
                distance_km: FAR,
            },
            0.0,
            1.0,
            t,
        )
    }

    /// 四个几何点。黄经差 ψ = 月 − 日。
    ///
    /// **注意 0.5 不等于 ψ=90°**：照亮 0.5 发生在**相位角** i = 90°，
    /// 而 i 与 ψ 相差一个"日月距离之比"的小量——ψ=90° 时 i ≈ 89.85°，
    /// 照亮 ≈ 0.5013。这是真实的几何，不是误差，所以这里分开测。
    #[test]
    fn cardinal_points() {
        // 朔：ψ = 0 → 全暗
        assert!(phase_at_elongation(0.0).illuminated.abs() < 1e-12);
        // 望：ψ = 180 → 全亮
        assert!((phase_at_elongation(180.0).illuminated - 1.0).abs() < 1e-12);
        // 上/下弦附近
        let p90 = phase_at_elongation(90.0);
        assert!((p90.illuminated - 0.5).abs() < 2e-3, "ψ=90 照亮 {:.6}", p90.illuminated);
        assert!((p90.phase_angle_deg - 89.85).abs() < 0.05, "i={}", p90.phase_angle_deg);

        // 让相位角恰为 90° 的那个 ψ：cos ψ = Δ/R ⇒ 照亮恰为 0.5
        let d = FAR;
        let r = AU_KM;
        let psi = (d / r).acos().to_degrees();
        let p = phase_at_elongation(psi);
        assert!(
            (p.illuminated - 0.5).abs() < 1e-9,
            "i=90 处照亮应为 0.5，得 {:.12}",
            p.illuminated
        );
        assert!((p.phase_angle_deg - 90.0).abs() < 1e-9);
    }

    /// 相位角的**独立**推法：用日月地的三角形按余弦定理算，与级数公式比。
    ///
    /// 这是真正独立的一遍——不是把 `phase_from_apparent` 的公式再抄一次。
    #[test]
    fn phase_angle_matches_law_of_cosines() {
        for i in 0..120 {
            let t = Instant::from_tt(T0 + i as f64 * 3.0);
            let ap = moon_geocentric(t);
            let sun_lon = ephemeris::sun_geometric_longitude(t);
            let r = ephemeris::sun_distance_au(t) * AU_KM;
            let p = phase(ap, sun_lon, ephemeris::sun_distance_au(t), t);

            // 地心直角坐标（几何量，与上面同口径）
            let (lam_m, bet_m) = (ap.lon_deg.to_radians(), ap.lat_deg.to_radians());
            let m = [
                ap.distance_km * bet_m.cos() * lam_m.cos(),
                ap.distance_km * bet_m.cos() * lam_m.sin(),
                ap.distance_km * bet_m.sin(),
            ];
            let lam_s = sun_lon.to_radians();
            let s = [r * lam_s.cos(), r * lam_s.sin(), 0.0];
            // 日–月距离
            let dsm = ((m[0] - s[0]).powi(2) + (m[1] - s[1]).powi(2) + (m[2] - s[2]).powi(2)).sqrt();
            // 余弦定理，顶点取在**月**上：它是"月→日"与"月→地"的夹角。
            // 相邻边 Δ(=月地) 与 dsm(=月日)，对边 r(=日地)。
            let cos_i = (ap.distance_km * ap.distance_km + dsm * dsm - r * r)
                / (2.0 * ap.distance_km * dsm);
            let i_law = cos_i.clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                (i_law - p.phase_angle_deg).abs() < 0.01,
                "t={}: 级数 {:.5}° vs 余弦定理 {:.5}°",
                t.tt_jd(),
                p.phase_angle_deg,
                i_law
            );
        }
    }

    /// 盈（朔→望）与亏（望→朔）的判定。
    #[test]
    fn waxing_flag() {
        assert!(phase_at_elongation(45.0).waxing); // 朔后，盈
        assert!(phase_at_elongation(90.0).waxing); // 上弦
        assert!(!phase_at_elongation(180.0).waxing); // 望
        assert!(!phase_at_elongation(270.0).waxing); // 下弦，亏
        assert!(!phase_at_elongation(315.0).waxing); // 朔前，亏
    }

    /// 月龄与 `new_moon` 一致：在朔处为 0，在望附近约半个朔望月。
    #[test]
    fn age_matches_new_moon() {
        for k in [0i64, 100, -50] {
            let nm = ephemeris::new_moon(k);
            assert!(age_days(nm).abs() < 1e-6, "k={k} age={}", age_days(nm));
            // 朔后 14 天，月龄应约为 14
            let t = Instant::from_tt(nm.tt_jd() + 14.0);
            let a = age_days(t);
            assert!((a - 14.0).abs() < 1e-6, "age={a}");
        }
        // 月龄总落在 [0, 朔望月) 内。朔望月本身在 29.27–29.83 天之间变化，
        // 所以上界取 29.9 而不是平均值。
        for i in 0..60 {
            let a = age_days(Instant::from_tt(T0 + i as f64 * 0.5));
            assert!((0.0..29.9).contains(&a), "age={a}");
        }
    }

    /// 视直径应在 29.3′ ~ 34.1′ 之间（近地点/远地点）。
    #[test]
    fn angular_diameter_range() {
        for i in 0..200 {
            let t = Instant::from_tt(T0 + i as f64);
            let p = phase(
                moon_geocentric(t),
                ephemeris::sun_geometric_longitude(t),
                ephemeris::sun_distance_au(t),
                t,
            );
            let arcmin = p.angular_diameter_deg * 60.0;
            assert!((29.0..34.5).contains(&arcmin), "视直径 {arcmin}′");
        }
    }

    /// 出没：两次独立方法互相印证。
    ///
    /// 方法一（这里）：对高度**二分求根**。
    /// 方法二：用"月球大约每小时移动 0.55°"做线性割线迭代，从粗略时刻收敛。
    #[test]
    fn rise_set_two_methods_agree() {
        let o = Observer {
            lat_deg: 39.9042,
            lon_deg: 116.4074,
            height_m: 44.0,
        };
        let mut checked = 0;
        for d in 0..40 {
            let t0 = Instant::from_tt(T0 + d as f64);
            let rs = rise_set(t0, &o);
            for (target, label) in [(rs.rise, "rise"), (rs.set, "set")] {
                let Some(t_root) = target else { continue };
                // 方法二：以根部邻近时刻出发做割线迭代
                let sol = secant_altitude(t_root, &o);
                let diff = (sol.tt_jd() - t_root.tt_jd()).abs() * 86400.0;
                assert!(diff < 1.0, "{label}: 两法相差 {diff} 秒");
                // 根部高度确实过零
                let h = altitude_and_hour_angle(equatorial(t_root), t_root, &o).0.to_degrees() * 3600.0;
                assert!(h.abs() < 0.05, "{label}: 根部高度 {h}\"");
                checked += 1;
            }
        }
        assert!(checked > 40, "样本太少 ({checked})");
    }

    /// 割线法：从 `guess` 附近用自己的高度公式独立收敛一遍。
    fn secant_altitude(guess: Instant, o: &Observer) -> Instant {
        let (mut a, mut b) = (guess.tt_jd() - 0.02, guess.tt_jd() + 0.02);
        for _ in 0..60 {
            let (fa, fb) = (
                altitude_and_hour_angle(equatorial(Instant::from_tt(a)), Instant::from_tt(a), o).0,
                altitude_and_hour_angle(equatorial(Instant::from_tt(b)), Instant::from_tt(b), o).0,
            );
            if (fb - fa).abs() < 1e-15 {
                break;
            }
            let c = b - fb * (b - a) / (fb - fa);
            a = b;
            b = c;
            if (b - a).abs() < 1e-10 {
                break;
            }
        }
        Instant::from_tt(b)
    }

    /// 高纬：有些日子就是不升不落，应当返回 None 而不是硬凑一个。
    #[test]
    fn polar_site_may_have_no_rise() {
        let o = Observer {
            lat_deg: 89.0,
            lon_deg: 0.0,
            height_m: 0.0,
        };
        let mut none_days = 0;
        for d in 0..28 {
            let rs = rise_set(Instant::from_tt(T0 + d as f64), &o);
            if rs.rise.is_none() {
                none_days += 1;
            }
        }
        assert!(none_days > 0, "极地总该有几天不升");
    }
}
