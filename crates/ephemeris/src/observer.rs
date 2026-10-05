//! 观测者、恒星时与地平坐标。
//!
//! 把"地心视赤道坐标"化为"从地面上某点看到的地平坐标"：
//!
//! 1. **恒星时** —— [`gmst`] / [`gast`]，把地球自转角接到天球；
//! 2. **站心视差** —— 观测者不在地心，矢量要减去观测者的地心位置；
//! 3. **周日光行差** —— 观测者随地球自转的速度带来的 ~0.32″ 偏移；
//! 4. **大气折射** —— 可选，见 [`Atmosphere`]。
//!
//! 整条链都在**当日真赤道/真分点**系里做（赤经的 x 轴取真分点，故用 GAST 而非
//! GMST）。这避免了任何岁差变换——地平坐标本来就定义在"当日"，不需要 J2000。
//!
//! # 口径
//!
//! 输入 [`Equatorial`] 必须已是**视**位置（见 [`crate::frames::Ecliptic::equatorial`]），
//! 且距离是**地心**距离。输出的方位角自**北**起向东为正。

use crate::angle::{norm360, D2R, R2D};
use crate::frames::{mean_obliquity, nutation, Equatorial};
use crate::time::{Instant, J2000};

/// WGS84 椭球：长半轴 km 与扁率。
const WGS84_A_KM: f64 = 6378.137;
const WGS84_F: f64 = 1.0 / 298.257223563;
/// 地球自转角速度，rad/s（ERFA 的 `ERFA_OM`）。
const EARTH_OMEGA: f64 = 7.292115e-5;
/// 光速，km/s。
const C_KM_S: f64 = 299792.458;

/// 地面观测者。纬度东经为正，海拔以米计（可为负）。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Observer {
    /// 大地纬度，度，北正。
    pub lat_deg: f64,
    /// 大地经度，度，东正。
    pub lon_deg: f64,
    /// 椭球高，米。
    pub height_m: f64,
}

/// 地平坐标。
///
/// **方位角在天顶没有定义**（几何退化，方位角公式的分母趋零），极点也一样。
/// 那附近返回的方位角是数值噪声，调用方应忽略；高度角不受影响。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Horizontal {
    /// 方位角，度，自**北**起向东为正，[0,360)。天顶附近无意义。
    pub azimuth_deg: f64,
    /// 地平高度，度。折射后（若启用）为**视**高度。
    pub altitude_deg: f64,
}

/// 大气条件，用于折射。取 [`Atmosphere::default`] 即标准大气。
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Atmosphere {
    /// 地面气压，hPa。
    pub pressure_hpa: f64,
    /// 地面气温，摄氏度。
    pub temperature_c: f64,
}

impl Default for Atmosphere {
    /// 标准大气：1010 hPa、10 °C（`erfa.refco` 的参考条件附近）。
    fn default() -> Self {
        Self {
            pressure_hpa: 1010.0,
            temperature_c: 10.0,
        }
    }
}

/// Greenwich **平**恒星时，度，[0,360)。
///
/// IAU 2006：`GMST = ERA(UT1) + 多项式(TT)`。ERA 是地球自转角。
/// ERA 的系数含 1.00273781191135448，即恒星日与平太阳日之比。
pub fn gmst(t: Instant) -> f64 {
    let t_c = (t.tt_jd() - J2000) / 36525.0;
    let era = core::f64::consts::TAU
        * (0.7790572732640 + 1.0027378119113545 * (t.ut1_jd() - J2000));
    // 角秒 → 弧度；这是 IAU 2006 的 GMST 与 ERA 之差
    let poly = 0.014506
        + 4612.156534 * t_c
        + 1.3915817 * t_c * t_c
        - 0.00000044 * t_c.powi(3)
        - 0.000029956 * t_c.powi(4)
        - 0.0000000368 * t_c.powi(5);
    norm360((era + poly / 3600.0 * D2R) * R2D)
}

/// Greenwich **真**恒星时，度，[0,360)。
///
/// `GAST = GMST + 分点差`。分点差按 IAU 2000 的定义（ERFA `ee00.c`：
/// `ee = Δψ·cos(ε_A) + eect00`）取**平**黄赤交角 ε_A，不是真交角。
///
/// 曾经写成 `true_obliquity`（ε_A + Δε）——错得小（3e-4″），但定义就是平交角。
/// 交角章动 Δε 进的是赤道坐标（见 [`Equatorial`] 那一侧），不进恒星时。
///
/// 补充项 `eect00`（±0.003″）本实现未含。
pub fn gast(t: Instant) -> f64 {
    let ee = nutation(t).dpsi_deg * mean_obliquity(t).to_radians().cos();
    norm360(gmst(t) + ee)
}

impl Observer {
    /// 观测者的地心位置矢量，**地球固连**系（x 轴指向经度 0），km。
    ///
    /// WGS84 椭球上的标准公式。这一步与时间无关，所以单独拆出来——
    /// 可以直接对 `erfa.gd2gc` 验证，不必先绕开恒星时。
    pub fn geocentric_fixed_km(&self) -> [f64; 3] {
        let phi = self.lat_deg * D2R;
        let lam = self.lon_deg * D2R;
        let h_km = self.height_m / 1000.0;
        let e2 = WGS84_F * (2.0 - WGS84_F);
        let (sphi, cphi) = phi.sin_cos();
        let n = WGS84_A_KM / (1.0 - e2 * sphi * sphi).sqrt();
        [
            (n + h_km) * cphi * lam.cos(),
            (n + h_km) * cphi * lam.sin(),
            (n * (1.0 - e2) + h_km) * sphi,
        ]
    }

    /// 观测者的地心位置矢量，**当日真赤道**系（x 轴取真分点），km。
    ///
    /// 参数是 **GAST**（度）而不是 [`Instant`]：GAST 要算一次 IAU 2000A 章动
    /// （~11 µs），而一次地平坐标换算里要用到它三次（观测者位置、观测者速度、时角）。
    /// 让调用方一次算好传进来，就不会重复。见 [`Equatorial::horizontal`]。
    fn geocentric_km(&self, gast_deg: f64) -> [f64; 3] {
        let [x, y, z] = self.geocentric_fixed_km();
        let th = gast_deg * D2R;
        let (sth, cth) = th.sin_cos();
        [x * cth - y * sth, x * sth + y * cth, z]
    }

    /// 观测者随地球自转的速度，km/s，同日坐标系。用于周日光行差。参数同 [`Observer::geocentric_km`]。
    fn velocity_km_s(&self, gast_deg: f64) -> [f64; 3] {
        let r = self.geocentric_km(gast_deg);
        // ω 沿 z 轴：v = ω × r
        [-EARTH_OMEGA * r[1], EARTH_OMEGA * r[0], 0.0]
    }
}

impl Equatorial {
    /// 视赤道坐标 → 站心地平坐标。
    ///
    /// 依次做：站心视差 → 周日光行差 → 时角 → 地平坐标 →（可选）大气折射。
    ///
    /// `atmosphere` 为 `Some` 时计入折射，输出的是**视**高度（即肉眼看到的）；
    /// 为 `None` 时是真空中的几何高度。这两者相差在地平附近可达 34′，
    /// 所以这个选择不是可有可无的装饰。
    pub fn horizontal(
        self,
        t: Instant,
        observer: &Observer,
        atmosphere: Option<Atmosphere>,
    ) -> Horizontal {
        self.horizontal_at(gast(t), observer, atmosphere)
    }

    /// 同上，但 GAST 由调用方给定。**不需要** `Instant`——GAST 之后的步骤都与时间无关。
    ///
    /// 求 GAST 要算一次 IAU 2000A 章动（~11 µs）。调用方若在同一个瞬间还要拿 GAST
    /// 做别的事（例如求时角），自己算一次传进来就省下一遍。`Observer` 内部那两个
    /// 辅助方法（地心矢量、观测者速度）也是同一个理由改收 GAST。
    pub fn horizontal_at(
        self,
        gast_deg: f64,
        observer: &Observer,
        atmosphere: Option<Atmosphere>,
    ) -> Horizontal {
        // 地心矢量（km）
        let (sa, ca) = (self.ra_deg * D2R).sin_cos();
        let (sd, cd) = (self.dec_deg * D2R).sin_cos();
        let d = self.distance_km;
        let geo = [d * cd * ca, d * cd * sa, d * sd];

        // 站心视差：减去观测者的地心位置
        let obs = observer.geocentric_km(gast_deg);
        let mut v = [geo[0] - obs[0], geo[1] - obs[1], geo[2] - obs[2]];

        // 周日光行差：观测者速度 ~0.46 km/s → ~0.32″
        let vel = observer.velocity_km_s(gast_deg);
        let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        for k in 0..3 {
            v[k] += vel[k] / C_KM_S * r;
        }

        // 站心赤道坐标
        let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let dec = (v[2] / r).clamp(-1.0, 1.0).asin();
        let ra = v[1].atan2(v[0]);

        // 地方真恒星时 → 时角
        let last = (gast_deg + observer.lon_deg) * D2R;
        let h = last - ra;

        // 地平坐标：方位角自北起向东为正
        let phi = observer.lat_deg * D2R;
        let (sphi, cphi) = phi.sin_cos();
        let (sdec, cdec) = dec.sin_cos();
        let (sh, ch) = h.sin_cos();
        let alt = (sphi * sdec + cphi * cdec * ch).clamp(-1.0, 1.0).asin();
        let az = norm360(
            (-cdec * sh).atan2(sdec * cphi - cdec * sphi * ch) * R2D,
        );

        let alt_deg = match atmosphere {
            None => alt * R2D,
            Some(a) => (alt + refraction(alt, &a)) * R2D,
        };
        Horizontal {
            azimuth_deg: az,
            altitude_deg: alt_deg,
        }
    }
}

/// 大气折射改正，弧度。输入为**真空**地平高度，返回需加上的量。
///
/// 模型取 SOFA/ERFA `refco` 的形式 `dZ = A·tan Z + B·tan³ Z`（Z 为视天顶距），
/// 但 A、B **不是**调用 ERFA，而是拟合出来的：
///
/// ```text
/// A = P·(c1/T + c2),  B = P·(d1/T + d2)     P 为 hPa，T 为开尔文
/// ```
///
/// 在 P∈[950,1050]、T∈[−20,40]°C 上复现 `erfa.refco` 的 A、B 到 **< 0.09″**
/// （见 `tools/observercheck.py`）。ERFA 的模型本身就与光线追迹相差 8 mas RMS
/// （SOFA 文档），所以这个拟合不是瓶颈。
///
/// **已知取舍**：把湿度固定为 0.5、波长固定为 0.574 µm（黄光）。湿度偏离
/// 0.5 会带来 ≤0.29″，波长偏到 0.4 µm 会带来 ≤1.3″（色散）。要做全色散需要
/// Edlén 的空气折射率公式，暂未做——月球理论本身是 2.9″，且折射还随天气变。
///
/// 由于 `refco` 的定义是"给**视**天顶距加 dZ 得到真空天顶距"，而我们是从真空
/// 高度出发求视高度，故需迭代：`Z ← Z − dZ(Z)`，三次即收敛到 1e-9 弧度。
fn refraction(altitude_rad: f64, atm: &Atmosphere) -> f64 {
    let z0 = core::f64::consts::FRAC_PI_2 - altitude_rad;
    if z0 >= core::f64::consts::FRAC_PI_2 {
        return 0.0; // 已在地平以下，折射模型不适用
    }
    let tk = atm.temperature_c + 273.15;
    let a = atm.pressure_hpa * (16.619485936 / tk - 1.388866696e-3) / 3600.0 * D2R;
    let b = atm.pressure_hpa * (4.223473655e-3 / tk - 7.909970381e-5) / 3600.0 * D2R;

    let mut z = z0;
    for _ in 0..3 {
        let tz = z.tan();
        z = z0 - (a * tz + b * tz * tz * tz);
    }
    z0 - z
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gmst_advances_by_sidereal_day() {
        // 一个平太阳日之后，GMST 前进 360.9856°，即多出约 0.9856°
        let a = gmst(Instant::from_tt(2451545.0));
        let b = gmst(Instant::from_tt(2451546.0));
        let d = norm360(b - a);
        assert!((d - 0.985647).abs() < 0.01, "ΔGMST={d}");
    }

    #[test]
    fn zenith_gives_altitude_90() {
        // 天顶方向是**大地**垂线方向，不是地心向径方向——两者最多差 11.5′
        // （40° 处约 0.19°），搞混了会得到一个看似"差一点"的错。
        let t = Instant::from_tt(2451545.0);
        let o = Observer {
            lat_deg: 40.0,
            lon_deg: 116.0,
            height_m: 50.0,
        };
        let eq = Equatorial {
            ra_deg: norm360(gast(t) + o.lon_deg),
            dec_deg: o.lat_deg,
            distance_km: 1.0e9, // 远到视差可忽略
        };
        let hz = eq.horizontal(t, &o, None);
        // 残差应只剩周日光行差 ~0.32″
        assert!(
            (hz.altitude_deg - 90.0).abs() * 3600.0 < 1.0,
            "alt={}",
            hz.altitude_deg
        );
        // 天顶处方位角**没有定义**（几何退化），故不检查 az。
    }

    #[test]
    fn equator_rising_star_is_due_east() {
        // 赤道上、赤纬 0 的天体在升起时刻应正东（时角 −90°）
        let t = Instant::from_tt(2451545.0);
        let o = Observer {
            lat_deg: 0.0,
            lon_deg: 0.0,
            height_m: 0.0,
        };
        // 让时角恰为 −90°：ra = GAST + lon + 90°
        let ra = norm360(gast(t) + 90.0);
        let eq = Equatorial {
            ra_deg: ra,
            dec_deg: 0.0,
            distance_km: 1.0e9, // 远到视差可忽略
        };
        let hz = eq.horizontal(t, &o, None);
        assert!((hz.azimuth_deg - 90.0).abs() < 0.01, "az={}", hz.azimuth_deg);
        assert!(hz.altitude_deg.abs() < 0.01, "alt={}", hz.altitude_deg);
    }

    #[test]
    fn refraction_grows_toward_horizon() {
        let atm = Atmosphere::default();
        let a90 = refraction(89.0_f64.to_radians(), &atm) * R2D * 3600.0;
        let a45 = refraction(45.0_f64.to_radians(), &atm) * R2D * 3600.0;
        let a10 = refraction(10.0_f64.to_radians(), &atm) * R2D * 3600.0;
        // 89° 处天顶距 1°，tan1°≈0.0175，A≈57.9″ ⇒ 约 1.0″
        assert!((a90 - 1.01).abs() < 0.05, "89° 折射 {a90}\"");
        assert!((a45 - 57.9).abs() < 1.0, "45° 折射 {a45}\"");
        assert!((a10 - 316.0).abs() < 3.0, "10° 折射 {a10}\"");
    }
}
