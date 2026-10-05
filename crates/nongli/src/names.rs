//! GB/T 33661-2017 第 6 章的命名方法：干支、生肖、月名、日名。

/// 天干。
pub(crate) const GAN: [&str; 10] = ["甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"];
/// 地支。
pub(crate) const ZHI: [&str; 12] = [
    "子", "丑", "寅", "卯", "辰", "巳", "午", "未", "申", "酉", "戌", "亥",
];
/// 生肖 (与地支一一对应)。
pub(crate) const ZODIAC: [&str; 12] = [
    "鼠", "牛", "虎", "兔", "龙", "蛇", "马", "羊", "猴", "鸡", "狗", "猪",
];
/// 二十四节气，自冬至起 (附录 A)。
pub const TERM_NAMES: [&str; 24] = [
    "冬至", "小寒", "大寒", "立春", "雨水", "惊蛰", "春分", "清明", "谷雨", "立夏",
    "小满", "芒种", "夏至", "小暑", "大暑", "立秋", "处暑", "白露", "秋分", "寒露",
    "霜降", "立冬", "小雪", "大雪",
];
/// 数序纪月 (6.2)。
pub(crate) const MONTH_NAMES: [&str; 12] = [
    "正月", "二月", "三月", "四月", "五月", "六月", "七月", "八月", "九月", "十月",
    "十一月", "十二月",
];

const DIGITS: [&str; 10] = ["一", "二", "三", "四", "五", "六", "七", "八", "九", "十"];

/// 六十干支名，`i` 取模 60，0 = 甲子。
pub(crate) fn ganzhi(i: i64) -> String {
    let i = i.rem_euclid(60) as usize;
    format!("{}{}", GAN[i % 10], ZHI[i % 12])
}

/// 数序纪日 (6.3.1)。第 21–29 日用"廿一…廿九"。
pub fn day_name(day: u32) -> String {
    match day {
        1..=10 => format!("初{}", DIGITS[(day - 1) as usize]),
        11..=19 => format!("十{}", DIGITS[(day - 11) as usize]),
        20 => "二十".to_string(),
        21..=29 => format!("廿{}", DIGITS[(day - 21) as usize]),
        30 => "三十".to_string(),
        _ => panic!("农历日超出 1..=30: {day}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ganzhi_cycle() {
        assert_eq!(ganzhi(0), "甲子");
        assert_eq!(ganzhi(1), "乙丑");
        assert_eq!(ganzhi(41), "乙巳");
        assert_eq!(ganzhi(59), "癸亥");
        assert_eq!(ganzhi(60), "甲子");
        assert_eq!(ganzhi(-1), "癸亥");
    }

    #[test]
    fn days() {
        assert_eq!(day_name(1), "初一");
        assert_eq!(day_name(10), "初十");
        assert_eq!(day_name(11), "十一");
        assert_eq!(day_name(20), "二十");
        assert_eq!(day_name(21), "廿一");
        assert_eq!(day_name(29), "廿九");
        assert_eq!(day_name(30), "三十");
    }
}
