//! 命令行：农历日期换算与年历打印。
//!
//! ```text
//! nongli 2026-02-17              # 公历 -> 农历
//! nongli 2026-02-17T15:30:00     # 带时刻 (北京时间)
//! nongli --utc 2026-02-17T07:30  # 输入按 UTC
//! nongli --julian 1500-01-01     # 输入按儒略历
//! nongli --year 2026             # 打印全年公农历对照 + 节气
//! nongli --sui 2026              # 打印 2026 岁(含闰月判定)
//! ```

use ephemeris::{jdn_from_ymd, new_moon, ymd_from_jdn, Calendar, DateTime};
use nongli::calendar::{new_moon_on_or_before, winter_solstice_jdn};
use nongli::{beijing_jdn, LunarDate, TERM_NAMES};

const WEEKDAYS: [&str; 7] = ["日", "一", "二", "三", "四", "五", "六"];

fn weekday_name(jdn: i64) -> &'static str {
    WEEKDAYS[(jdn + 1).rem_euclid(7) as usize]
}

fn print_date(d: &LunarDate) {
    let (y, m, dd) = d.gregorian();
    println!(
        "公历 {y:04}-{m:02}-{dd:02} (星期{})",
        weekday_name(d.jdn)
    );
    println!(
        "{d}   [{}{}年 {}{}]",
        d.ganzhi_year(),
        d.zodiac(),
        d.month_name(),
        d.day_name()
    );
    println!(
        "干支  {}年 {}日",
        d.ganzhi_year(),
        d.ganzhi_day()
    );
    if let Some(t) = d.solar_term() {
        println!("节气  {t}");
    }
    // 下一节气
    let next = nongli::solar_terms_between(d.jdn, d.jdn + 20);
    if let Some(&(j, name)) = next.iter().find(|&&(j, _)| j > d.jdn) {
        let dt = j - d.jdn;
        println!("下一节气  {name} ({} 天后)", dt);
    }
}

fn parse_input(s: &str) -> Option<DateTime> {
    let (date, time) = match s.split_once(['T', ' ']) {
        Some((a, b)) => (a, Some(b)),
        None => (s, None),
    };
    let mut it = date.split('-');
    let y: i64 = if date.starts_with('-') {
        it.next(); // 开头的空串
        -it.next()?.parse::<i64>().ok()?
    } else {
        it.next()?.parse::<i64>().ok()?
    };
    let m: u32 = it.next()?.parse().ok()?;
    let d: u32 = it.next()?.parse().ok()?;
    let (h, mi, se) = match time {
        None => (0, 0, 0.0),
        Some(t) => {
            let p: Vec<&str> = t.split(':').collect();
            (
                p.first().and_then(|x| x.parse().ok()).unwrap_or(0),
                p.get(1).and_then(|x| x.parse().ok()).unwrap_or(0),
                p.get(2).and_then(|x| x.parse().ok()).unwrap_or(0.0),
            )
        }
    };
    Some(DateTime::new(y, m, d, h, mi, se))
}

fn print_year(year: i64) {
    println!("===== {year} 年公历/农历对照 =====");
    let mut jdn = jdn_from_ymd(year, 1, 1, Calendar::Gregorian);
    let end = jdn_from_ymd(year, 12, 31, Calendar::Gregorian);
    while jdn <= end {
        let (_, m, d) = ymd_from_jdn(jdn, Calendar::Gregorian);
        let l = LunarDate::from_jdn(jdn);
        let term = l.solar_term().map(|t| format!("  ★{t}")).unwrap_or_default();
        let year_part = if l.day == 1 {
            format!("农历{}年", l.ganzhi_year())
        } else {
            "农历".to_string()
        };
        println!(
            "{year:04}-{m:02}-{d:02}  {year_part}{}{}{}",
            l.month_name(),
            l.day_name(),
            term
        );
        jdn += 1;
    }
}

/// 三伏。规则见 `nongli::sanfu` 的文档。
fn print_sanfu(year: i64) {
    let s = nongli::sanfu(year);
    let d = |j: i64| {
        let (yy, mm, dd) = ymd_from_jdn(j, Calendar::Gregorian);
        format!("{yy:04}-{mm:02}-{dd:02}")
    };
    println!("===== {year} 年三伏 =====");
    println!("初伏  {} 起  10 天", d(s.first_jdn));
    println!(
        "中伏  {} 起  {} 天",
        d(s.middle_jdn),
        s.middle_days()
    );
    println!("末伏  {} 起  10 天", d(s.last_jdn));
    println!("出伏  {}", d(s.end_jdn()));
}

fn print_sui(year: i64) {
    println!("===== 以 {year} 年冬至所在月为十一月的\"岁\" =====");
    let month11 = new_moon_on_or_before(winter_solstice_jdn(year));
    let next11 = new_moon_on_or_before(winter_solstice_jdn(year + 1));
    let nm = next11.index - month11.index;
    let (_, wm, wd) = ymd_from_jdn(winter_solstice_jdn(year), Calendar::Gregorian);
    println!("冬至落在 {year}-{wm:02}-{wd:02}");
    println!("本岁共 {nm} 个月{}", if nm == 13 { "，需要置闰" } else { "，为平年" });
    for i in 0..nm {
        let jd = beijing_jdn(new_moon(month11.index + i));
        let (yy, mm, dd) = ymd_from_jdn(jd, Calendar::Gregorian);
        // 岁内月序：0 = 十一月，之后按数序递增（闰月重复前一月名）
        let l = LunarDate::from_jdn(jd);
        println!(
            "  [{i:2}] {yy:04}-{mm:02}-{dd:02}  {}{}",
            l.month_name(),
            if l.leap { " (闰)" } else { "" }
        );
    }
    println!("下一个十一月起于日编号 {}", next11.jdn);
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        print_help();
        return;
    }
    let mut utc = false;
    let mut julian = false;
    let mut positional: Vec<String> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print_help();
                return;
            }
            "--utc" => utc = true,
            "--julian" => julian = true,
            "--year" | "-y" => {
                i += 1;
                let y: i64 = args.get(i).and_then(|s| s.parse().ok()).expect("需要年份");
                print_year(y);
                return;
            }
            "--sui" => {
                i += 1;
                let y: i64 = args.get(i).and_then(|s| s.parse().ok()).expect("需要年份");
                print_sui(y);
                return;
            }
            "--sanfu" => {
                i += 1;
                let y: i64 = args.get(i).and_then(|s| s.parse().ok()).expect("需要年份");
                print_sanfu(y);
                return;
            }
            "--terms" => {
                i += 1;
                let y: i64 = args.get(i).and_then(|s| s.parse().ok()).expect("需要年份");
                let lo = jdn_from_ymd(y, 1, 1, Calendar::Gregorian);
                for (j, name) in nongli::solar_terms_between(lo, lo + 370) {
                    let (yy, mm, dd) = ymd_from_jdn(j, Calendar::Gregorian);
                    println!("{yy:04}-{mm:02}-{dd:02}  {name}");
                }
                return;
            }
            a if a.starts_with("--") => {
                eprintln!("未知选项: {a}");
                std::process::exit(2);
            }
            a => positional.push(a.to_string()),
        }
        i += 1;
    }
    let cal = if julian {
        Calendar::Julian
    } else {
        Calendar::Gregorian
    };
    for s in &positional {
        let dt = parse_input(s).unwrap_or_else(|| {
            eprintln!("无法解析日期: {s}");
            std::process::exit(2);
        });
        let l = if utc {
            LunarDate::from_utc(dt)
        } else {
            LunarDate::from_datetime(dt, cal)
        };
        print_date(&l);
    }
}

/// `--help` 的正文；末尾的二十四节气表由 [`print_help`] 追加。
const HELP: &str = "nongli —— 中国农历 (GB/T 33661-2017)

用法:
  nongli <公历日期>[T时刻] ...   换算农历（默认北京时间）
  nongli --year  <年>            打印全年公历/农历对照与节气
  nongli --sui   <年>            打印该岁各月与置闰
  nongli --sanfu <年>            打印该年三伏
  nongli --terms <年>            打印全年二十四节气

选项:
  --utc     输入按 UTC 解释
  --julian  输入按儒略历解释（输出仍为公历，见下例）
  -h, --help  显示本帮助

示例:
  nongli 2026-02-17              农历、干支、生肖；当天有节气也会列出
  nongli 2024-04-04              当天节气：清明
  nongli 2026-02-17T23:59:59     带时刻，按北京时间
  nongli --utc 2026-02-17T17:00  按 UTC 输入（北京 02-18 01:00）
  nongli --julian 1500-01-01     按儒略历输入
  nongli 2026-02-17 2033-12-22   一次换算多个日期
  nongli -1974-06-01             负年份为天文纪年（0 = 公元前 1 年）
  nongli --year 2026             全年公历/农历对照，节气标 ★
  nongli --terms 2026            只列全年二十四节气
  nongli --sui 2033              该岁各月与置闰判定

二十四节气: ";

fn print_help() {
    print!("{HELP}");
    println!("{}", TERM_NAMES.join(" "));
}

