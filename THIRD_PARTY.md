# 第三方内容与许可

本仓库收录了以下第三方数据与系数。**它们各自适用其原许可**，声明保留在对应目录里。
这些还都是宽松许可（BSD-3 / MIT），所以本仓库整体选哪个许可都不受它们限制——
只要保住这些声明。

| 目录 | 内容 | 许可 |
|---|---|---|
| `third_party/erfa/` | IAU 2000A 章动系数（ERFA `nut00a.c`） | **BSD 3-clause**（NumFOCUS 基金会）|
| `third_party/elp2000-82b/` | ELP2000-82B 月球理论系数 | **MIT** |
| `third_party/vsop87/` | VSOP87D / VSOP87B 地球级数 | CDS 目录数据，引用 Bretagnon & Francou (1988) |
| `third_party/eclipse-catalog/` | 日食/月食目录 | NASA GSFC（Espenak & Meeus），仅用于验收 |

## 未收录的

* **GB/T 33661-2017《农历的编算和颁行》原文**（PDF 与 OCR 文本）——标准文本由 SAC
  持有版权，重新分发不合适。代码按条款号引用并实现（规则本身是事实）；
  需要原文请自行从标准发布渠道获取。见 `.gitignore`。
* **JPL 星历内核**（`kernels/de421.bsp`）——体积与许可都不适合进版本库，
  取法写在 `tools/truthcheck.py` 头部。
