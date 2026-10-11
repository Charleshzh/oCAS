# Rothstein–Trager 代数数根（二次扩域片段）算法规格

**版本**：oCAS 0.30.0（Wave 1 基层 + Wave 2 塔层/证书）
**状态**：**已实现**（2026-10，工作树）：§7 清单第 1/3/4/6 项完成；第 2 项
（`gcd_anf` 提 pub）按实现决策**不做**——基层直接用
`DenseUnivariatePolynomial<AlgebraicNumberField>` 的通用欧几里得 gcd
（`extended_gcd_poly`，`rt_algebraic.rs`），塔层用自写的 `KAPoly` gcd
（`logpart.rs`），`ocas-poly` API 无变更。第 5 项（ocas-tests 护栏 + 1892
回归）随 W3 收尾。
**实现偏差记录**（相对本文）：

1. 基层 `rothstein_trager_algebraic` 在 **`chain::parts_active()` 作用域内抑制**
   （`rational.rs`）：代数路径会把 parts 内层残项升级为含 `atan`/`log` 的完整答案，
   分部续算无法积分这些超越项而烧尽链预算（`rubi-01646` 回归，0.29.0 E2 机制在
   更强发射面上的复现）。抑制后恢复 0.29.0 的诚实残项形状。
2. 双二次四次的实根情形（`Δ₂ > 0`，`s±` 同号实数）**未实现**（`rt_algebraic.rs`
   只收 `Δ₂ < 0` 的复 4-循环；题库无需求证据）。
3. 实测发现并已修复的**既有潜伏错案**（与代数根无关）：`linear_form` 把 `x·x`
   读作 `1·x + 0`，幂法则对 `(x·x)²` 发射错误答案 `(x·x)³/3`（塔 proptest 的
   证书门随机命中）；修复为重复变量因子返回 `None`。

**前置规格**：[LOGPART_RT_SPEC_CN.md](LOGPART_RT_SPEC_CN.md)（0.29.0 的 ℚ 根片段；
本文是其「无理/虚常量根 → 0.30.0」拒绝分支的承接）
**核对来源**：

- SymPy 1.14.0：`ocas-tests/.venv/Lib/site-packages/sympy/integrals/rationaltools.py`
  （`log_to_atan` 279–325、`log_to_real` 343–445）、`integrals/risch.py`（`residue_reduce` 1258–1343）。
- oCAS 现状：`ocas-calc/src/integral/rational.rs`（`rothstein_trager` 424–475、
  `rational_roots` 512–554、`complete_square` 329–422）、`ocas-calc/src/integral/logpart.rs`、
  `ocas-calc/src/integral/certify.rs`、`ocas-poly/src/factor/algebraic.rs`（`gcd_anf` 273）、
  `ocas-poly/src/gcd.rs`（`gcd` 98）、`ocas-domain/src/algebraic.rs`（`AlgebraicNumberField` 388）。
- 本文全部数值例子（`R(z)`、根、`v`、完备性）已由 SymPy 1.14 精确计算核验（2026-10-01，
  Wave 0 探针，脚本过程见 §6）。

---

## 1. 问题与片段边界

0.29.0 的 Rothstein–Trager 对数部分只接受 `R(z)` 的**有理**根；`rational.rs:458-461`
与 `logpart.rs` 对「抽出全部有理一次因子后余因子非常量」一律诚实拒绝。被拒的典型族
（0.30.0 验收族，探针实测全部回退）：

| 族 | `R(z)`（SymPy 核验） | 根 | 所需扩张 |
|---|---|---|---|
| `1/(x²−2)` | `1−8z²` | `±√2/4` | ℚ(√2)，实 |
| `1/(1−3x²+x⁴)` | `(20z²−10z+1)(20z²+10z+1)` | `±1/4 ± √5/20` | ℚ(√5)，实 |
| `1/(1+x⁴)` | `256z⁴+1` | `(±√2 ± i√2)/8` | ℚ(√2, i)，复 |
| `(1+x²)/(1+4x²+x⁴)` | `4(24z²+1)²` | `±i√6/12`（二重） | ℚ(i√6)，复 |
| `1/(1−x⁴+x⁸)` | `(2304z⁴+1)²` | 四重复根 | ℚ(√6, i)，复 |
| `exp(x)/(1+exp(2x))`（塔层，`rubi-00174`） | `4z²+1` | `±i/2` | ℚ(i)，复 |

**0.30.0 片段**（划界即设计）：

1. `R(z)` 的不可约因子次数 ∈ {1, 2}，外加**双二次四次因子** `Az⁴+Bz²+C`
   （`z²` 经二次公式落 ℚ(√Δ)，再各开方得四个根；总扩张次数 ≤ 4）。
2. 其余（次数 3、非双二次四次、≥5）→ 继续诚实拒绝，**不引入 `RootSum` 记号**
   （SymPy 对 `exp(x)/(1+exp(2x))` 亦只给 `RootSum`；oCAS 的实值化 `atan(exp(x))`
   反而更强，但代价是显式根，这正是本片段的定位）。
3. 塔层（`logpart.rs`）只做**常数域根**：要求 `R(z)` 的各 `KElem` 系数有公共
   `KElem` 因子使商落入 `ℚ[z]`（逐系数 `ratio_constant` 对首项系数判定），
   之后与基层共用同一套根提取；0.29.0 的候选筛不变，仅把「完备性」从
   「ℚ 一次因子」推广为「{1,2} 次因子 + 双二次四次因子」。
4. 重根按不同值取一次（0.29.0 先例；`C`、`D` 族的平方因子先取无平方部分）。

## 2. 算法（基层 `ℚ(x)`，`rational.rs`）

### 2.1 流程

```
输入：简分数 a/d（d 首一无平方，gcd(a,d)=1，deg a < deg d；调用点已保证）
1. R(z) = resultant_x(a − z·d′, d) ∈ ℚ[z]            （现有插值路径，不动）
2. R 在 ℚ 上因式分解（现有 rational_roots 的 primitive.factor() 已给出全部不可约因子）
   - 一次因子 → 现有 ℚ 根路径（不动）
   - 二次因子 m(z)（判别式 Δ = B²−4AC）：
       Δ > 0：α = (−B+√Δ)/(2A)，§2.2 实根对
       Δ < 0：α = (−B+i√|Δ|)/(2A)，§2.3 复根对 → 实值化
   - 双二次四次因子 Az⁴+Bz²+C（B²−4AC 非平方）：§2.4
   - 其他 → 整个对数部分拒绝（与 0.29.0 相同的诚实回退；
     不做「部分因子积分 + 残项」——残项份额本身就要在扩域里算 gcd，半套输出无意义）
3. 每个被接受的因子：扩域 gcd（§2.5）得 v_α ∈ ℚ(α)[x]，Πv == d 与留数恒等式自证（§2.6）
4. 发射：实根对 → 带 √Δ 字面量的 log 项；复根对 → log(u²+w²) 与 atan 项（§2.3 公式）
```

### 2.2 实二次因子（Δ > 0）

两根 `α± = (−B±√Δ)/(2A)`。各算 `v± = gcd(a − α±·d′, d) ∈ ℚ(√Δ)[x]`，发射

```
α+·log(v+) + α−·log(v−)
```

`√Δ` 字面量沿用现有 `sqrt_positive_rational` 的 `√(p·q)/q` 形状
（`rational.rs:563-579`）。`v±` 的次数 = 共享留数的 `d` 根数（一般 1，
>1 时按 α 的多项式系数原样发射，公式不变）。

**核验例（B 族）**：`1/(1−3x²+x⁴)`，`c = 1/4−√5/20`，
`v = x − (√5+1)/2`（SymPy `gcd(..., extension=√5)` 核验），四根 `Πv == d` ✓。

### 2.3 复二次因子（Δ < 0）：共轭对实值化

根 `α = p + qi`（`p = −B/(2A)`、`q = √|Δ|/(2A)`，`q > 0`），共轭 `ᾱ = p − qi`。
算 `v = gcd(a − α·d′, d) ∈ ℚ(α)[x]`，把 `v` 的系数按 `α ↦ p + qi` 拆出实虚部：
`v(x) = u(x) + i·w(x)`，`u, w ∈ ℚ(√…)[x]`（二次情形下 `u, w ∈ ℚ[x]`，
因为 `v` 的系数是 `a₀ + a₁α` 形，`a₀, a₁ ∈ ℚ`）。共轭对合并：

```
α·log(u + i·w) + ᾱ·log(u − i·w) = p·log(u² + w²) − 2q·atan(w/u)
```

（`log(u+iw) = ½log(u²+w²) + i·atan(w/u)`，取实部即得；atan 参数 `w/u` 的
分支差一个逐区间常量，形式导数精确——这是微分域证书的语义范围，实值分支
说明按既有约定写入文档。）

**核验例（C 族）**：`(1+x²)/(1+4x²+x⁴)`，`R = 4(24z²+1)²` 取无平方因子
`24z²+1`，`α = i√6/12`（`p = 0`、`q = √6/12`）：`v = gcd(a − α·d′, d)`，
实值化后 `u²+w²` 为常数（log 项消失），只剩 `atan` 项——与 SymPy 的
纯 atan 答案 `√6(2atan(√6x/6) + 2atan(√6x³/6+5√6x/6))/12` 同族
（atan 参数形状可差逐区间常量，证书判定等价性）。

### 2.4 双二次四次因子 `Az⁴ + Bz² + C`

`z² = s± = (−B±√Δ₂)/(2A)`，`Δ₂ = B²−4AC`。四个根 `z = ±√s±`：

- `s±` 一正一负或同号实数 → 根为两实两虚或全实，按 §2.2/§2.3 逐根处理；
- `s±` 共轭复数（`Δ₂ < 0`，如 `256z⁴+1`：`s± = ±i/16`）→
  `√s± = (u ± iv)` 显式形：`u = √((|s|+Re s)/2)`、`v = √((|s|−Re s)/2)`
  （`s = ±i/16`：`u = v = √2/8`）。四个根 `±u±iv` 配成两对共轭
  `(u+iv, u−iv)`、`(−u+iv, −u−iv)`。

**gcd 只算一次**：添加四次域 ℚ(α)（`α` 为 `Az⁴+Bz²+C` 任一根），
`v_α = gcd(a − α·d′, d) ∈ ℚ(α)[x]`（单次扩域 gcd，deg α = 4），
再把 `α` 的四个显式根形逐个代入 `v_α` 的系数取实虚部、按对实值化。
这避免嵌套域 `AlgebraicExtension<AlgebraicNumberField>` 的欧几里得链。

**核验例（A 族）**：`1/(1+x⁴)`，`α = √2(1+i)/8`，
`v_α = x + √2/2 + i√2/2`（SymPy `extension=[√2, I]` 核验），
四根 `Πv == d` ✓；两对实值化给出
`(√2/8)·log((x+√2/2)²+1/2) − (√2/8)·log((x−√2/2)²+1/2) − (√2/4)·atan((√2/2)/(x+√2/2)) − (√2/4)·atan((√2/2)/(x−√2/2))`，
与 SymPy 的 `√2·log(x²∓√2x+1)/8`、`√2·atan(√2x±1)/4` 逐项等价
（log 参数恒等；atan 参数差逐区间常量）。

**D 族**（`1/(1−x⁴+x⁸)`，`R = (2304z⁴+1)²`）同构：无平方四次因子
`2304z⁴+1 = (48z²−4√6z+1)(48z²+4√6z+1)`（ℚ(√6) 上），复根实值化后
log+atan 混合形（与 SymPy 的 √6 log/atan 答案同族）。

### 2.5 扩域 gcd 的落地

`v = gcd(a − α·d′, d)`，全部系数在 `K = ℚ(α)`：

- 资产：`DenseUnivariatePolynomial<AlgebraicNumberField>` 的域运算现成
  （`AlgebraicExtension<RationalDomain>` 已实现 `Domain` + `EuclideanDomain`，
  `ocas-domain/src/algebraic.rs:300-385`；`ocas-poly` 的 Trager 路径
  `factor/algebraic.rs` 已在该类型上跑通分解，0.17.0 起）。
- **不用** `gcd.rs:98` 的子结式 PRS：其 `.expect("subresultant … division is
  exact")` 的整环精确除法语义在域上是浪费；域上普通带余除法即精确，
  用 `dense.rs:522` 的 `extended_gcd_poly` 取首一 gcd（与 `logpart.rs`
  的 `KPoly::resultant`「域上 Euclidean 递归」同一原则，LOGPART 规格 §4.2）。
- `factor/algebraic.rs:273` 的 `gcd_anf` 现为 `pub(crate)`：0.30.0 将其
  提升为 `pub`（或经 `factor` 模块 re-export），供 `ocas-calc` 复用——
  这是「复用 Trager 资产」路线图条目的具体落点。
- 构造输入：`a − α·d′` 的系数 = `aᵢ − α·d′ᵢ`（`aᵢ, d′ᵢ ∈ ℚ`），
  即 `AlgebraicElement` 的二系数形 `aᵢ + (−1)·d′ᵢ·α`。

### 2.6 自证闭环（与 0.29.0 §1.4 同构，扩到扩域）

每个不可约因子 m 处理后验证（任一失败 → 整个对数部分拒绝）：

```
(i)   Π_{m 的各根} v_α ·（其余因子的 v 积）== d       （首一，逐系数）
(ii)  a·Πv == d·Σ c·Dv·(Π/v)                          （交叉相乘，ℚ(α) 内精确）
(iii) 实值化后的 atan/log 发射形式另行满足：
      顶层 tower_difference_is_zero / 证书层（§4）终检
```

(ii) 在 ℚ(α) 上验证即以 α 的多项式系数逐位比对；实值化等价性不由
(ii) 覆盖（它在 ℂ 上成立），由数值 oracle 与证书层兜住——
**这正是 D2（证书根式层）必须先于/随 Wave 1 落地的原因**。

## 3. 塔层（`logpart.rs`）

1. `R(z)` 的 `KElem` 系数判定：全体系数对首项系数做 `ratio_constant`
   （`tower/elem.rs:27-48`）——全为有理数比则归约为 `ℚ[z]` 情形，
   走 §2 的同一根提取；否则拒绝（非常量系数或符号常量，不在本片段）。
2. `v = gcd(a1 − α·D d1, d1)` 需要 `k(α)[t]` 上的 gcd：新增
   `KAlphaElem = Vec<KElem>`（`α` 的次数 < deg m 的系数向量）+
   逐分量算术 mod m(α)（deg m ≤ 2 先行；双二次四次塔层无题库需求，不做）。
   **复二次 m(z) = Az²+C（B=0）** 是 `α = i·r`（`r ∈ ℚ(√…)`）的特例——
   `rubi-00174`（`R = 4z²+1`）落在这里：`v = t + i`（SymPy 核验），
   实值化 `−atan(1/t)`（与 `atan(t) − π/2` 差常量，形式导数精确）。
3. 实值化公式同 §2.3；`u, w` 为 `KElem` → atom 走现有 `KElem` 发射路径。

## 4. 证书根式层（`certify.rs` W3，`CertMethod::RadicalField` 首次启用）

**动机**：Wave 1/2 的答案含 `√Δ` 常量字面量与 `log(u²+w²)`；现有初等域层把
常量根式当独立生成元（`certify.rs:548-552`），看不到 `√2·√2 = 2`，
`D(F)−f` 无法判零 → 新答案将全体落入 `notinfield`/`budget`，
`certified_rate` 不进反退。根式层是本波的**硬前置**（计划 D2）。

**设计**（与模块文档 21–24 行的预留一致）：

```
1. 收集 D(F)−f 中全体非常量基根式 Sᵢ^(pᵢ/qᵢ)（含常量基根式 2^(1/2)——
   常量情形现有层已添加为生成元但无关系，根式层为其补关系）；
2. 每根式添加生成元 rᵢ，关系 rᵢ^qᵢ = Sᵢ；qᵢ 取各 pᵢ/qᵢ 的分母；
   同一基的不同分数幂共用同一 rᵢ（pᵢ>1 为 rᵢ^pᵢ）；
3. 关系集对各自 rᵢ 首一 => 三角列：把差表达式转换为以
   {x, 原子, 常量, rᵢ} 为变元的 Sparse 有理函数后，对分子按任一固定
   顺序逐 rᵢ 取模余式（单变元首一关系，普通多项式余式即得正规形，
   无需完整 Gröbner）；
4. 判零：分子正规形 == 0 且分母正规形 != 0。
```

可靠性：关系 `rᵢ^qᵢ = Sᵢ` 是精确恒等式，约化只会把真零判零、
绝不把非零判零（假阴性只来自预算）。预算：`MAX_CERT_GENS` 同样约束
`rᵢ` 个数；`qᵢ ≤ 8`、根式基节点数上限随层记录。

## 5. 接入点与输出形状

| 位置 | 改动 |
|---|---|
| `rational.rs::rothstein_trager` | `fully_split == false` 分支（458–461）改为：收集非一次不可约因子，逐因子走 §2.2–2.4；任何因子越界 → 原样诚实回退 |
| `rational.rs::log_part` | 不动（Case 1/2 快路径优先；RT 只对 deg d ≥ 3 进入） |
| `logpart.rs::rothstein_trager` | ℚ 根完备性判定泛化 + §3 塔层 `KAlphaElem` gcd |
| `certify.rs` | 新增根式层；`CertMethod::RadicalField` 开始产出（绑定层 `"radical"` 字符串已存在，无需 API 变更） |
| `ocas-poly/src/factor/algebraic.rs` | `gcd_anf` 提升为 `pub` |

输出形状：新答案含 `d^(1/2)` 常量字面量、`log(u²+w²)`、`atan(w/u)`；
均不自匹配 `rational`/`risch` 入口模式（幂等不变式），
其中 `atan` 实形保持实值（不引入 `I`——`build.rs:124` 的 `I` 拒绝不变）。

## 6. 测试例（全部经 SymPy 1.14 精确核验，2026-10-01）

| # | 输入 | 期望 | 关键中间量 |
|---|---|---|---|
| 1 | `1/(x²−2)` | 已解（Case 2，回归保护） | — |
| 2 | `1/(1−3x²+x⁴)` | 四个带 √5 的 log | `R=(20z²±10z+1)`；`c=1/4−√5/20`，`v=x−(√5+1)/2` |
| 3 | `1/(1+x⁴)` | 两 log + 两 atan（√2） | `R=256z⁴+1`；`α=√2(1+i)/8`，`v=x+√2/2+i√2/2`；`Πv==d` ✓ |
| 4 | `(1+x²)/(1+4x²+x⁴)` | 纯 atan（√6） | `R=4(24z²+1)²`；`p=0` ⇒ log 项消失 |
| 5 | `1/(1−x⁴+x⁸)` | log+atan（√6） | `R=(2304z⁴+1)²`；`2304z⁴+1=(48z²±4√6z+1)` |
| 6 | `1/(x⁴+1)²` | 有理部 + 3/4·(例 3) | Hermite 现有；残项链内解析现有 |
| 7 | `exp(x)/(1+exp(2x))`（塔层） | `atan(exp(x))` 形 | `R=4z²+1`；`v=t+i`；实值化 `−atan(1/t)` |
| 8 | `1/(x³+x+1)` 类（R 含三次不可约因子） | **拒绝**（诚实回退） | 片段外 |
| 9 | 探针回归：`1/(x²+5)`、`1/(x⁶(5+x²))`、`1/(x²−1)²` 等 | 不变 | 现有测试全绿 |

每例三验：证书（结构/域/根式层之一）+ 数值 oracle + SymPy 对拍。

## 7. 实现检查清单（0.30.0 Wave 1–2）

1. `certify.rs` 根式层（§4）+ 单测（含 `1/sqrt(1-x²)`→`asin(x)` 现拒绝例的放行）。
2. `ocas-poly`：`gcd_anf` 提 `pub` + rustdoc 示例。
3. `rational.rs`：§2.1 流程 + 二次/双二次根提取 + 实值化 + §2.6 自证。
4. `logpart.rs`：§3 塔层常数根推广 + `KAlphaElem`（deg ≤ 2）。
5. 测试：§6 全表 + 历史错案护栏（`integral_certify.rs` 增补）+ 1892 探针回归。
6. 明确不做：`RootSum` 记号；非双二次四次与 ≥3 次因子；塔层双二次；
   嵌套域 `AlgebraicExtension<AlgebraicExtension<…>>`；符号常量塔
   （Wave 0 决策：与本波扩域机械**不同构**——那 60 例是 Const(k) 推理
   问题而非数值代数数问题，**0.30.0 no-go**，另立项评估）。
