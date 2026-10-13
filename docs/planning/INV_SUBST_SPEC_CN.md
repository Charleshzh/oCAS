# 反函数代换引擎（inverse-substitution engine）算法规格

**版本**：oCAS 0.30.0（Wave 5）
**状态**：实现规格，不修改任何 Rust 代码（`ocas-calc/` 下的 W4 实现由另一分支
并行进行，本文只写接口约定）
**前置规格**：[MULTI_RADICAL_SPEC_CN.md](MULTI_RADICAL_SPEC_CN.md)（W4 多根式基
归约；§2.7 的接口约定与其配套）、[ALGEBRAIC_RT_SPEC_CN.md](ALGEBRAIC_RT_SPEC_CN.md)
（诚实拒绝与自证闭环体例、parts 作用域教训）
**核对来源**：

- 归因数据：`ocas-tests/data/invtrig_attribution_030.csv`（Rubi 1892 题库
  `inverse-trig-hyper` 桶 **138 例未解**，基线 8/146 解；分类器见
  `ocas-tests/scripts/attribute_030.py:129-156`——regex 分桶，量级参考，
  边界噪声在 §1.2 注明）。
- 路线图：`docs/planning/ROADMAP_CN.md` §0.30.0「反函数代换引擎：多项式权 ×
  `(a+b·f(ax+b))^k`（Rubi 4.x–5.x 思路的定理化子集）」；
  `GAP_ANALYSIS_CN.md` 0.27.3 段（该桶 94.5% 未解归因于「无此引擎」）。
- oCAS 现状（能力与接口均以源码为准）：
  `ocas-calc/src/integral/inverse_trig.rs`（2405 行，0.27.1 起：M1 核导数幂法则
  `(a+b·F(σx))^k·C·R^e`、`1 ≤ k ≤ MAX_POWER_K=6`；M2 inv-hyp 代换 `t=F(arg)`，
  仅 `asinh/acosh`、`1 ≤ k ≤ MAX_SUBST_K=4`、多项式权总次数 ≤ 4、至多一个
  伴随根式 `(s·E)^(m/2)`、`m∈{1,3}`、仅零截距参数；M3 仿射参数裸形；
  P1 逆复合消去、P2 分支平移核、P3 伴随被开方有理幂；预算 `MAX_NODES=200`/
  `MAX_SUBST_NODES=300`/`MAX_REDUCED_NODES=600`/`MAX_HYP_EXP=8`；
  入口 `integrate_inverse_trig`，2405 行文件内 `InvFun` 六头表、
  `family_kernel_of_u`/`family_kernel_exp`/`family_sign` 族表）、
  `exp_log.rs` E4（`exp(k·F(u))` 代数化，机制无关但共享族判定）、
  `halfpower.rs`（片层因子 `σ` 先例：形式导数为零的逐区间常量，模块文档
  89–99 行）、`binomial.rs`（代换→回积分→回代→残留拒绝的流水线先例）、
  `certify.rs`（结构/初等域/根式三层证书；`MAX_CERT_NODES=20_000`、
  `MAX_CERT_GENS=16`）、`chain.rs`（`parts_active()` 作用域、
  `MAX_CHAIN_ENTRIES=256`、`MAX_RESIDUAL_ENTRIES=128`）、`mod.rs` 链阶段序
  （`… sqrt_quadratic → binomial → exp_log → inverse_trig → trig_reduce →
  trig_reduction → trig_kernel → expand_retry → halfpower → elliptic →
  heuristic`；W4 的 `multi_radical` 将插在 `binomial` 与 `exp_log` 之间）。
- 本文全部数值例子已由 SymPy 1.14 精确核验（2026-11，`ocas-tests/.venv`，
  `uv run python`；核验方式逐例标注于 §5）。

---

## 1. 问题与片段边界

### 1.1 问题

被积函数含反函数头 `F(w)`（`F ∈ {asin, acos, atan, asinh, acosh, atanh,
acsc, asec, acoth, asech, acsch}`，`w` 仿射）乘以多项式/根式权。现有
`inverse_trig.rs` 只覆盖三个狭窄片段：M1 要求**导数形状的核因子**精确在场
（`(a+b·F(σx))^k·C·E(σx)^e`，其余因子一律拒绝）；M2 只做 `asinh/acosh` 两族
且多项式权 ≤ 4 次；M3 只做无权的 `F(a+b·x)` 裸形。于是
`x²·asin(c·x)`、`(a+b·asin(c·x))²`、`x³·asin(c·x)/√(1−c²x²)` 这类 Rubi 4.x–5.x
的主力族整体回退——138 例、占桶 94.5%，是 0.30.0 最大的单机制缺口。

### 1.2 138 例的形状分布（子类×特征交叉表）

`class` 列为多标签（`族:特征+特征…`；族按首个反函数头分 `invtrig`/`invhyp`）。
特征计数（一例多标签，合计 > 138）：

| 特征 | 含义 | invtrig | invhyp | 合计 |
|---|---|---|---|---|
| `zero-int-arg` | 参数为 `c·x`（零截距） | 57 | 44 | **101** |
| `kernel-power` | `(a+b·F)^k` 核幂在场 | 23 | 29 | **52** |
| `radical` | 伴随根式/分数幂在场 | 29 | 18 | **47** |
| `bare-power` | 裸 `F^k`（k≠1） | 14 | 14 | **28** |
| `affine-arg` | 参数为 `c+d·x`（带截距） | 10 | 13 | **23** |
| `F-in-denom` | F 出现在分母 | 6 | 7 | **13** |
| `nonlinear-arg` | 参数含 `x^n`（n≥2） | 6 | 6 | **12** |
| `bare` | 无特征（指数参数等） | 0 | 2 | **2** |

完整子类表（31 类，按数量降序；题号去 `rubi-` 前缀，列前 8 个）：

| 数量 | 子类 | 代表题号 | 归属（§1.2 末尾） |
|---|---|---|---|
| 15 | invtrig:zero-int-arg | 00167, 00244, 00306, 00441, 00489, 00637, 00642, 00709 | W5a（00441/01786 属 W5c 倒数头；00244 含非伴随分母，见下） |
| 11 | invtrig:radical+zero-int-arg | 00091, 00228, 00269, 00459, 00550, 00747, 00821, 00950 | W5b（伴随）/拒绝（非伴随） |
| 9 | invtrig:kernel-power+radical+zero-int-arg | 00062, 00320, 00553, 00581, 00873, 01472, 01563, 01666 | W5b；00062/00873/01563 另涉 W4 双仿射基（§2.7） |
| 9 | invhyp:zero-int-arg | 00071, 00251, 00529, 00535, 00635, 01317, 01376, 01585 | W5a（00251 F 在分母→拒绝；00635/01317/01376 倒数头→W5c） |
| 9 | invhyp:kernel-power+zero-int-arg | 00081, 00172, 00592, 01225, 01259, 01467, 01487, 01677 | W5a（00081/01225/01259/01487 分母不配对→polylog 拒绝，§2.5） |
| 8 | invhyp:bare-power+zero-int-arg | 00164, 00308, 00465, 01274, 01347, 01500, 01624, 01745 | 负幂/高幂为主→多为拒绝（§2.5） |
| 7 | invtrig:bare-power+radical+zero-int-arg | 00149, 00526, 00700, 00811, 01236, 01242, 01402 | W5b（00811/01402）与拒绝（00149/01236/01242 负半幂） |
| 7 | invhyp:kernel-power+affine-arg | 00181, 00225, 00229, 00372, 00781, 00886, 01044 | W5a（00372/00781 倒数头→W5c） |
| 6 | invtrig:bare-power+zero-int-arg | 00193, 00554, 00555, 00688, 01692, 01839 | W5a（00554/00555）与拒绝（01692 符号幂 n） |
| 6 | invtrig:kernel-power+zero-int-arg | 01094, 01132, 01148, 01205, 01218, 01518 | W5a（01132/01218 倒数头→W5c） |
| 5 | invhyp:bare-power+radical+zero-int-arg | 00299, 00610, 00720, 01409, 01571 | W5b（01409）与拒绝（负半幂族） |
| 5 | invhyp:radical+zero-int-arg | 00408, 00544, 00689, 01304, 01733 | W5b（01733 核验例 CR3）；00408 含 `i`→拒绝 |
| 5 | invtrig:affine-arg | 00511, 00782, 00951, 01564, 01708 | W5a/W5c（00511/01708 倒数头核验例 A4；00782/01564 倒数参数 `c/x`） |
| 4 | invtrig:nonlinear-arg | 00152, 00444, 00761, 01886 | 拒绝（§7）；部分经 `u=x²` 预代换后重评（W5c 可选） |
| 4 | invhyp:kernel-power+radical+zero-int-arg | 00426, 00569, 00837, 01830 | 符号幂 m/n→拒绝（预算外形状，§7） |
| 3 | invhyp:kernel-power+nonlinear-arg | 00273, 00383, 00850 | 拒绝 |
| 3 | invtrig:kernel-power+affine-arg | 00351, 01643, 01746 | W5a（01643 的 k=5/2→拒绝，核验例 R3） |
| 2 | invhyp:affine-arg | 00026, 00531 | 00531 核验例 A3→polylog 拒绝；00026→W5c |
| 2 | invhyp:kernel-power+F-in-denom+zero-int-arg | 00088, 01028 | 拒绝（核验例 R1 族） |
| 2 | invtrig:kernel-power+F-in-denom+zero-int-arg | 00364, 01301 | 拒绝 |
| 2 | invhyp:kernel-power+F-in-denom+affine-arg | 00431, 00458 | 拒绝 |
| 2 | invhyp:bare | 00717, 01681 | 拒绝（指数参数→polylog，0.31.0） |
| 2 | invtrig:kernel-power+F-in-denom+nonlinear-arg | 01702, 01730 | 拒绝 |
| 1×9 | 其余 9 类各 1 例 | 00205, 00302, 00425, 00486, 01184, 01194, 01207, 01322, 01388, 01460 | 逐类见归属：00425→W5b；00486 完全平方基→P 类坍缩后重评；余多为拒绝 |
| **138** | | | |

**标签噪声声明**（分类器为 regex，与 `attribute_030.py` 头注一致）：
`F-in-denom` 只识别 F 落在 `/(...)` 括号内的写法，`x^4/((1−a²x²)^4·atanh(a·x))`
（00251）这类 F 在分母乘积因子里的漏标——特征计数是下界；归属列已逐例复核
修正明显错标。**归属汇总**：

| 归属 | 含义 | 例数（估） |
|---|---|---|
| **W5a** 零截距/仿射代换族（无根式，§2.1–2.5） | 6 主头、`k` 整数、多项式/配对有理权 | 候选 ~45，预估解锁 **30–40** |
| **W5b** 伴随根式吸收（§2.2） | 根式全为族核 `s·E(w)` 的半幂 | 候选 ~30，预估解锁 **15–25** |
| **W5c** 可选扩展 | 5 个倒数头（~18 例涉）、`u=x²` 预代换子集 | 不计入验收数 |
| **拒绝** | F-in-denom（13）、nonlinear-arg（12）、bare（2）、符号幂、含 `i`、polylog 需求 | ~45，诚实回退 |

### 1.3 片段总原则（划界即设计）

1. **定理覆盖的参数化族**：引擎的每一条接受路径都是一个定理族（「对任意
   关于 `x` 的多项式权与族核半幂，`u=F(w)` 把被积函数化为 `P(u)·T(u)`，
   `T` 为三角/双曲多项式」），**禁止逐例规则枚举**——这是与 Rubi 路线
   划清界限的核心：Rubi 4.x–5.x 把该族展开成数百条按 `m,k,e` 分案的规则，
   本引擎只实现代换 + 归约 + 回积分三段通用机械，初等性由回积分的
   残项检查**涌现**而非逐案预判。
2. **诚实拒绝**：前提不满足、t-形回积分残留 `Integral`、证书不通过——任一
   发生即返回 `None`，链继续下行；绝不输出部分猜测。
3. **确定性预算**：新机制自带预算常量（§3.2），单阶段总工作量有界；代换
   扫描与回积分按节点记账（防「扫描免费、回积分爆炸」）。
4. **输出必须过符号证书**：`certify.rs` 三层（结构/初等域/根式）。回代后
   答案只含原始头（`F(w)`、根式、`log`），`D(F)−f` 在初等域层对
   `asin` 等头列独立生成元即可判零；伴随根式归根式层。证书拒绝 ⇒ 不得发射。
5. **幂等不变式**：输出不再匹配本模块入口模式——答案中 `F(w)` 只出现在
   与原式**同参数**的位置且幂次不增（核幂 `k` 的答案含 `F^{k−1}` 以下与
   `F^k·x` 形；§2.4），重入链时 M1 先 claim 同族形状，不产生循环
   （`MAX_CHAIN_ENTRIES` 内终止，用回归测试钉死）。
6. **`chain::parts_active()` 作用域纪律**（`rubi-01646` 教训，ALGEBRAIC_RT
   实现偏差记录第 1 条）：本机制在 parts 作用域内只准输出完整答案或
   `None`；代换要么整体成功要么整体拒绝，不把残项升级为含 `atan`/`log`
   的完整答案。
7. **实值分支语义单独文档化**（§4）：主值区间、片层因子、定义域门全部
   显式记录；不引入 `I`（`tower/build.rs:124` 不变式）。
8. **不改动 Risch 塔**：全部工作在原子层完成；t-形回积分经
   `integrate_raw` 回链，复用 `trig_reduce`/`trig_reduction` 既有能力，
   不新增塔生成元。

## 2. 算法设计

### 2.1 代换构造（`u = F(w)` 型，核心定理族）

**定理（反函数代换，三角/双曲有理化）**：设 `F` 为下表之头，`w = ρ+σ·x`
（`σ ≢ 0`，符号允许），`G = F⁻¹`，则被积函数中

```
x → (G(u) − ρ)/σ,        dx → G′(u)/σ · du,        F(w) → u,
(a + b·F(w))^k → (a + b·u)^k        （k 整数时仍为 u 的多项式）
```

| F | `x=G(u)`（零截距 `w=σx`） | `dx/du`（×σ⁻¹） | 族核 `E(w)` | `E(w)^(m/2) →` |
|---|---|---|---|---|
| asin | `sin(u)/σ` | `cos u` | `1−σ²x²` | `cos^m u` |
| acos | `cos(u)/σ` | `−sin u` | `1−σ²x²` | `sin^m u` |
| atan | `tan(u)/σ` | `sec²u` | `1+σ²x²` | `sec^{2m} u` |
| asinh | `sinh(u)/σ` | `cosh u` | `1+σ²x²` | `cosh^m u` |
| acosh | `cosh(u)/σ` | `sinh u` | `σ²x²−1` | `sinh^m u` |
| atanh | `tanh(u)/σ` | `sech²u` | `1−σ²x²` | `sech^{2m} u` |
| acsc（W5c） | `csc(u)/σ` | `−csc u·cot u` | `σ²x²−1` | `cot^m u·σ` |
| asec（W5c） | `sec(u)/σ` | `sec u·tan u` | `σ²x²−1` | `tan^m u·σ` |
| acoth（W5c） | `coth(u)/σ` | `−csch²u` | `1−σ²x²` | `csch^{2m} u·σ` |
| asech（W5c） | `sech(u)/σ` | `−sech u·tanh u` | `1−σ²x²` | `tanh^m u·σ` |
| acsch（W5c） | `csch(u)/σ` | `−csch u·coth u` | `1+σ²x²` | `coth^m u·σ` |

（一般仿射 `w=ρ+σx`：`x=(G(u)−ρ)/σ` 展开后 `G(u)` 的多项式权不变质。）
族核列即 M1 既有 `family_kernel_of_u` 的同一六族表（`E(u)`：`asin/acos/atanh`
为 `1−u²`，`atan/asinh` 为 `1+u²`，`acosh` 为 `u²−1`），**复用不复制**。
`E(w)^{m/2} = s^{m/2}·(E 在 x=G(u) 处的值)^{m/2}` 化为 `G′(u)` 的幂乘**片层
因子 σ**（§4.2），符号归一化 `s` 沿用 inverse_trig 模块文档 61–64 行的
既有约定（数值 `s` 必须为正，符号 `s` 假定正）。

**初等性边界的诚实刻画**（不做预判，只列后果）：t-形为
`P(u)·T(u)`，`T` 为 `sin/cos`（三角族）或 `sinh/cosh`（双曲族）的有理式。
多角归约后逐项为 `u^j·H(μ·u)`（`H ∈ {sin,cos,sinh,cosh}`）——`j` 次分部
初等——或残留 `u^j·tan u`、`u^j·sech u`、`u^j·log(sin/cosh u)` 等
**polylog 形**（0.31.0 领土）。后者经 `integrate_raw` 回积分必然残留
`Integral`，机制据此拒绝——**初等性由残项检查涌现**，不为它写预判规则
（核验例 R1/A3/K4：SymPy 对 polylog 形同样回退，§5）。

### 2.2 导数因子吸收与伴随根式（W5b 的 47 例根式）

- 分子权：`x^m → ((G(u)−ρ)/σ)^m`，展开为 `G(u)` 的多项式
  （`m ≤ MAX_IS_POLY_DEG`，预算）。
- **`dx` 的 `G′(u)` 因子是唯一自动到场的导数因子**——这正是代换相对 M1
  的增益：M1 要求被积函数**自带**精确核因子 `E^e`，代换族不需要。
- 伴随根式：`(s·E(w))^{m/2}`（`m` 奇偶任意，`|m| ≤ MAX_IS_RAD_M`）化为
  `G′(u)^m·σ`；**非伴随根式一律拒绝**（如 `atan` 族的 `sqrt(1+u²)` 半幂、
  与 `w` 无关的 `(d+e·x²)^{m/2}`——后者代换后为 `√(d+e·G²)`，是椭圆领土，
  诚实回退）。伴随判定复用 P3 的 `radicand_scale` 机制（同一参数 `w`、
  同一族核、常量归一化商）。
- 与 M1 的关系：M1 是本族 `T(u)=1`（恰好只剩 `(a+bu)^k`）的退化情形，
  保持先 claim（成本更低、答案形状更规范）；本机制只在 M1/P3 拒绝后运行。

### 2.3 t-形回积分与多角归约

```
1. 代换得 t-形 P(u)·T(u)（节点 ≤ MAX_IS_SUBST_NODES）；
2. 三角族：T 中 sin/cos 幂积经 trig_reduce_products 化多角和
   （复用，预算内）；tan/sec 偶幂经 tan²=sec²−1 归约为 sec 幂多项式；
   双曲族：复用 M2 既有 cosh^A·sinh^B → Σc_m·H(m·u) 指数二项展开
   （MAX_HYP_EXP=8 不变）；
3. 展开项经 integrate_expanded_terms / integrate_raw 回积分：
   - u^j·H(μu) 族由链内 parts（深度 ≤ 2）+ M2 同款表式闭式收尾
     （三角类比：trig_reduction T3 的 u^j·sin/cos 表）；
   - 任何 Integral 残留 → 整体 None（polylog 形在此被诚实拦下）；
4. contains_integral 终检 + 幂等形状检查（§1.3-5）。
```

### 2.4 回代

`u → F(w)` 逐节点替换（`replace_symbol`/`replace_atom` 先例）。答案形状：

- `x·(a+b·F(w))^k`、`(a+b·F(w))^{k−j}`（j≥1）、`F(w)` 一次；
- 伴随根式 `E(w)^{1/2}` 字面量（如 `√(1−c²x²)`，来自 `cos u`/`sin u`
  的回代：`cos(F(w)) = √E(w)`——仅当 `m` 为奇数时需要，且必须乘回片层
  因子 σ，§4.2）；
- `log`、`atan`（t-形有理残项经有理后端产生，如 `atan(tanh(u/2))` 形
  回代为 `atan(tanh(F(w)/2))`）。
- **禁止** `F(G(·))` 复合残留（`asin(sin u)` 形）——回代必须在发射前
  完成，这是 P1 既有护栏的同一纪律（§8 风险 2）。

### 2.5 拒绝条件（任一即 `None`）

1. 参数 `w` 非仿射（`linear_form` 拒绝；重复变量因子返回 `None` 的
   0.30.0 已修复语义）；`σ ≢ 0` 判零失败（`collect_terms` 折叠后判不了
   即拒绝，诚实）。
2. `k` 非整数或 `k ∉ [−0, MAX_IS_K]`（负幂/F-in-denom：核验例 R1/R2/R3
   证明其 t-形为 `P(u)^{−k}·T(u)`——`u` 幂在分母，非初等族，0.31.0
   polylog 评估）。
3. 多项式权次数 > `MAX_IS_POLY_DEG`、节点/预算超支。
4. 非伴随根式、多个不同基根式（→ §2.7 接口）、根式含 `i`（00408）。
5. 符号幂（`m`/`n`/`p` 指数：00426/00569/00734/01692 等 ~10 例）——
   与 0.29.0 exp-log 桶符号指数同一纪律：不可判定即拒绝。
6. t-形回积分残留 `Integral`（polylog/椭圆涌现于此）。
7. 证书三层全部拒绝。

### 2.6 与 `inverse_trig.rs` 既有机制的分工（扩展 vs 新模块）

**建议：扩展 `inverse_trig.rs`，新增机制 M4（反函数代换），不新建模块。**
理由：

1. **共享族表**：M4 的代换表（§2.1）与 M1 的核/指数/符号表、M2 的
   `G`/`G′` 构造是同一组数据的三个侧面；新模块要么复制要么把
   `InvFun`/`family_*` 提为 `pub(crate)` 造成不必要的 API 面扩张。
2. **链位置已正确**：`inverse_trig` 阶段在 `exp_log` 之后、`trig_reduce`
   之前——正是 M4 应在的位置（先于 heuristic 防分部烧预算，后于初等
   引擎保既有答案形状的 first claim）。新模块需要在 `mod.rs` 新挂一个
   相邻阶段，徒增链序推理负担。
3. **M2 是先例**：M2 本身就是「受限的反函数代换」（asinh/acosh 两族），
   M4 是它向 6（11）头的推广；模块文档的 M1–M3+P1–P3 编号体系直接续
   M4，认知连续。
4. **回积分基建共享**：`integrate_expanded_terms`、`MAX_CANCEL_BUDGET`
   记账、`pick_subst_symbol` 均在本文件/同模块可达。

**代价与缓解**：文件已 2405 行；M4 预计 +600–800 行。缓解：M4 代码以
独立 section 组织（同既有 `// ====` 分节惯例），模块文档头部加 M4 段；
若实测超过 ~3200 行再拆 `inverse_trig/subst.rs` 子模块（纯机械移动，
不影响链）。`integrate_inverse_trig` 的分派序改为
`P1 → M3 → M1 → P3 → P2 → M2 → **M4**`（M4 最通用、最贵，放最后）。

### 2.7 与 W4 `multi_radical.rs` 的接口约定（另一 agent 实现中，只写契约）

1. **链序**：`multi_radical`（binomial 后）先于 `inverse_trig`。
   契约：`multi_radical` 的基收集对**含反函数头的基一律视为不可归约并
   整体拒绝**（其片段是纯代数基）；M4 对**非伴随根式一律拒绝**（§2.2）。
   两者互不抢形状，无双 claim。
2. **交接形状**：00062/00873/01563 这类「双仿射基根式 × `(a+b·asin)^k`」
   的正确顺序是 **M4 先**——`u=asin(cx)` 后 `(d+cdx)^{1/2}(e−cex)^{3/2}`
   化为 `(1±sin u)` 的幂（三角多项式，已非根式），W4 根本不需要看到它们。
   若 M4 因其他原因拒绝，multi_radical 早已因反函数头拒绝——整体诚实回退。
3. **回积分的反向调用**：M4 的 t-形若残留纯根式项（不应发生——伴随根式
   已全部化为 `G′` 幂；出现即 bug），拒绝而非路由。

## 3. 接入点与预算常量

### 3.1 接入点

`inverse_trig.rs` 新增 `pub(crate) fn integrate_inv_subst(...)`（M4），
挂入 `integrate_inverse_trig` 分派序末尾（§2.6）；`mod.rs` 链不动。
**parts 纪律**：M4 在 `chain::parts_active()` 作用域内直接返回 `None`
（代换+回积分整体要么成功要么拒绝，与 ALGEBRAIC_RT 偏差记录第 1 条
同款抑制，避免 parts 内层残项升级烧链）。

### 3.2 预算常量（新增，全部确定性）

| 常量 | 值 | 约束对象 |
|---|---|---|
| `MAX_IS_NODES` | 200 | 输入被积函数节点数（沿用模块 `MAX_NODES`） |
| `MAX_IS_K` | 6 | 核幂 `(a+b·F)^k` 指数上界（对齐 `MAX_POWER_K`） |
| `MAX_IS_POLY_DEG` | 6 | 多项式权总次数（M2 为 4；三角族多角展开更贵，放宽一挡试） |
| `MAX_IS_RAD_M` | 5 | 伴随根式半幂 `|m|` 上界 |
| `MAX_IS_SUBST_NODES` | 400 | 代换后 t-形节点数 |
| `MAX_IS_REDUCED_NODES` | 600 | 多角归约后回积分节点数（沿用 `MAX_REDUCED_NODES`） |
| `MAX_IS_TRIG_EXP` | 8 | 三角族 sin/cos 幂积指数和（对齐 `MAX_HYP_EXP`） |
| `MAX_IS_WORK` | 100 000 | 阶段总工作量（对齐 `MAX_SQRT_WORK` 记账模型） |

### 3.3 输出形状

见 §2.4；均不自匹配 M1/P3/M4 入口模式（`F` 幂次不增、无新增根式基、
无 `F(G(·))` 复合），幂等不变式由回归测试钉死（§5 W2）。

## 4. 实值分支语义（单独成节）

### 4.1 主值区间（发射答案的成立区间）

| F | 主值 | 参数定义域 | 备注 |
|---|---|---|---|
| asin | `[−π/2, π/2]` | `\|w\| ≤ 1` | 区间内部严格成立，端点连续延拓 |
| acos | `[0, π]` | `\|w\| ≤ 1` | `G′ = −sin u ≤ 0`，符号进片层因子 |
| atan | `(−π/2, π/2)` | ℝ | 无片层问题（`sec² > 0`） |
| asinh | ℝ | ℝ | 双射，无分支 |
| acosh | `[0, ∞)` | `w ≥ 1` | `w ≤ −1` 侧**不做**（实值纪律） |
| atanh | ℝ | `\|w\| < 1` | `\|w\|>1` 侧是 acoth 领土，不做 |

符号系数下定义域不可判定时：照常发射（与 M1/M3 的既有约定一致——
形式导数精确是证书的语义范围，实值区间文档化由本节承担）。

### 4.2 片层因子

`E(w)^{m/2} = G′(u)^m·σ`，`σ = G′(u)·E(G(u))^{−1/2}`（`m` 奇时在场）：
逐区间 `±1` 常量、形式导数为零，`halfpower.rs` 模块文档 89–99 行的同一
机制。发射规则：有数值符号证据（能对目标区间判定 `G′(u)` 定号）折为
`±1`；否则保留 `σ` 节点并注释级记录；`σ` 参与证书（根式层关系
`r² = E(w)` 三角列判零，W3 已落地）。

### 4.3 证书看不到代换

证书认证的是最终 `D(F)−f ≡ 0`（各头以原始参数 `w` 出现）——代换只影响
求解路径，不影响认证强度（「数值证据 + 符号证书」双门分工的又一实例）。

## 5. 测试例（全部经 SymPy 1.14 核验，2026-11；核验方式逐例标注）

核验记号：**[D]** = `simplify(diff(F,x)−f)==0`；**[F]** = SymPy `integrate`
直接给闭式（对照组）；**[N]** = SymPy 回退（机制增益证据）；**[T]** =
代换 t-形经 SymPy 精确积分验证。

| # | 族 | 输入 | 期望/中间量 | 核验 |
|---|---|---|---|---|
| Z1 | W5a 零截距·三角 | `x^2*asin(x)` | `x³asin x/3 + x²√(1−x²)/9 + 2√(1−x²)/9` | [F][D] |
| Z2 | W5a 零截距·双曲 | `x^2*asinh(x)` | `x³asinh x/3 − x²√(x²+1)/9 + 2√(x²+1)/9` | [F][D] |
| Z3 | W5a 零截距·acos | `x*acos(x)` | `x²acos x/2 − x√(1−x²)/4 − acos x/4` | [F][D] |
| K1 | W5a 核幂·三角 | `(a+b*asin(x))^2` | `a²x + 2ab(x·asin x + √(1−x²)) + b²(x·asin²x − 2x + 2√(1−x²)asin x)` | [F][D] |
| K2 | W5a 核幂×多项式权 | `x^3*(1+2*asin(x))^2`（rubi-01205 数值形） | 闭式存在（asin²/asin/常数混合） | [F][D] |
| K3 | W5a 核幂·双曲 | `x*(a+b*asinh(c*x))^2` 族 | t-形 `(a+bu)²·sinh u·cosh u/c²` → 多角展开初等 | [T] |
| T3 | W5a bare-power | `x*atan(x)^2` | `x²atan²x/2 − x·atan x + log(x²+1)/2 + atan²x/2` | [F][D] |
| B2 | W5a bare-power·acos | `x^2*acos(x)^2`（rubi-00554 数值形） | 闭式存在 | [F][D] |
| H2 | W5a M1 回归 | `atanh(x)^2/(1-x^2)` | `atanh³x/3`（M1 先 claim，幂等回归） | [F][D] |
| A1 | W5a 仿射 | `(a+b*asinh(c+d*x))^2` | M2 已覆盖——**回归保护**（first claim 不变） | [F][D] |
| A4 | W5c 倒数头仿射 | `x^3*acsc(a+b*x)`（rubi-00511） | `u=acsc` ⇒ `u·(csc u−a)³(−csc u cot u)/b⁴`，一次分部后 `∫P(csc u)du` 初等；SymPy 回退 | [N][T] |
| CR1 | W5b 伴随根式 | `x*asin(x)/sqrt(1-x^2)` | `x − √(1−x²)·asin x` | [F][D] |
| CR2 | W5b 伴随根式×权 | `x^3*asin(x)/sqrt(1-x^2)`（rubi-01568 族） | t-形 `u·sin³u`；闭式存在；SymPy 给闭式 | [F][D] |
| CR3 | W5b atanh 伴随 | `x*(1-x^2)^(3/2)*atanh(x)`（rubi-01733 数值形） | t-形 `u·tanh u·sech u` → `−u·sech u + 2·atan(tanh(u/2))` 族；**SymPy 回退** | [N][T] |
| R1 | 拒绝·F-in-denom | `1/(1+asin(x))`（rubi-01301 族） | t-形 `cos u/(1+u)` → 非初等；**拒绝**；SymPy 回退 | [N] |
| R2 | 拒绝·负幂 bare | `1/acos(x)^3` | 非初等；**拒绝**；SymPy 回退 | [N] |
| R3 | 拒绝·非整数 k | `(1+asin(x))^(-5/2)`（rubi-01643 族） | 非初等（Fresnel 形）；**拒绝**；SymPy 回退 | [N] |
| R4 | 拒绝·polylog 涌现 | `atanh(1+x)/(2+2*x)`（rubi-00531） | t-形 `u·csch(2u)/…` → polylog；**拒绝**；SymPy 回退 | [N] |
| R5 | 拒绝·非伴随根式 | `x*(a+b*asinh(x))^2/(1+x^2)`（rubi-01225 族） | t-形 `(a+bu)²·tanh u` → polylog；**拒绝**；SymPy 回退 | [N] |
| R6 | 拒绝·非线性参数 | `a+b*asin(1+d*x^2)`（rubi-00152） | 既有拒绝护栏，M4 不得抢走 | 链回归 |
| W1 | first-claim 回归 | M1/M2/M3/P3 既有测试全绿 + §1.2 已解 8 例 | 答案形状逐字节不变 | 链回归 |
| W2 | 幂等 | Z1–CR3 全部答案再入链 | 无 `Integral` 回退、无循环（`MAX_CHAIN_ENTRIES` 内终止） | 链回归 |

每例三验（与既有口径一致）：符号证书（结构/初等域/根式层）+ 数值
oracle（`ocas-tests/src/integral_eval.rs`）+ SymPy 对拍（本表）。

## 6. 实现检查清单（0.30.0 Wave 5）

1. M4 骨架：`integrate_inv_subst` 签名、预算常量（§3.2）、模块文档 M4 段
   （代换表、初等性边界、拒绝面、分支语义指针）。
2. 代换构造器：12 头 `G`/`G′` 表（W5a 先接 6 主头，W5c 接倒数头）、
   仿射参数 `w=ρ+σx` 反解、`dx` 因子吸收、`(a+b·F)^k` 整数幂重表。
3. 伴随根式吸收（`radicand_scale` 复用）+ 片层因子 σ（§4.2）+ 非伴随
   根式拒绝。
4. t-形归约：三角族 sin/cos 幂积多角化（复用 `trig_reduce_products`）、
   tan/sec 偶幂归约、双曲族复用 M2 展开；`u^j·H(μu)` 表式收尾
   （M2 双曲表的三角类比）。
5. 回代 + `contains_integral` 终检 + 幂等形状断言；`parts_active` 抑制。
6. 接入 `integrate_inverse_trig` 分派序（M4 最末）；`mod.rs` 不动。
7. 与 W4 接口契约的断言测试：含反函数头的多基形状 multi_radical 拒绝、
   M4 非伴随根式拒绝（§2.7）。
8. 测试：§5 全表（Z/K/T3/B2/H2/A/CR/R/W 各组）+ 1892 探针回归 +
   历史错案护栏（`integral_certify.rs` 增补：把 `asin(sin u)` 复合残留
   或片层因子漏乘的伪答案必须被证书拒绝）。
9. 文档：mdBook 积分章补反函数代换一节（含 §4 分支语义）；CHANGELOG
   0.30.0 条目。
10. 验收口径：inverse-trig-hyper 桶 solved 8/146 → **≥48/146**
    （W5a 30–40 + W5b 15–25 的中值下沿），0 回归、证书门保持。

## 7. 明确不做

1. **逐例规则枚举**——Rubi 4.x–5.x 把本族展开为数百条按指数分案的规则；
   本引擎只实现代换/归约/回积分三段通用机械（§1.3-1）。
2. **F-in-denom 与负核幂族**（13+ 例；核验例 R1/R2/R3）——t-形为
   `P(u)^{−k}·T(u)`，非初等（Si/Ci/逆正切积分 Ti₂ 领土），诚实拒绝。
3. **polylog 需求族**（核验例 R4/R5：tanh/tan 奇幂带 u 权、`csch` 带
   u 权、分母不配对的核幂）——0.31.0 的 `Li₂` 落地后重评，本波不引入
   新特殊函数头。
4. **非线性参数** `F(c·x^n)`（12 例）——`x=G(u)^{1/n}` 重新引入根式，
   不闭合；`u=x²` 预代换子集（如 00444/01886 的 `/x^11`、`/x^4` 配对）
   列 W5c 可选评估，不进验收数。
5. **符号幂**（`x^m`、`(f·x)^m`、`F^n` 指数含符号，~10 例）——同 0.29.0
   符号指数纪律。
6. **acot 头**（题库本桶未出现）与**复值分支**（`w<−1` 的 acosh、
   `|w|>1` 的 atanh 等另一侧）——实值输出纪律（§4.1）。
7. **`RootSum`/隐式记号**——同 ALGEBRAIC_RT §7，输出一律显式。
8. **改动 Risch 塔**（§1.3-8）与 **`mod.rs` 链序**（§3.1）。

## 8. 风险与实测陷阱预判

1. **三角回代后的实值等价性是第一风险**。`cos u → √(1−c²x²)` 只在
   `cos u ≥ 0` 的片上成立；漏乘片层因子 σ 会在 `asin` 主值区间外差一个
   符号（`halfpower.rs` 的既有教训）。对策：§4.2 发射纪律 + 证书双门 +
   清单 8 的错案护栏。数值 oracle 的采样区间必须覆盖 σ 变号的两侧
   （`integral_eval.rs` 的域表沿用 `acosh` 既有先例）。
2. **证书对 `asin(sin u)` 类形状的判定**。若回代未完成就发射，
   `D(F)−f` 会含 `asin(sin(·))` 复合——初等域层把 `asin`/`sin` 列独立
   生成元后**无法**判零（`asin(sin u) = u` 是分支恒等式不是域恒等式），
   正确答案反被证书拒杀（假阴性）。对策：§2.4 的回代完备性终检
   （结构检查：答案中不得存在 `F(G(·))` 复合，先于证书运行）；
   P1 的消去带门是同款先例。
3. **与 W4 根式机制的联动顺序**。双仿射基×核幂形状（00062 族）若被
   multi_radical 先处理会把根式 Möbius 化但留下 `F(w)` 无法消去——
   §2.7 的契约（multi_radical 见反函数头即拒绝）是硬前提，W4 实现落地
   后第一件事是跑清单 7 的接口断言；若契约被违反，症状是 00062 族
   从「M4 可解」退为「双拒」，1892 回归会立刻显现。
4. **t-形回积分的链预算**。多角展开后项数 ×`u^j` 分部深度可能烧
   `MAX_CHAIN_ENTRIES`——`MAX_IS_WORK` 按 MULTI_RADICAL 风险 2 的记账
   模型把展开与回积分都计入，且 M4 在 parts 作用域内整体抑制（§3.1）。
5. **符号系数的 `σ ≢ 0` 与 `s` 定号判定**是解锁前提（题库几乎全带符号
   系数）：`collect_terms` 折叠后判零，判不了即拒绝；实测时首先验证
   （MULTI_RADICAL 风险 5 同款）。
6. **`linear_form` 复用必须过 0.30.0 已修复版本**（`x·x` 型重复变量
   因子返回 `None`；ALGEBRAIC_RT 偏差记录第 3 条）——参数仿射判定与
   多项式权解析都走它。
7. **M4 最后分派的代价**：M4 匹配器要对每个含反函数头的被积函数跑
   一遍代换构造尝试；匹配必须先做廉价的结构闸（头计数、参数仿射、
   `k` 整数、根式伴随判定）再进代换，保证拒绝路径 O(节点数) 有界——
   0.27.1 的 inverse_trig 超时归因（`declines_hang_attributed_shapes_
   quickly` 护栏）提醒：拒绝快与解得对同等重要。
