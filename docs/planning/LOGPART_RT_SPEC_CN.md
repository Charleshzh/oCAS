# Risch 各塔层对数部分一般结构定理（Rothstein–Trager 结式方法）算法规格

**版本**：oCAS 0.29.0（ℚ 根片段；代数数根推迟到 0.30.0）
**状态**：实现规格，不修改任何 Rust 代码
**核对来源**（本文所有公式均逐行核对于下列文件，行号以此为准）：

- SymPy 1.14.0：`ocas-tests/.venv/Lib/site-packages/sympy/integrals/risch.py`（下称 `risch.py`）、`ocas-tests/.venv/Lib/site-packages/sympy/integrals/rationaltools.py`（下称 `rationaltools.py`）。
- oCAS 现状：`ocas-calc/src/integral/risch.rs`、`ocas-calc/src/integral/rational.rs`、`ocas-calc/src/tower/elem.rs`、`ocas-calc/src/tower/build.rs`、`ocas-poly/src/resultant.rs`。

**对任务前提的一处更正**：`log_to_atan` / `log_to_real` **不在** `risch.py` 中；`risch.py:1280` 只有一行 TODO 注释 `# TODO: Use log_to_atan() from rationaltools.py`。两者实际位于 `rationaltools.py:279-325`（`log_to_atan`）与 `rationaltools.py:343-445`（`log_to_real`），见 §3。

---

## 1. 经典 Rothstein–Trager 算法

### 1.1 输入与前置条件

设第 `level` 层塔域为 `k(t)`，`k = ℚ(x, t₁, …, t_{level−1})`，`t` 为超越 `log` 或 `exp` 生成元，`D` 为塔导子。对数部分的输入是 Hermite 约化后的**简分数（simple）** `f = a1/d1`，满足：

| 条件 | 说明 | oCAS 保证来源 |
|---|---|---|
| `d1` 平方自由 | Hermite 余量分母 | `hermite_tower` 经 `square_free` 因子构造，`risch.rs:306-307`；`square_free` 见 `elem.rs:548-573` |
| `d1` 首一 | 便于验证恒等式 | `square_free` 内部 `monic()`（`elem.rs:553`、`511-517`）；首项系数为 1 还推出 `log` 层 `deg D d1 ≤ deg d1 − 1` |
| `gcd(a1, d1) = 1` | 简分数定义 | `KRat::new` 构造时约分（`elem.rs:592-606`）；实现时加一次 `d1.gcd(&a1).is_one()` 断言 |
| `deg a1 < deg d1` | 真分式 | 多项式部分已被 `div_rem` 分离（`risch.rs:193`） |
| `a1 ≠ 0`，`deg d1 ≥ 1` | 否则无对数部分 | 调用点判零（`risch.rs:248`） |

**输出**：`∫ a1/d1 dx` 的对数部分 `Σᵢ cᵢ·log(vᵢ)`，其中 `cᵢ` 为常量、`vᵢ ∈ k[t]` 首一；余部 `a1/d1 − D(Σ cᵢ log vᵢ)` 落入 `k`（log 层）或继续走多项式部分（exp 层恒等 `D log v` 为有理函数，余部恰好为 `k` 中元素，见例 1、3）。

### 1.2 结式的构造

```
q(z, t) = a1 − z·D d1        ∈ (k[z])[t]
R(z)    = resultant_t(q, d1) ∈ k[z]
```

逐行核对（SymPy `residue_reduce`，`risch.py:1258-1343`）：

- `Dd = derivation(d, DE)`、`q = a - pz*Dd`：`risch.py:1296-1297`；
- 结式：`risch.py:1299-1302`——`if Dd.degree(DE.t) <= d.degree(DE.t): r, R = d.resultant(q, includePRS=True) else: r, R = q.resultant(d, includePRS=True)`。注意 SymPy 的参数序是 `resultant(d, q)` 或按次数交换为 `resultant(q, d)`，与教科书 `Res(a − z·Dd, d)` 相差符号因子 `(−1)^{deg·deg}`；**根的集合不受整体符号影响**，oCAS 只需固定一种约定并记录（见 §4.3）；
- 计算前先约化：`a, d = a.cancel(d, include=True)`（`risch.py:1286`）、`d` 首一化（`risch.py:1287`）、`_, a = a.div(d)` 取余（`risch.py:1292`）——即 §1.1 的前置条件；
- `recognize_log_derivative` 中同一结式用于判定（`risch.py:1242-1245`）：`q = a - pz*Dd; r, _ = d.resultant(q, includePRS=True)`。

**次数界**：`deg_z R ≤ deg_t(a1 − z·D d1) ≤ max(deg a1, deg D d1) ≤ deg d1`（log 层 `D d1` 次数 ≤ `deg d1 − 1`，exp 层 `t ∤ d1` 时 `deg D d1 = deg d1`）。故 `deg d1 + 1` 个插值点必够（§4.4）。

### 1.3 根的求法与 `vᵢ` 的构造

定理（Bronstein *Symbolic Integration I* Thm 5.6.x；Rothstein–Trager）：`f` 在 `d1` 各根处的留数恰好是 `R(z)` 的全部根；对数部分初等 ⟺ **所有根都是微分域的常量**，且此时

```
∫ a1/d1 = Σᵢ cᵢ·log(vᵢ),   vᵢ = gcd(a1 − cᵢ·D d1, d1)   （对 R 的每个不同根 cᵢ 各取一次）
```

实现要点：

1. **根按不同值取，不计重数**。`R` 不必平方自由（例 3：`R = (1−2z)²`）；重根只贡献一个 `log` 项。SymPy 用 `splitfactor_sqf(r, DE, coefficientD=True, z=z)` 做平方自由分解（`risch.py:1309`），oCAS 在 ℚ 根片段里直接去重即可。
2. `cᵢ ∈ ℚ`（0.29.0 片段）时 `gcd(a1 − cᵢ·D d1, d1)` 是 **`k[t]` 上的普通 gcd**，用现成 `KPoly::gcd`（`elem.rs:496-508`），无需代数数运算。
3. `R(c)=0` ⟹ `deg gcd ≥ 1`（结式为零 ⟺ 有公共因子）；仍保留 `v.degree() > 0` 守卫，与基底层 `rational.rs:466-468` 一致。
4. SymPy 的经典版本判定函数 `recognize_log_derivative`（`risch.py:1228-1256`）：结式后对每个平方自由因子取 `real_roots` 并要求 `all(j.is_Integer ...)`（`risch.py:1249-1255`）——这是 SymPy 中「根必须是常量」的最严格实例（对数导数识别要求整数根）；一般 `residue_reduce` 的常量判常见 §2.2。

### 1.4 验证恒等式（自证闭环）

算完 `(cᵢ, vᵢ)` 后必须验证（与 0.28.0 起「不发未认证答案」的硬性不变量一致，`risch.rs:95-112`）：

```
(i)   Πᵢ vᵢ == d1                     （首一多项式逐系数 eq_cross）
(ii)  a1·Πᵢ vᵢ == d1·Σᵢ cᵢ·D vᵢ·Π_{j≠i} vⱼ     （k[t] 中交叉相乘等式，即 a1/d1 == Σ cᵢ·D vᵢ/vᵢ）
(iii) D cᵢ == 0                        （ℚ 根片段自动成立）
```

(i)(ii) 均在 `k[t]` 内精确判定：`KPoly` 系数 `KElem` 用 `eq_cross`（`elem.rs:177-182`）。任一失败 → 拒绝（返回 `None` 走 `integral_fallback_atom`，`risch.rs:744-757`），绝不可发错答案。顶层 `tower_difference_is_zero`（`risch.rs:119-132`）是最后一道网。

### 1.5 「根必须是常量」的判定

- 初等性判据：`R(z)` 的**每一个**根都是常量域 `Const_D(k)` 中的元素。存在非常量根 ⟹ `∫ a1/d1` 在 `k(t)` 上非初等（诚实拒绝，见例 5）；存在无理/虚常量根 ⟹ 初等但需要代数扩张（0.30.0，见例 4）。
- 0.29.0 的操作化判定：只接受 `cᵢ ∈ ℚ` 且 `R` 的有理根**完备**（抽出全部有理一次因子后余因子为常量）。完备性无法证明时一律拒绝（§2.3）。

---

## 2. ℚ 根片段（0.29.0）

### 2.1 问题

`R(z) ∈ k[z]`，`k = ℚ(x, t₁, …, t_{n−1})`：系数一般是含下层生成元的有理函数（例 2 的 `R = (1−z²)/x²`，系数 `±1/x² ∉ ℚ`）。需要从中**判定并提取有理常数根**，并保证完备性。

### 2.2 SymPy 怎么处理（逐行核对）

SymPy 的一般路径**从不求根**，而是惰性 `RootSum`：

- `residue_reduce` 返回 `(s_i, S_i)` 对，答案为 `g = Σ RootSum(s_i, λz: z·log(S_i(z,t)))`（docstring `risch.py:1271-1278`；组装见 `residue_reduce_to_basic`，`risch.py:1346-1355`：`sum(RootSum(a[0].as_poly(z), Lambda(i, i*log(a[1].as_expr()).subs({z: i})...)))`）。
- 常量性判定：`r = Poly(r, z)`（`risch.py:1308`）→ `splitfactor_sqf(r, DE, coefficientD=True, z=z)`（`risch.py:1309`；`splitfactor_sqf` 本体 `risch.py:994-1029`，其中 `kkinv = [z]` 见 `1008-1010`，`Si` 为 `pi` 与其系数导子的 gcd 见 `1018-1020`）→
  `b = not any(cancel(i.as_expr()).has(DE.t, z) for i, _ in Np)`（`risch.py:1341`）。`coefficientD=True` 的导子只导 `k[z]` 的系数：系数全为常量的因子整除自身导数（导数为 0），全部落入 special 侧 `S`；凡系数含非常量的因子必在 normal 侧 `N` 留下真因子，使 `b = False`。即 **`b` ⟺ `R` 的每个平方自由因子的系数全在常量域**（⟺ 根全为常量）。
- 显式求根只出现在判定函数 `recognize_log_derivative`（`risch.py:1249-1255`：`real_roots` + `is_Integer`）与有理函数前端 `rationaltools.py`（`_get_real_roots`，`rationaltools.py:328-340`；`log_to_real`，`rationaltools.py:388-443`）。

oCAS 0.29.0 没有 `RootSum` 记号，因此必须显式求根——但只在 ℚ 内。

### 2.3 oCAS 候选提取算法（规格）

输入 `R(z) = Σᵢ rᵢ·zⁱ`，`rᵢ ∈ KElem`（`KElem = Sparse/Sparse`，系数域 ℚ，`elem.rs:93-100`）。

**步骤 1（快路径）**：若所有 `rᵢ.as_rational()`（`elem.rs:184-187`）成功，则 `R ∈ ℚ[z]`：转为 `DPoly` 复用基底层现有 `rational_roots`（`rational.rs:512-554`；整数化 → `primitive_part().factor()` → 一次因子给根、`fully_split` 标志）。要求 `fully_split`，否则拒绝。

**步骤 2（一般 k）**：有理根定理在 `ℤ[x, t₁, …][z]` 上的形式。

1. 把各 `rᵢ = nᵢ/dᵢ` 通分：乘 `M = Πⱼ dⱼ`（无需多元 lcm——`KElem` 本就不做多元 gcd，`elem.rs:9-12`），再乘所有系数的 ℚ 分母之 lcm，得 `N(z) = Σᵢ aᵢ·zⁱ ∈ ℤ[x, t₁, …][z]`（`aᵢ` 为整系数 `Sparse`）。
2. 设 `a₀`（常数项）、`aₘ`（首项）的**整数含量** `g₀ = gcd(a₀ 的所有整数系数)`、`gₘ = gcd(aₘ 的所有整数系数)`。任何最简有理根 `c = p/q`（`p,q ∈ ℤ`）满足 `p | a₀`、`q | aₘ`（于 `ℤ[x,…]` 中），故 `p | g₀`、`q | gₘ`。**候选集有限**：`{±p/q : p | g₀, q | gₘ}`。
3. 对每个候选 `c`，精确代入检验 `N(c) == 0`：`KElem` 求值（`Sparse` 代入常数）后看分子是否为零多项式（零检测可靠，`elem.rs:9-12`）。
4. **完备性**：每确认一个有理根 `c`，用综合除法把 `N(z)`（视为系数在 `k` 的 `z` 多项式）除以 `(z − c)`。全部有理根抽出后，**余因子必须是不含 `z` 的常量**；否则拒绝——余因子的根或为无理/虚常量（0.30.0），或为非常量（非初等），0.29.0 无法区分也不必区分，统一诚实拒绝（对应基底层 `fully_split == false` → 未求值 `Integral` 回退，`rational.rs:458-461`）。

**拒绝策略**（诚实、可回退）：步骤 1/2 任何一步失败 → 该对数部分不产生 `log` 项，调用点退回 `integral_fallback_atom`（`risch.rs:252-255`），与现状行为一致；**绝不**输出部分猜测的封闭形式。

---

## 3. Lazard–Rioboo 改进

**区别（一句话）**：经典 Rothstein–Trager 对每个根 `cᵢ` 在扩域 `k(cᵢ)[t]` 上各做一次 `gcd`（需要代数数运算）；Lazard–Rioboo(-Trager) 用**同一条子结式 PRS** 中次数与根次数匹配的余式 `S_i(z, t)` 代替逐根 `gcd`，全部计算留在 `k` 上，结果以 `RootSum` 打包。SymPy 两侧都是这个变体：`residue_reduce` 自述 *"Lazard-Rioboo-Rothstein-Trager resultant reduction"*（`risch.py:1260`，PRS 取自 `resultant(..., includePRS=True)`，`risch.py:1300/1302`，按 `R_map[i.degree()]` 取余式，`risch.py:1304-1306`）；有理函数前端 `ratint_logpart` 自述 *"Lazard-Rioboo-Trager algorithm"*（`rationaltools.py:190`，同样 `includePRS=True`，`rationaltools.py:231`）。因此 **SymPy 用的是 Lazard-Rioboo-Trager，不是经典逐根 gcd 的 Rothstein-Trager**。

**0.29.0 建议：不需要**。ℚ 根片段下 `vᵢ = gcd(a1 − cᵢ·D d1, d1)` 本就在 `k[t]` 内（`KPoly::gcd` 现成），逐根 gcd 成本可忽略；LRT 的收益（避免代数数 gcd）要等 0.30.0 引入代数数根时才兑现，届时还需配套 `RootSum` 记号或代数数域，再评估不迟。

附：`log_to_atan`（`rationaltools.py:279-325`，`2*atan(u)` 递归消去复对数 `I·log((f+Ig)/(f−Ig))`）与 `log_to_real`（`rationaltools.py:343-445`，`t → u+Iv` 分裂、共轭根配对、实根过滤）服务于**实域输出美化**（复对数 → 实 `log`/`atan`），属于基底层有理积分的事后处理；oCAS 基底层已对二次分母直接发 `atan`（`rational.rs` Case 2），0.29.0 塔层规格不涉及。

---

## 4. 到 oCAS 类型的映射

### 4.1 为什么 `ocas-poly` 的现成 resultant 不能直接用

`ocas-poly/src/resultant.rs:40-162` 是 `DenseUnivariatePolynomial<D: EuclideanDomain>` 上的 Brown 子结式 PRS：依赖整环上的**伪余式 + beta 精确除法**（`resultant.rs:104-126`，两处 `.expect("subresultant ... division is exact")` 见 `92`、`123`）。`KPoly` 的系数是 `KElem`（多元有理函数，非 `EuclideanDomain` 类型），无法实例化该泛型；且域上根本不需要伪除法。

### 4.2 推荐：`KPoly::resultant`（域上 Euclidean 递归）

系数域 `k` 是**域**，普通带余除法即精确（`KPoly::div_rem`，`elem.rs:469-493`，首项系数相除走 `KElem::div` 的域逆元），无零除问题：

```
resultant(a, b):  // a, b ∈ k[t]，返回 k
  if a.is_zero() || b.is_zero() → 0                      // 守卫，正常路径不可达
  if deg b == 0 → b₀^(deg a)                             // KElem::pow（elem.rs:236）
  (q, r) = a.div_rem(b)                                  // 精确除法（elem.rs:469）
  if r.is_zero() → (deg b ≥ 1 ? 0 : lc(b)^(deg a))
  else → (−1)^(deg a·deg b) · lc(b)^(deg a − deg r) · resultant(b, r)
```

要点：

- **零除**：`KElem::div` 返回 `Option`（`elem.rs:227-229`），仅当除数为零失败；`KPoly` 的 trim 不变量保证 `lc()` 非零（`elem.rs:275-283`、`358-363`）。`KElem` 的零检测看分子，可靠（`elem.rs:9-12`）。
- **分数膨胀**：`KElem` 不做多元约分（`elem.rs:9-12`），PRS 中间分数会膨胀。0.29.0 规模下接受；若成为瓶颈，每步把余式首一化并补记账：`r = lc(r)·monic(r)` ⟹ `Res(b, r) = lc(r)^(deg b − deg r)·Res(b, monic(r))`。
- **性能备选**：欲消膨胀可移植 `resultant.rs` 的子结式 PRS——域上 beta 除法同样精确（子结式定理），但 0.29.0 不建议先做。

### 4.3 符号约定

固定 `R(z) = resultant_t(a1 − z·D d1, d1)`（第一参数为 `q`，与 SymPy `recognize_log_derivative` 的 `d.resultant(q)` 序不同，差 `(−1)^{deg d·deg q}`，根集合不变）。`KPoly::resultant` 内部交换约定：若 `deg a < deg b`，`Res(a,b) = (−1)^{deg a·deg b}·Res(b,a)`——与 `ocas-poly` 现有实现一致（`resultant.rs:43-54`）。常量情形：`Res(c₀, d₀) = 1`（0×0 Sylvester 行列式约定，现有测试 `resultant.rs:221-227`）。**根的提取对 `R` 的整体符号与含量免疫**（§2.3 步骤 2 通分时已乘任意非零因子），但 `vᵢ` 一律取首一。

### 4.4 `R(z)` 的具体构造：插值路径（推荐，对照基底层）

不引入 `z` 作为新多项式变量，改为 `KElem` 值插值——完全镜像基底层 `rothstein_trager`（`rational.rs:424-475`）：

```
n = deg d1
for j in 0..=n:                               // n+1 个点（次数界见 §1.2）
    u_j = a1.sub(&dd1.mul_kelem(&KElem::from_rational(j)))   // a1 − j·D d1
    y_j = u_j.resultant(&d1)                  // KElem
R(z) = lagrange_interpolate(z 视为新变元, [(j, y_j)])   // 节点为有理数，分母 ℚ 精确
```

- 节点 `z_j = j` 处 `deg_t u_j` 可能掉次（如 `z=0` 且 `deg D d1 > deg a1`）：~~不影响正确性——结式是系数的多项式，特化与求值可交换~~ **更正（0.29.0 实现证伪）**：形式次模板下结果式与特化在掉次节点**不可交换**——`a1 = 1, d1 = 1 + t, D d1 = t` 在节点 0 的朴素结果式是 `Res(1, 1+t) = 1`，而形式多项式 `R(z) = −(1+z)` 要求 `R(0) = −1`（Sylvester 矩阵的模板次数随掉次而变）。实现按 `max(deg a1, deg D d1)` 的形式次模板求值并**跳过掉次节点**（至多一个：首项系数之比），其余节点插值恢复；`logpart.rs` 的 `rt_resultant_polynomial` 与单元测试 `rt_exp_level_single_root` 是该更正的落点。
- Lagrange 分母只含有理节点差，精确（对照 `rational.rs:482-510`）。插值结果是「系数为 `KElem` 的 `z` 多项式」，建议表示为 `Vec<KElem>`（稠密、升幂）+ 辅助例程（求值、综合除法、乘以 `(z−c)`）。
- 备选路径（0.30.0 再议）：定义 `KZPoly`（系数 `KElem` 的 `z` 多项式）并把 §4.2 递归直接跑在 `(k[z])[t]` 上，一次算出符号 `R(z)`；0.29.0 代码量不划算。

### 4.5 接入点与输出形状

替换/扩展 `integrate_level` 的对数部分块（`risch.rs:247-256`）：

```
if !a1.is_zero() {
    let dd1 = tower_diff_kpoly(&d1, &tower.gens[..level - 1], &tgen.dt);   // build.rs:87-94
    if let Some(c) = kpoly_scalar_multiple(&a1, &dd1) {                     // 恒等式快路径保留
        out.logs.push((c, d1.kelem()));
    } else {
        match rothstein_trager_tower(&a1, &d1, &dd1) {   // 本规格 §1–§4，含 §1.4 验证
            Some(logs) if 完整 → out.logs.extend(logs),
            _ → out.extras.push(integral_fallback_atom(...)),              // 现状回退
        }
    }
}
```

输出形状零改动：`LevelResult.logs: Vec<(Rational, KElem)>`（`risch.rs:147`）本来就允许任意有理 `cᵢ` 与 `KElem` 的 `vᵢ`（`KPoly::kelem()` 收敛为场元素，`elem.rs:339-348`）。

---

## 5. 测试例

记号：塔 `ℚ(x, t)`，`D x = 1`。每个例子给出 `a1, d1, D d1, R(z), cᵢ, vᵢ` 与最终答案；`R` 按 §4.3 约定 `resultant_t(a1 − z·D d1, d1)` 计算（二元 Sylvester：`Res(a₁t+a₀, b₁t+b₀) = a₁b₀ − a₀b₁`；`Res(u, Π(t−βⱼ)) = lc^…·Π u(βⱼ)`，差号不影晌根）。

### 例 1（exp 层，单有理根）：`∫ dx/(1 + exp(x)) = x − log(1 + exp(x))`

`t = exp(x)`，`D t = t`；`f = 1/(1+t)` 已为简分式（Hermite 平凡）。

- `a1 = 1`，`d1 = 1 + t`（首一），`D d1 = t`。
- `q = 1 − z·t`；`R(z) = Res_t(1 − zt, 1 + t) = (−z)·1 − 1·1 = −(1 + z)`。
- 根 `c₁ = −1 ∈ ℚ` ✓；完备（一次式）。
- `v₁ = gcd(a1 − c₁·D d1, d1) = gcd(1 + t, 1 + t) = 1 + t`。
- 验证：`D(−log(1+t)) = −t/(1+t) = 1/(1+t) − 1` ⟹ 余部 `1 ∈ k`，`∫1 dx = x`。
- **答案**：`x − log(1 + exp(x))`。
- 恒等式快路径：`a1 == c·D d1` 即 `1 == c·t`，次数不符 → **失败，必须走结式**。

### 例 2（log 层，两个有理根；R 系数含下层生成元）：`∫ dx/(x·log(x)·(log(x)+1)) = log(log x) − log(log x + 1)`

`t = log(x)`，`D t = 1/x`；`f = (1/x)/(t² + t)`，`d1 = t² + t = t(t+1)` 平方自由、`gcd(d1, D d1) = 1`。

- `a1 = 1/x`，`d1 = t² + t`，`D d1 = (2t + 1)/x`。
- `q = (1 − z)/x − (2z/x)·t`；`d1` 的根 `t = 0, −1`：
  `R(z) = q(0)·q(−1) = ((1−z)/x)·((1−z+2z)/x) = (1 − z²)/x²`。
- **§2.3 步骤 2**：通分乘 `x²` 得 `N(z) = 1 − z² ∈ ℤ[x][z]`，`a₀ = 1, aₘ = −1` ⟹ 候选 `±1`；代入确认两根；除尽后余常量 ⟹ 完备 ✓。（步骤 1 快路径此处不适用：系数 `±1/x² ∉ ℚ`。）
- `c₁ = 1`：`q(1) = −2t/x`，`v₁ = gcd(−2t/x, t²+t) = t`；`c₂ = −1`：`q(−1) = (2t+2)/x`，`v₂ = t + 1`。
- 验证 (i) `t·(t+1) = d1` ✓；(ii) `D(log t − log(t+1)) = (1/x)(1/t − 1/(t+1)) = (1/x)/(t²+t)` ✓。
- **答案**：`log(log x) − log(log x + 1)`。
- 恒等式路径：`1/x` 不是 `(2t+1)/x` 的常数倍 → **必须走结式**。

### 例 3（R 含二重 z 根，根按不同值取一次）：`∫ dx/(exp(2x) − 1) = ½·log(exp(2x) − 1) − x`

`t = exp(x)`，`f = 1/(t² − 1)`；`d1 = t² − 1` 平方自由，`D d1 = 2t²`，`gcd(t²−1, 2t²) = 1`。

- `a1 = 1`，`q = 1 − 2z·t²`；`d1` 的根 `t = ±1`：`R(z) = q(1)·q(−1) = (1 − 2z)²`。
- **二重根 `c₁ = 1/2` 只取一次**（`R` 不必平方自由；§1.3 要点 1）。§2.3：`N = (1−2z)²`，候选 `p|1, q|4`：`±1, ±1/2, ±1/4`；确认 `1/2`；两次综合除法后余常量 ⟹ 完备。
- `v₁ = gcd(1 − t², t² − 1) = t² − 1`（首一）。
- 验证：`½·D(t²−1)/(t²−1) = t²/(t²−1) = 1 + 1/(t²−1)` ⟹ 余部 `−1`，`∫(−1)dx = −x`。
- **答案**：`½·log(exp(2x) − 1) − x`。

### 例 4（无理数根 → 0.29.0 拒绝）：`∫ dx/(x² − 2)`

基底层 `ℚ(x)`（同一规格，`rational.rs` 已有此行为）。

- `a1 = 1`，`d1 = x² − 2`，`d1' = 2x`；`R(z) = Res(1 − 2zx, x² − 2) = (1 − 2√2·z)(1 + 2√2·z) = 1 − 8z²`。
- 根 `c = ±1/(2√2) = ±√2/4 ∉ ℚ`：`rational_roots` 的 `fully_split = false`（`rational.rs:544-552`）→ 回退未求值 `Integral`（`rational.rs:458-461`）。
- **期望行为**：0.29.0 诚实拒绝（真答案 `(√2/4)·log((x−√2)/(x+√2))` 留给 0.30.0 代数数根）。

### 例 5（非常量根 → 拒绝，非初等）：`∫ dx/(exp(x) + x)`

`t = exp(x)`；`d1 = t + x` 平方自由，`D d1 = t + 1`，`gcd(t+x, t+1) = 1`；`a1 = 1`。

- `q = 1 − z(t + 1)`；`R(z) = Res_t(q, t + x) = (−z)·x − (1 − z)·1 = z(1 − x) − 1`。
- 根 `c = 1/(1 − x) ∈ ℚ(x) ∖ ℚ`：**§2.3 步骤 2** 的候选集：`N = z(1−x) − 1 ∈ ℤ[x][z]`，`a₀ = −1, aₘ = 1 − x`，`g₀ = gₘ = 1` ⟹ 候选 `±1`；代入 `R(1) = −x ≠ 0`、`R(−1) = x − 2 ≠ 0` ⟹ **无有理根**；抽出零个因子后余因子仍含 `z` ⟹ 拒绝 ✓。
- 数学上一致：`c` 对 `x` 非常量（`D c = 1/(1−x)² ≠ 0`）⟹ 积分非初等。
- **期望行为**：返回 `None` / 未求值 `Integral`。

### 例 6（含符号系数的行为）：`∫ dx/(1 − exp(x)) = x − log(exp(x) − 1)`

`t = exp(x)`；首一化 `d1 = t − 1`（`a1 = −1`，符号随首一化翻转），`D d1 = t`。

- `q = −1 − z·t`；`R(z) = Res_t(−1 − zt, t − 1) = (−z)·(−1) − (−1)·1 = z + 1`。
- 根 `c₁ = −1`；`v₁ = gcd(−1 + t, t − 1) = t − 1`。
- 验证：`−D(t−1)/(t−1) = −t/(t−1) = 1/(1−t) − 1 = f − 1` ⟹ 余部 `1` → `x`。
- **答案**：`x − log(exp(x) − 1)`。
- **行为说明**：(a) `d1` 首一化会把符号搬进 `a1`，`R` 可能整体变号（`Res(−u, v) = (−1)^{deg v}Res(u, v)`），根集合不变；(b) `cᵢ` 可为负或分数（例 1 `−1`、例 3 `1/2`），`LevelResult.logs` 的 `(Rational, KElem)` 已覆盖；(c) `vᵢ` 一律首一输出，保证 `log` 参数形状规范（对照 SymPy `monic()`，`risch.py:1313/1327`）。

### 例 7（恒等式已够的对照组）：`∫ dx/(x·log x) = log(log x)`

`t = log(x)`；`a1 = 1/x`，`d1 = t`，`D d1 = 1/x`。`a1 == 1·D d1` ⟹ 快路径直接给 `(1, t)`。若走结式：`R(z) = Res_t((1−z)/x, t) = (1−z)/x`，根 `1`，`v = gcd(0, t) = t = d1`——**结式包含恒等式**（见 §6）。

---

## 6. 与现有恒等式快速路径的关系

现有快路径（`risch.rs:247-256`：`kpoly_scalar_multiple(&a1, &dd1)`，`risch.rs:393-413`）处理的恰是 **`R(z)` 只有一个不同根 `c` 且 `v = d1`** 的特例：

```
a1 == c·D d1  ⟺  a1/d1 == c·D d1/d1  ⟹  R(z) = (c − z)^{deg d1}·Res(D d1, d1)，唯一根 c，v = gcd(0, d1) = d1
```

- **恒等式足够**：单根、全分母单一 `log`——例 7，以及底层 `∫ D u/u` 型。
- **恒等式不够、必须结式**：(a) `a1` 与 `D d1` 次数不符但结式有有理根（例 1、6：`1` vs `t`）；(b) 多个不同留数（例 2：`±1`）；(c) 留数为分数（例 3：`1/2`）；(d) 判定非初等/越界片段并**诚实拒绝**（例 4、5）——这一点快路径做不到区分，只能笼统回退，结式路径给出精确归因。
- **保留策略**：恒等式是 O(次数) 的比例检测，结式是 `deg d1 + 1` 次 `KPoly::resultant`；保留恒等式作为前置快路径，失败再进 §1–§4 流程，两路径输出统一经 §1.4 验证与顶层 `tower_difference_is_zero`（`risch.rs:108-111`）认证。

---

## 7. 实现检查清单（0.29.0）

1. `KPoly::resultant`（§4.2）+ 单元测试（对照 `ocas-poly` 测试值与本文例 1–3 的 `R`）。
2. `R(z)` 插值构造（§4.4，系数 `KElem` 的稠密 `z` 多项式）。
3. ℚ 根提取：步骤 1 复用 `rational_roots`；步骤 2 候选 + 精确代入 + 综合除法完备性（§2.3）。
4. 逐根 `vᵢ = gcd(a1 − cᵢ·D d1, d1)` + §1.4 三项验证，失败即回退。
5. 接入 `integrate_level`（§4.5），保留恒等式快路径。
6. 测试：例 1–6 全过（例 4、5 断言拒绝/回退），例 7 断言仍走通；`cargo fmt/clippy/test` 按 CLAUDE.md 清单。
7. 明确不做：Lazard–Rioboo（§3）、代数数根（0.30.0）、`RootSum` 记号。
