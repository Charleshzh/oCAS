# 多项式因式分解

oCAS 实现了整数和素有限域上一元与二元多项式的因式分解。本章介绍相关算法、适用范围以及公共 API。

---

## 适用范围

当前支持的因式分解包括：

| 域 | 多项式 | 算法 |
|---|---|---|
| $\mathbb{Z}[x]$ | 一元 | 无平方分解、Berlekamp–Zassenhaus、Hensel 提升 |
| $\mathbb{F}_p[x]$ | 一元 | 无平方分解、Berlekamp |
| $\mathbb{Z}[x,y]$ | 二元（关于 $x$ 首一） | Wang Hensel 提升 |
| $\mathbb{F}_p[x,y]$ | 二元（关于 $x$ 首一） | 有限域上的 Hensel 提升 |
| $\mathbb{Z}[x_1,\dots,x_n]$ | 任意多元 | Wang EEZ Hensel 提升 + 首项系数预处理 + Zassenhaus 重组 |
| $\mathbb{F}_p[x_1,\dots,x_n]$ | 任意多元 | EEZ Hensel 提升（含特征 $p$ 的 $p$ 次幂处理；非常数 LC 走域版 Wang 预处理） |
| $\mathbb{Q}(\alpha)[x]$ | 一元 | Trager 算法：范数 + $\mathbb{Q}$ 上分解 + $\mathrm{GF}(p^d)$ 模 GCD |

自 0.16.0 起支持多于两个变量的任意多元因式分解。非首一的一元多项式经首项系数变换后同样可分解。自 0.16.1 起，$\mathbb{Z}$ 上的多元路径通过 p-adic 系数 Hensel 提升（骨架插值 Diophantine）完整支持非常数首项系数强加；二元非常数 LC 输入自动转入 EEZ 路径。自 0.17.0 起支持代数数域上的一元因式分解（Trager 算法）。

---

## 有限域上的一元因式分解

对于素域 $\mathbb{F}_p$，oCAS 使用 **Berlekamp 算法**。首先将多项式化为无平方形式，然后计算 Frobenius 矩阵 $Q - I$ 的核空间。核空间的每个基向量通过与非平凡元素的 gcd 给出因子分解。

```rust
use ocas_domain::{FiniteField, Integer};
use ocas_poly::DenseUnivariatePolynomial;

let p = FiniteField::new(Integer::from(101));
let f = DenseUnivariatePolynomial::<FiniteField>::from_coeffs(p.clone(), vec![
    p.element(1), // 常数项
    p.element(0),
    p.element(1), // x^2
]);
let factors = f.factor();
```

返回结果是 `(因子, 重数)` 的列表。经过无平方分解后，有限域上的重数总是 1。

---

## 整数上的一元因式分解

对于 $\mathbb{Z}[x]$，oCAS 将**无平方分解**与 **Berlekamp–Zassenhaus 型 Hensel 提升**结合使用：

1. 计算内容并约化为本原多项式。
2. 计算无平方分解。
3. 选取一个小素数 $p$，使得模 $p$ 约化后仍无平方且次数与原多项式相同。
4. 使用 Berlekamp 在模 $p$ 下分解。
5. 通过 Hensel 提升将模 $p$ 因子提升为 $\mathbb{Z}[x]$ 上的因子。

```rust
use ocas_domain::IntegerDomain;
use ocas_poly::DenseUnivariatePolynomial;

let f = DenseUnivariatePolynomial::<IntegerDomain>::from_coeffs(
    IntegerDomain,
    vec![1.into(), 0.into(), 1.into()], // x^2 + 1
);
let factors = f.factor();
```

---

## 整数上的二元因式分解

`ocas-poly` 使用 **Wang Hensel 提升** 对 $\mathbb{Z}[x,y]$ 上的二元多项式进行因式分解，假设多项式关于 $x$ 首一。

算法流程：

1. 选取赋值点 $y = \alpha$，使得一元像 $f(x, \alpha)$ 无平方且不可约因子数最少。
2. 在 $\mathbb{Z}[x]$ 上分解该一元像。
3. 通过修正 $f$ 在 $y = \alpha$ 处的 Taylor 系数，将一元因子提升回二元因子。

修正步骤使用一元因子的有理 Bézout 系数，然后通过整数除法重建整系数修正项。如果重建失败（余式非整）或提升后的乘积与原多项式不符，实现会尝试下一个候选赋值点，最终回退为将该多项式视为不可约。

```rust
use ocas_domain::IntegerDomain;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::Lex;

type MPoly = SparseMultivariatePolynomial<IntegerDomain, Lex>;

let domain = IntegerDomain;
let f = MPoly::from_terms(domain, 2, vec![
    (vec![3, 0], 1.into()),  // x^3
    (vec![2, 1], 1.into()),  // x^2*y
    (vec![2, 0], 2.into()),  // 2*x^2
    (vec![1, 1], 1.into()),  // x*y
    (vec![1, 0], 1.into()),  // x
    (vec![0, 2], 1.into()),  // y^2
    (vec![0, 1], 3.into()),  // 3*y
    (vec![0, 0], 2.into()),  // 2
]);

let factors = f.factor();
// factors 包含 (x^2 + y + 1, 1) 和 (x + y + 2, 1)
```

---

## 有限域上的二元因式分解

在 $\mathbb{F}_p[x,y]$ 上，使用相同的 Hensel 提升框架，但所有运算直接在有限域中进行。Bézout 系数通过有限域上的 gcd 计算，所有修正项都保持在域中，因此无需整系数重建。

```rust
use ocas_domain::FiniteField;
use ocas_poly::SparseMultivariatePolynomial;
use ocas_poly::sparse::Lex;

let domain = FiniteField::new(101.into());
type FpPoly = SparseMultivariatePolynomial<FiniteField, Lex>;

let f = FpPoly::from_terms(domain.clone(), 2, vec![
    (vec![2, 0], 1.into()), // x^2
    (vec![0, 1], 1.into()), // y
    (vec![0, 0], 1.into()), // 1
]);
let factors = f.factor();
```

---

## C/C++ 多项式 API

C 绑定为二元整数和有限域多项式提供不透明句柄。多项式可通过 ASCII 字符串创建、因式分解、打印和释放。

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;
    OcasPolyZ *f = ocas_poly_z_create("x^2 + y + 1", &err);
    if (f == NULL) {
        fprintf(stderr, "parse error: %s\n", ocas_error_last_message());
        return 1;
    }

    OcasPolyFactorArray factors = {0};
    int rc = ocas_poly_z_factor(f, &factors, &err);
    if (rc != OCAS_OK) {
        fprintf(stderr, "factor error: %s\n", ocas_error_last_message());
        ocas_poly_z_free(f);
        return 1;
    }

    printf("factors: %zu\n", factors.len);
    for (size_t i = 0; i < factors.len; ++i) {
        OcasPolyZ *factor = (OcasPolyZ *)factors.factors[i].poly;
        char *s = ocas_poly_z_to_string(factor, &err);
        printf("  %s^%zu\n", s, factors.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_z_free(factor);
    }
    ocas_poly_factor_array_free(&factors);
    ocas_poly_z_free(f);
    return 0;
}
```

有限域版本使用 `OcasPolyFp`、`ocas_poly_fp_create` 和 `ocas_poly_fp_factor`。因子数组中存储的是 `void*` 多项式句柄，调用者必须先将其转换为正确的具体类型，再进行打印或释放。

---

## 任意多元因式分解（Wang EEZ）

自 0.16.0 起，oCAS 支持 $\mathbb{Z}$ 与 $\mathbb{F}_p$ 上任意变量数的多项式因式分解，核心是 **Wang 的 EEZ（Evaluation and EZ-lifting）算法**：

1. **无平方分解**（Yun）：对主变量 $x_1$ 求偏导，用 $n$ 元 GCD（稠密递归求值–插值）剥离重因子；在特征 $p$ 下处理 $p$ 次幂（收缩/展开）。
2. **求值点采样**：把副变量 $x_2,\dots,x_n$ 代入样本点，使一元像次数不变且无平方。
3. **Wang 首项系数预处理**：分解首项系数 $\ell(x_2,\dots)$，用其在样本点的两两互素整数像把各因子分配到各一元因子，重建真实首项系数 $\ell_i$。
4. **逐变量 Hensel 提升**：通过多元 Diophantine 方程把一元因子沿理想 $(x_k - a_k)$ 逐变量提升回多元。
5. **Zassenhaus 重组**：当提升的因子多于真实因子时，枚举子集组合并试除，得到不可约因子。

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::sparse::{Lex, SparseMultivariatePolynomial};

// (x + y + z)(x - y + 2z)
let f = SparseMultivariatePolynomial::<IntegerDomain, Lex>::from_terms(
    IntegerDomain, 3,
    vec![
        (vec![2, 0, 0], Integer::from(1)),
        (vec![1, 0, 1], Integer::from(3)),
        (vec![0, 1, 1], Integer::from(1)),
        (vec![0, 2, 0], Integer::from(-1)),
        (vec![0, 1, 0], Integer::from(0)),
        (vec![0, 0, 2], Integer::from(2)),
    ],
);
let factors = f.factor();
```

返回 `(因子, 重数)` 列表，因子按本原且首项系数为正规范化。

---

## 代数数域上的因式分解（Trager）

自 0.17.0 起，oCAS 支持代数数域 $K = \mathbb{Q}(\alpha)$ 上一元多项式的因式分解，其中 $\alpha$ 由 $\mathbb{Q}$ 上的首一极小多项式给出。域类型位于 `ocas-domain`：`AlgebraicExtension<D>`（$D = \mathbb{Q}$ 时为 `AlgebraicNumberField`）；同一类型在 $D = \mathbb{F}_p$ 时表示 $\mathrm{GF}(p^d)$。

分解器实现 **Trager 算法**：

1. **无平方分解**（Yun）：在 $K$ 上进行，GCD 使用下述模方法而非伪余式序列。
2. **范数与平移**：对 $s = 0, 1, 2, \dots$ 计算 $g(x) = f(x - s\alpha)$ 的范数——通过求值–插值计算标量结式 $\operatorname{Res}_\alpha(m, g)$——直到范数在 $\mathbb{Q}$ 上无平方（用小素数模检验；接受条件是精确的）。
3. **有理数域分解**：无平方范数在 $\mathbb{Z}$ 上分解（Hensel 路径）。
4. **数域上的模 GCD**：把每个有理因子映到 $K[x]$，用模方法求其与 $g$ 的 GCD——在 $m \bmod p$ 不可约的素数上映到 $\mathrm{GF}(p^d)$，CRT 合并首一模 GCD，有理重构系数，试除验证。
5. **回代** $x \mapsto x + s\alpha$ 并首一化。

当 $f$ 的所有系数都是有理常数时走**有理快速通道**：先在 $\mathbb{Q}$ 上分解，仅把 $\mathbb{Q}$ 不可约因子送入 Trager 范数流程。

```rust
use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;

// ℚ(√2)：极小多项式 α² − 2。
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
// x² − 2 = (x − √2)(x + √2) 在 ℚ(√2) 上分裂。
let f = DenseUnivariatePolynomial::from_coeffs(
    field.clone(),
    vec![field.from_base(Rational::new(-2, 1)), field.zero(), field.one()],
);
let factors = f.factor();
assert_eq!(factors.len(), 2);
```

基准（criterion 组 `poly_factor_anf`）：$\mathbb{Q}(\sqrt2)$、$\mathbb{Q}(\sqrt[3]{2})$ 与分圆域 $\mathbb{Q}(\zeta_5)$ 上次数 ≤ 12 的多项式分解均远低于 100 ms。

---

## 限制与未来工作

- $\mathbb{Z}$ 上的多元路径自 0.16.1 起已完整支持非常数首项系数强加（p-adic 系数 Hensel 提升 + 骨架插值 Diophantine）。
- 二元 Wang Hensel 实现的原有"首项系数须为常数"限制已由 0.16.0 的任意多元路径覆盖（走 EEZ）；二元非常数 LC 输入现自动转入 EEZ 路径。
- 自 0.16.2 起 $\mathbb{F}_p$ 多元路径同样支持非常数首项系数（域版 Wang LC 预处理 `wang_reconstruct_lcoeffs_fp` + `eez_lift_imposed`）。
- 代数数域因式分解（0.17.0）目前**仅限一元**；扩域上的多元因式分解与解析器的 `RootOf`/根式语法留待后续版本。极小多项式的不可约性由调用方保证（模可约时环中存在零因子，求逆会失败）。
- 自适应赋值点搜索已扩大（ℤ：7 → 25，$\mathbb{F}_p$：8 → 32）；对于非常稀疏或高度特殊的多项式，可能仍需进一步扩大搜索范围或额外的 Diophantine 插值遍次。

这些限制已在项目路线图中跟踪，并会随着代数内核的成熟逐步解除。
