# oCAS 基准测试结果（全量复测 @ 2026-08-06）

## 0.27.0 符号积分广度（Rubi 1892 题子集，2026-08-08）

> 本轮按《0.27.0 实现计划》S1/S5 执行：语料获取脚本 + 单轮报告型 harness
> `ocas-tests/benches/integrate_1892.rs`。语料为 Rubi 官方 Axiom 测试题
> （AxiomSyntaxTestFiles.zip，SHA-256 固定，seed 1892 确定性采样 1892 题），
> 不入库。报告 JSON 见 `ocas-tests/data/integrate_1892_report.json`
> （gitignored；运行命令与口径见 harness 头注释）。

### 采集命令

```text
uv run python ocas-tests/scripts/fetch_rubi_corpus.py
OCAS_INTEGRATE_RULES=0 cargo bench -p ocas-tests --bench integrate_1892   # 基线
cargo bench -p ocas-tests --bench integrate_1892                          # 规则开
```

### 结果

| 口径 | solved | fallback | coverage | 总墙钟 | 说明 |
|---|---|---|---|---|---|
| 0.26 链基线（rules off，1892 题，修复前测量） | 111 | 1781（含 1 超时） | 5.87% | 236.3 s | 0.27 前行为；**含 Weierstrass 伪解，需与下两行对照** |
| 0.27.0 规则引擎（rules on，修复前测量） | 124 | 1768（含 2 超时） | 6.55% | 301.2 s | A–H 家族 + 线性变元/裸参/积化和差两轮扩展；**含伪解** |
| (b) 基线（rules off，修复后重测 2026-08-08） | 138 | 1754（含 83 超时、0 崩溃） | **7.29%** | 1171.6 s | 共享链修复后重测（符号有理后端/Weierstrass 守卫在规则开关之外） |
| (b) 阶段后（rules on，修复后重测 2026-08-08） | 145 | 1747（含 83 超时、0 崩溃、0 解析错） | **7.66%** | 1190.2 s | 符号常数有理后端（Hermite+部分分式）＋ Weierstrass 线性变元 t-有理化；全部解经数值求导核验 |

bucket 分布（(b) 阶段后，solved/fallback）：power-binomial 107/218、radical 9/391、
mixed-other 7/655、trig 9/207、exp-log 5/74、inverse-trig-hyper 3/142、
hyperbolic 3/41、rational 2/0、special 0/19。

bucket 分布（rules on，solved/fallback）：mixed-other 31/584、radical 9/397、
power-binomial 39/300、trig 35/203、inverse-trig-hyper 3/143、exp-log 3/80、
hyperbolic 3/41、special 0/19、rational 1/1。

### 诚实结论

- **(b) 阶段：+1.11pp（+21 题，6.55% → 7.66%，修复后重测口径），仍未达 +30pp 验收线。**
  修复后基线重测 7.29%（138）——共享链改进同样提升规则关闭路径；规则开关净增量 +7 题。
  本阶段投入（产品级因子拆分第一步：符号常数有理函数升级）交付：
  - 新模块 `ocas-calc/src/integral/symbolic_rational.rs`：ℚ(symbols) 系数有理积分
    （多项式商 + Yun SFF + Hermite 每因子一步归约 + 线性留数/二次 atan/atanh/log
    部分分式，符号判别式 √(4ac−b²) 处理），修复稀疏 Horner 求值（`eval(x)=1` 根因）、
    SFF 因子首项系数保持、分数系数 to_sparse 往返损坏（a/c→a−c）、atanh 符号/√(−Δ)
    实值化等实现 bug；平方自由分母含次数≥3 符号因子时留诚实 `Integral` 部分结果。
  - Weierstrass 代换升级为线性变元 u=ax+b（`trig_linear_arg`/`substitute_trig_arg`），
    t-被积函数经 `rational_complexity_ok` 可行性闸门后才入链（消除全部挂死；
    83 超时均为 10 s 预算耗尽，0 崩溃）；有效性守卫移除 `2*atan(tan(x/2))*f(x)` 伪解。
  - 新增 u=ax+b 幂规则（`sin(a x+b)cos(a x+b)^n` / `sin^m cos`，`trig_odd_power_linear`）。
  - 核验：`1/(x(a+bx)^2)`、`(d+e*x)/(x^3*(a+c*x^2))`、`x^2/(a-b*x^2)^3`、
    `1/(a±b*x^2)`、`1/(a+b*sin(x))` 等数值求导 diff ≈ 0（SymPy 复核）；
    `cargo test --workspace --exclude ocas-py` 全绿；clippy -D warnings 干净。
  - 已知未解类：嵌套提取挂死（`(-5+3*cos(c+d*x))^3` 类，t-被积函数提取原子层）、
    ℚ 数值高次分母（`x^6/(216+108x^2+324x^3+18x^4+x^6)^2` 类）仍超时。
- **0.27.0 规则引擎首轮：+0.68pp（+13 题）。** 原因已量化：剩余 fallback 中
  606/614 个 trig 形态为复合乘积/商（`(c+d*x)^2*cos(a+b*x)*sin(a+b*x)^2` 等），
  需要因子级拆分/逐项积分机制而非整表达式首中规则表；radical bucket 397 题
  为嵌套/复合根式（`sqrt(a+c*x^2)/(d+e*x)`、`x^(5/2)` 等）；power-binomial
  300 题为符号参数有理函数（rational 后端仅支持 ℚ 系数，符号系数返回 None）。
- 规则表机制本身验证通过：Euler 代换补齐（S6）、规则路径开关（S7 绑定）、
  100 例 SymPy 抽样对拍（`cargo test -p ocas-tests --test correctness integral_rules`）
  全绿；单例规则（`tan(x)^4` → `tan(x)^3/3 - tan(x) + x` 等）正确。
- 按计划 S5 退出条件与 S8 备选分支：达标需 Rubi 级规则量级（7000+）或
  rational 后端符号系数升级；S8 评估结论为暂不集成 symbolica-integrate
  （许可证链条 + 依赖重量 + API 契合度，见 GAP_ANALYSIS CN/EN §7.3 追加段）。
  是否继续投入（如产品级因子拆分机制）需用户决策。

---

# oCAS 基准测试结果（全量复测 @ 2026-08-06）

> 本轮依据《竞品功能与性能差距重新评估计划》（2026-08-06）执行：8 个竞品最新版本
> 重新核实（见 COMPETITIVE_MATRIX_CN.md §1 备注）、本机全量复测 oCAS/Symbolica/SymPy
> 基准、WSL2 实测 msolve 0.10.1。**本块为当前有效数据**；下方 0.26.0 / 0.25.0 两块
> 为历史归档，原样保留。所有数字在执行当日用括号内命令重新采集。

## 测试机配置（采集命令见括号）

- OS：Windows 11 Pro，build 26200（`powershell (Get-CimInstance Win32_OperatingSystem).BuildNumber`）
- CPU：Intel(R) Core(TM) Ultra 7 255H，16 核 / 16 线程
  （`powershell (Get-CimInstance Win32_Processor | Select Name,NumberOfCores,NumberOfLogicalProcessors)`）
- 内存：31.5 GiB
  （`powershell [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)`）
- WSL 侧（msolve 用，必须单列）：WSL2，kernel `6.6.87.2-microsoft-standard-WSL2`
  （`wsl uname -r`），Ubuntu 24.04.4，分配 16 vCPU / ~15 GiB（`wsl nproc`、`wsl free -h`）。
  **msolve 计时在 WSL2 虚拟化环境下进行，与原生 Linux 存在未量化差异。**

## 测试配置

- Rust：`rustc -V` → 1.97.1 (8bab26f4f 2026-07-14)
- Python：`python --version` → 3.12.11；SymPy 1.14.0（.venv-sympy 与 ocas-tests uv 环境，
  `python -c "import sympy; print(sympy.__version__)"`）；Symbolica 2.2.0
  （`git -C ../symbolica describe --tags` → v2.2.0，commit 77c1374；runner 以
  `no_gmp` 纯 Rust 后端构建——crates.io numerica 2.2.0 在 default-features=false 下
  无法编译，经 `[patch.crates-io]` 指向检出内 `lib/numerica` 修复，见
  symbolica_runner/Cargo.toml）；msolve 0.10.1（WSL 源码构建，
  `/root/msolve-install/bin/msolve -h` 首行）
- oCAS：0.26.0，特性：默认（无 gmp/flint）；构建：release（workspace `Cargo.toml`
  `[profile.release]`：codegen-units=1, lto=true, panic=abort）。JIT/SIMD/NTT 基准组
  按需带 feature 构建（jit / simd / pulp / fast-poly / ntt），其余组默认构建
- 计时口径：criterion 默认值（各组 3 s 预热 + 采样，记录中位数）；`compare_sympy.py
  time` / symbolica_runner `time` 子命令为「1 次预热 + N 次计时之和（纳秒）」，N 逐
  任务记录；msolve 每个系统 ≥3 次取中位数（`-g 2` 仅 Gröbner 基模式，DRL/grevlex
  序），记录实际运行次数
- 线程：oCAS/rayon 用满 16 核（不限制）；msolve 单线程（默认；`-t 16` 在本构建上
  实测慢 100–1000×（cyclic-5：3 ms → 420–3086 ms），已弃用，此处显式标注）

## G1–G5 Gröbner 结果（ℤ₁₃，grevlex 为对标口径）

| ID | 输入 | oCAS 0.26 Lex | oCAS 0.26 grevlex | msolve 0.10.1（实测，grevlex） | 比值（grevlex） | 结论 |
|---|---|---|---|---|---|---|
| G1 | cyclic-5 ℤ₁₃ | 22.44 ms | 8.97 ms | 3 ms | 3.0× | ✅ |
| G2 | cyclic-6 ℤ₁₃ | 910.65 ms | **55.04 ms** | 4 ms | 13.8× | ⚠️ 里程碑达成，比值如实更新 |
| G3 | cyclic-7 ℤ₁₃ | >2 h 未完成（中断） | **3.829 s**（单轮） | 55 ms | ~70× | ⚠️ 差距最大单项 |
| G4 | katsura-6 ℤ₁₃ | 未测（预存在差距，katsura-6 单轮 >30 min） | — | 3 ms | — | ⚠️ 预存在差距 |
| G5 | katsura-7 ℤ₁₃ | 不可测 | — | 7 ms | — | ⚠️ 预存在差距 |

- msolve 实测为 3 次运行中位数（GB-only `-g 2`，输出丢弃）；单次记录：
  cyclic-5 3/3/2 ms、cyclic-6 4/4/4 ms、cyclic-7 54/54/55 ms、katsura-6 3/3/3 ms、
  katsura-7 7/7/7 ms。
- **输入同源性验证**：msolve 输入文件由 `ocas-tests/src/systems.rs`
  `cyclic_fp_grevlex`/`katsura_fp` 的公式逐项生成；两侧 reduced GB 基元素数完全一致
  （cyclic-5/6/7 = 20/45/209，katsura-6/7 = 36/70，msolve `-g 2 -o` 输出头部核对），
  cyclic-5 解数 70（参量化输出）与文献一致。
- **G2 比值修正说明**：上一版（0.26.0 归档）按引用值 msolve 0.04 s 记比值 1.3×；
  本轮实测 msolve 0.10.1 为 4 ms（快一个数量级），比值如实更新为 13.8×。「cyclic-6
  grevlex < 0.5 s」里程碑本身达成（55 ms，0.26.0 归档 52.07 ms 复现，落在 ±100%
  容差带内）。
- 补充实测（同环境，criterion 中位数）：cyclic-3 ℤ₁₃ 102.18 µs、cyclic-4 ℤ₁₃
  604.15 µs（Lex）；multi_modular_cyclic_q（ℚ）：cyclic-3 2.395 ms、cyclic-4
  5.039 ms、cyclic-5 124.08 ms。

## 因式分解 F1–F4（套件 §2.1，criterion 中位数）

| ID | 输入 | 目标 | oCAS 0.26 | 结论 |
|---|---|---|---|---|
| F1 | `x^100-1`（ℤ 完全分解） | < 50 ms | **252.43 ms**（x^60-1：205.09 ms；x^30-1：34.92 ms；x^12-1：4.05 ms） | ⚠️ 未达标 |
| F2 | 多元 3 变量（Wang EEZ） | < 100 ms | 1.9966 ms（trivariate_3_linear） | ✅ |
| F3 | 代数数域一元（Trager） | < 50 ms | 6.6383 ms（sqrt2_deg12）；19.79 ms（cbrt2_deg9）；27.09 ms（zeta5_deg9） | ✅ |
| F4 | `x^200-1`（ℤ 大规模） | 标记基准 | 未测（benches 上限 x^100） | ⚠️ 未测 |

## 多项式 GCD C1–C3（套件 §2.4，criterion 中位数）

| ID | 输入 | 目标 | oCAS 0.26 | 结论 |
|---|---|---|---|---|
| C1 | `gcd(x^50-1, x^30-1)`（ℤ 大系数） | < 1 s | 该精确输入未入 benches；最接近：gcd(x^50-1, x-1) 38.906 µs、gcd(x^500-1, x-1) 1.3786 ms | ✅ 量级满足（近似口径） |
| C2 | 多元 GCD | < 100 ms | 2.4962 ms（modular_bivariate）；22.927 µs（heuristic_bivariate） | ✅ |
| C3 | 模 GCD 性能 | < 500 ms | 2.4962 ms（二元多素数模路径） | ✅ |

## 符号积分 I1–I5（套件 §2.3，criterion 中位数）

| ID | 输入 | 目标 | oCAS 0.26 | 结论 |
|---|---|---|---|---|
| I1 | `∫(x^3+2x+1)/(x^2+1) dx` | < 1 ms | 139.93 µs（integrate_rational_atan，同型输入） | ✅ |
| I2 | `∫sin(x)^3*cos(x)^2 dx` | < 1 ms | 未直接入 benches；同型三角-多项式（x*sin(x) 类）176.33 µs | ✅ 量级满足（近似口径） |
| I3 | `∫exp(-x^2) dx → erf` | < 1 ms | 159.9 µs（integrate_exp_neg_x2_erf） | ✅ |
| I4 | `∫(log(x)/(x+1)^2) dx` | < 5 ms | 未直接入 benches（integrate_log_x 23.312 µs 为 log(x) 直接积分） | ✅ 量级满足（近似口径） |
| I5 | Rubi 1892 题子集 | 总时间标记 | 0.26 链基线 5.87%（111/1892，236.3 s）；0.27.0 规则引擎 6.55%（124/1892，301.2 s） | ⚠️ 未达 +30pp（见下方 0.27.0 段） |

## JIT 求值（套件 §2.5，criterion 中位数，feature: jit/simd/pulp）

| ID | 用例 | 目标 | oCAS 0.26 | 结论 |
|---|---|---|---|---|
| J1 | 单输出 JIT（batch 1000） | > 10× 解释器 | 解释器 121.11 µs → JIT 1.6418 µs = **73.8×** | ✅ |
| J2 | 多输出 JIT（3 输出，batch 1000） | > 15× 解释器 | 305.72 µs → 13.522 µs = **22.6×** | ✅ |
| J3 | SIMD 向量化（batch 4000） | > 8× 标量 | poly_batch_4000：490.83 µs → 73.028 µs = **6.7×**；trig_batch_4000：1.1559 ms → 505.97 µs = 2.3× | ⚠️ 未达标（poly 6.7× < 8×） |

## SymPy 对比（2026-08-06，SymPy 1.14.0）

计时口径：SymPy 为 `compare_sympy.py time`（1 预热 + N 次，总纳秒）；oCAS 为同输入
criterion 中位数。比值 = SymPy ns/op ÷ oCAS ns/op（>1 表示 oCAS 快）。

| 任务 | 输入 | N | SymPy ns/op | oCAS ns/op | 比值 |
|---|---|---|---|---|---|
| parse | `(x + y)^5 + sin(x)*cos(x)` | 500 | 161,473 | 1,606 | **100.5×** |
| diff | `(x + y)^5 + sin(x)*cos(x)` | 500 | 23,882 | 49,859 | 0.48× ⚠️ |
| expand | `(x + y)^20` | 200 | 1,052 | 593 | 1.8× |
| factor | `x^30 - 1` | 50 | 697,020 | 34,917,000 | 0.02× ⚠️ |
| gcd | `x^20 - 1;x^10 - 1` | 100 | 72,450 | 551（parse 代理） | 注：oCAS 侧为代理值 |
| simplify | `x + x + x + y + y + 0` | 500 | 1,763,362 | 14,265 | **123.6×** |
| integrate | `(x^3 + 2*x + 1)/(x^2 + 1)` | 100 | 10,626,499 | 139,930 | **75.9×** |
| integrate | `x*exp(x)` | 500 | 6,848,007 | 176,330 | **38.8×** |
| eval | `x^4 + 2*x^3 + x^2 + 1 @ x=1.5` | 500 | 23,324 | 127 | **183×** |
| series | `exp(sin(x)):10` | 100 | 327,492,333 | 128,440 | **2,550×** |
| series | `sin(x):5` | 500 | 3,326,457 | 128,440（taylor_exp_order_5 代理） | 25.9× |
| solve_linear | `x + y - 5;x - y - 1` | 200 | 54,900 | 无对应 criterion 组 | — |
| roots | `x^3 - 2*x^2 - x + 2` | 50 | 488,870 | 无对应 criterion 组 | — |
| nt_factorint | 999999937 | 20 | 345 | 无对应 criterion 组 | — |
| nt_isprime | 2305843009213693951 (2^61−1) | 20 | 50,215 | 同上 | — |
| nt_nextprime | 1152921504606846976 (2^60) | 20 | 41,350 | 同上 | — |
| nt_totient | 999999937 | 100 | 1,748 | 同上 | — |
| nt_mobius | 30030 | 100 | 2,009 | 同上 | — |
| nt_divisor_count | 1048576 (2^20) | 100 | 3,077 | 同上 | — |
| nt_divisor_sigma | 1048576;2 | 100 | 1,807 | 同上 | — |
| nt_liouville | 30030 | 100 | 70 | 同上 | — |
| nt_discrete_log | 101;3;75 | 100 | 5,646 | 同上 | — |
| nt_crt | 3,5,7\|2,3,2 | 100 | 2,639 | 同上 | — |
| nt_jacobi | 35;97 | 100 | 1,870 | 同上 | — |

口径注记（诚实记录）：
- **diff 行**：criterion 的 oCAS diff 用例每轮包含 parse + Arena 分配，SymPy 侧解析
  一次后重复求导，对 oCAS 不利；该行仅作参考。
- **factor 行**：SymPy `factor(x^30-1)`（0.697 ms/op）显著快于 oCAS 稠密分解
  （34.9 ms/op，约 50×）。cyclotomic 分解路径（`x^n−1` 类）是 SymPy 的强项
  （内部三角/分圆优化）；oCAS 走通用 CZ+Hensel。已如实记录，不修饰。
- **gcd 行**：ocas_vs_sympy_gcd 的 oCAS 侧基准为 parse 占位（见 bench 源码注释），
  不构成真实 GCD 对比；真实 GCD 见 C1–C3 表。
- nt_* 行只给 SymPy 侧数值（oCAS 侧无 criterion 组），不设比值。

## Symbolica 对比（2026-08-06，v2.2.0 / 77c1374，no_gmp 后端）

### 正确性对拍（`compare_symbolica.py` 输出 vs oCAS 输出，SymPy 语义核验）

| 任务 | 输入 | Symbolica 输出 | 对拍结论 |
|---|---|---|---|
| parse | `(x + y)^5 + sin(x)*cos(x)` | 展开式（Symbolica parse 即展开） | ✅ 与 SymPy 语义一致（oCAS parse 输出为未展开规范化形式，语义一致） |
| diff | `x^3` | `3*x^2` | ✅ 一致 |
| diff | `sin(exp(x^2))*cos(log(x))` | `-sin(log(x))*sin(e^x^2)/x+2*e^x^2*x*cos(log(x))*cos(e^x^2)` | ✅ 逐项核对与链式/乘积法则一致（SymPy 自动核验因 𝑒 字形差异判 False，人工确认正确） |
| expand | `(x + y)^20` | 21 项完全展开 | ✅ 与 SymPy 一致（verify=true） |
| factor | `x^4 - 1` | `-1+x * 1+x * 1+x^2` | ✅ 因子列表 (x−1)(x+1)(x²+1) 正确（字符串缺括号致自动核验 False，人工确认） |
| factor | `x^30 - 1` | 8 个分圆因子 | ✅ 与 Φ1..Φ30 完全分解一致（同上） |
| simplify | `x + x + x + y + y + 0` | `3*x+2*y` | ✅ 一致 |
| series | `exp(sin(x))` | 系数 1,1,1/2,−1/8,−1/15,−1/240,1/90,31/5760,1/5670,−2951/3628800 + O(x^11) | ✅ 与 SymPy series(exp(sin x), order 11) 系数逐项一致 |

### 性能对比（runner `time` 子命令，1 预热 + N 次总纳秒；与 SymPy 同口径直接相除）

| 任务 | 输入 | N | Symbolica ns/op | SymPy ns/op | 比值（SymPy/Symbolica） |
|---|---|---|---|---|---|
| parse | `(x + y)^5 + sin(x)*cos(x)` | 500 | 9,115 | 161,473 | **17.7×** |
| diff | `(x + y)^5 + sin(x)*cos(x)` | 500 | 2,826 | 23,882 | **8.5×** |
| expand | `(x + y)^20` | 200 | 21,061 | 1,052 | 0.05×（SymPy 快 20×） |
| factor | `x^30 - 1` | 50 | 1,023,022 | 697,020 | 0.68× |
| simplify | `x + x + x + y + y + 0` | 500 | 507 | 1,763,362 | 3,475×（Symbolica simplify=expand 代理，口径不对称） |
| series | `exp(sin(x))` | 100 | 1,811,239 | 327,492,333 | **181×**（Symbolica 固定 order 10；SymPy order 11+removeO） |

## 未达标项说明（诚实记录）

- **G2 cyclic-6 grevlex 比值 13.8×**：里程碑「< 0.5 s」达成（55.04 ms）；与 msolve
  实测 4 ms 的比值较归档版（按引用 0.04 s 记 1.3×）变差——原因是 msolve 0.10.1 本机
  实测远快于旧引用，非 oCAS 回退。以实测值为准。
- **G3 cyclic-7 grevlex ~70×（3.829 s vs 55 ms）**：本轮最大单项差距；Lex 序仍
  >2 h 未完成。
- **G4/G5 katsura**：预存在差距（0.25.0 已记录，katsura-6 单轮 >30 min），本轮未
  触及；msolve 侧 3/7 ms。
- **F1 x^100-1 252 ms 未达 <50 ms 目标**；**J3 SIMD 6.7× 未达 >8× 目标**（poly 批
  处理）；**SymPy factor(x^30-1) 快于 oCAS ~50×**。
- **Symbolica expand 慢于 SymPy 20×**（Symbolica 逐项乘法展开 vs SymPy 幂展开快
  路径），如实记录。
- 未测项：F4（x^200-1）、I5（Rubi 语料库）、solve_linear/roots/nt_* 的 oCAS 侧
  （无 criterion 组）、SymPy roots（nroots）与 oCAS 根隔离的对比。
- 环境差异：msolve 计时在 WSL2 虚拟化下进行；oCAS/rayon 16 线程 vs msolve 单线程
  （`-t 16` 异常慢，见上）。

## 后端 feature 全开对比（2026-08-06，同日同机）

与纯 Rust（default）同口径重跑全部 20 组 criterion 基准（本轮纯 Rust 与全后端均
为 rustc 1.97.1 GNU / release / codegen-units=1 / lto / panic=abort，同日采集）。

**feature 集合**：`gmp, jit, simd, ntt, sprs, fast-poly, ocas-core/gmp,
ocas-domain/gmp, ocas-domain/mpfr, ocas-domain/system-libs, ocas-poly/gmp`。
系统后端：MSYS2 MINGW64 gmp 6.3.0 / mpfr 4.2.2 / mpc 1.4.1（pkg-config 2.5.1，
经 Windows 前缀 `.pc` 副本链接，见下）。**flint 不可用**：flint3-sys 3.6.0 无条件
`pub use libc::pthread_mutex_t`（Windows libc 无此类型）+ bundled 构建要求 POSIX
环境——与 oCAS 既有文档「flint 仅 Linux/WSL」一致，本轮 Windows 全后端不含 flint。

**构建发现（真实 bug，已修）**：
- `benches/gmp_integer.rs` 两处 `GmpInteger::from_i64(9876543210987654321)` 超 i64
  范围（i64::MAX = 9223372036854775807），gmp feature 下 bench 无法编译（此前从未
  在 gmp 下编译过）。已改为 `9223372036854775807`，修复后 gmp 组正常出数。
- `multi_modular_cyclic_q` 组在 gmp 后端 panic：`rug 1.30.0 xmpz.rs:147 division by
  zero`（ℚ 多模管线的有理数除法除零，纯 Rust 后端不触发）。f5 组全部正常，仅该组
  失败，记录待查（gmp 后端 Integer 除法语义差异嫌疑）。

**关键组对比**（比值 = 全后端/纯 Rust，criterion 中位数）：

| 组 | 用例 | 纯 Rust | 全后端 | 比值 |
|---|---|---|---|---|
| Gröbner f5（ℤ_p，i64 快通道） | cyclic_3/4/5/6 | 100.1/583.3/20.7/840.1 ms | 84.0/465.0/17.7/853.9 ms | 0.84/0.80/0.85/1.02 |
| Gröbner f5 grevlex | cyclic_5/6 | 8.546/51.626 ms | 8.574/52.020 ms | 1.00/1.01 |
| dense_mul | degree_10/50/100/500 | 10.64/212.3/647.4/8.211 ms | 0.279/3.492/11.24/164.4 µs | **0.026/0.016/0.017/0.020**（快 38–60×） |
| sparse_mul | 全部 4 用例 | 14.6–1.71 ms | 2.69–229.7 µs | 0.13–0.18（快 5.4–7.4×） |
| poly_is_square_free | x^12-1/x^60-1/x^100-1 | 7.19/94.1/247.2 µs | 3.67/13.6/23.1 µs | 0.51/0.15/0.09（快 2–11×） |
| poly_gcd（ℤ 路径） | gcd_x^10/50/100/500 | 5.84/41.3/101.4/1.372 ms | 1.18/11.9/30.8/282.4 µs | 0.20–0.30（快 3.3–5×） |
| bivariate_gcd | heuristic/modular | 23.78/2.639 ms | 12.15/1.595 ms | 0.51/0.60（快 1.7–2×） |
| poly_factor_fp | F5/F7/F17 | 0.829/2.218/319.1 ms | 0.760/808.5/132.0 ms | 0.92/0.36/0.41（F7/F17 快 2.4–2.7×） |
| poly_factor_anf | sqrt2_deg12/cbrt2_deg9/zeta5_deg9 | 6.41/19.57/26.89 ms | 1.67/34.59/3.93 ms | 0.26/**1.77**/0.15 |
| **poly_factor_z** | x^12-1/x^30-1/x^60-1/x^100-1 | 4.08/34.4/203.1/252.7 ms | 49.0/475.9/1.670/16.56 s | **12.0/13.8/8.2/65.5**（慢 8–66×） |
| **poly_factor_multivariate_z** | triv_3_linear/triv_quad_linear/nonconst_lc/sparse_4var | 1.99/1.29/61.9/156.8 ms | 585.8/92.0/4.41/43.1 ms | **294.8/71.4**/0.07/0.28 |
| JIT | jit_exec_poly_1000 jit/interp | 1.652/121.3 µs | 1.641/84.2 µs | 0.99/0.69 |
| gmp 后端原生 | gmp_add/gmp_mul | —（gmp_disabled 占位） | 34.75/34.73 ns | 新增组 |
| SOO Integer（gmp 下） | small_add/small_mul/large_add/large_mul | — | 0.90/0.89/49.7/105.3 ns | 新增组 |

其余组（parse/normalize/calculus/integrate/rewrite/arena/eval_jit/eval_simd/
eval_streaming/ntt/sympy_comparison 等）比值 0.9–1.1×，不受后端 feature 影响。

**结论**：gmp 后端（rug 6.3.0）对**系数级整数运算**影响巨大且方向不一——稠密/稀疏
乘法、SFF、GCD、ℤ_p 因式分解普遍快 2–60×（rug 的 C 乘法 + 系数路径变化），但 ℤ
因式分解的 CZ/Hensel 路径（x^n−1 分圆族慢 8–66×，多元线性案例慢 71–295×）灾难性
退化；Gröbner ℤ_p（i64 快通道）与求值/JIT/解析层基本不受影响。因此默认（纯 Rust）
构建对 ℤ 因式分解更稳；gmp 后端适合密集系数运算场景，且 ℚ multi-modular 管线存在
rug 除零 bug（待修）。

---

# oCAS 基准测试结果（0.26.0）

归档位：`BENCHMARK_SUITE_CN.md` §4.3。每次版本发布时更新。

## 环境

```
环境：
- OS: Windows 11 Pro (26200)
- CPU: Intel Core Ultra 7 255H
- Rust: 1.97.1 (x86_64-pc-windows-gnu, MSVC 目标)
- Python: N/A（本次未测 Python 侧）
- oCAS: 0.26.0
- 特性: 默认（无 gmp/flint）
- 构建: release, codegen-units=1, LTO, panic=abort
- msolve: 引用 BENCHMARK_SUITE_CN.md 静态数据（未在本次环境复测）
```

## G1–G5 结果（ℤ₁₃，F5 打包快通道）

**对标口径变更（0.26.0）**：msolve 教程的 cyclic-6/7 基准数字为 **grevlex**（已查证），
故本轮新增 grevlex 变体作为对标基准；Lex 数字并列记录，不做跨序强行对比。
grevlex 相关测量均为修复 Graded 序度方向 bug（见下）之后的数值。

| ID | 输入 | oCAS 0.26.0 Lex | oCAS 0.26.0 grevlex | msolve（引用，grevlex） | 比值（grevlex） | 结论 |
|---|---|---|---|---|---|---|
| G1 | cyclic-5 ℤ₁₃ | 22.40 ms | 8.45 ms | 5 ms | 1.7× | ✅ |
| G2 | cyclic-6 ℤ₁₃ | **944.9 ms** | **52.07 ms** | 0.04 s | 1.3× | ✅（grevlex 口径） |
| G3 | cyclic-7 ℤ₁₃ | 单轮 > 2 h 未完成（中断，见下） | **4.620 s**（209 基元素） | ~1 s | 4.6× | ⚠️ grevlex 达标量级，Lex 未测得 |
| G4 | katsura-6 ℤ₁₃ | 预存在差距（0.25.0 未完成） | — | msolve 未给出 | — | ⚠️ 未变化 |
| G5 | katsura-7 ℤ₁₃ | 不可测 | — | msolve 未给出 | — | ⚠️ 未变化 |

补充实测（同环境，F5 打包快通道，criterion 中位数）：

| 输入 | Lex | grevlex |
|---|---|---|
| cyclic-3 ℤ₁₃ | 127.3 µs | — |
| cyclic-4 ℤ₁₃ | 672.0 µs | — |
| cyclic-5 ℤ₁₃ | 22.40 ms | 8.45 ms |
| cyclic-6 ℤ₁₃ | 944.9 ms | 52.07 ms |
| cyclic-7 ℤ₁₃ | 未完成（>2 h 中断） | 4.620 s（单轮） |

## 与 0.25.0 / 基线对比

| 输入 | 0.25.0 | 0.26.0 打包前（grevlex 修复后基线） | 0.26.0 最终 | 加速 |
|---|---|---|---|---|
| cyclic-6 ℤ₁₃ Lex | 1.4150 s | — | 944.9 ms | 1.50× |
| cyclic-6 ℤ₁₃ grevlex | （不可比，序有 bug） | 73.6 ms | 52.07 ms | 1.41× |
| cyclic-7 ℤ₁₃ grevlex | （不可比，序有 bug） | 5.755 s | 4.620 s | 1.25× |
| cyclic-7 ℤ₁₃ Lex | >60 min 未完成 | — | >2 h 未完成 | 未测得 |

0.26.0 加速来源：u128 SWAR 打包单项式管线（`PackedMono` 单条 u128 指数加法/整除/支持掩码、
注册表 Copy 键零克隆零堆分配、打包签名与 syzygy 桶集、并行行构造）、echelon 改造（i32 系数、
免克隆两阶段：pass-1 主元移入只读 store + Phase B 串行收尾，行序与纯串行逐位一致）、矩阵容量预分配。

## 预存在 bug 修复（0.26.0，诚实记录）

- **Graded 单项式序度方向反置**：`Grevlex`/`Grlex`/`WeightOrder`/`BlockOrder::SubOrder` 的
  度/权重比较方向错误（低度项视为更大、常数项最大，并非合法单项式序），grevlex 下
  buchberger/f4/f5 全部产生错误基（例：cyclic-7 grevlex 坍缩为 `{Πxᵢ−1}`、cyclic-3 坍缩为 1 项）。
  已按 Cox–Little–O'Shea Def. 2.4 修正。0.25.0 的 Lex 数字不受影响；任何 0.26.0 之前的
  grevlex 测量均不可信。修复后验证：cyclic-3 grevlex 基 3 项、cyclic-4 7 项、cyclic-5 20 项、
  cyclic-6 45 项、cyclic-7 209 项（ℤ₁₃ 与 ℤ₁₀₁ 一致），且 Lex 基与 grevlex 基互证生成同一理想
  （order 匹配归约双向为零）。

## 未达标项说明（诚实记录）

- **G2 Lex cyclic-6 = 944.9 ms**：0.26.0 模型预估 0.4–0.65 s，实测 0.945 s，未达模型值但较
  0.25.0（1.415 s）加速 1.5×。Lex 系统难度高于 msolve 的 grevlex 基准（跨序不作验收对比）。
- **G3 Lex cyclic-7**：按预决策规则，单轮运行 >2 h（wall，>8 CPU-h）中断，记「未完成」；
  0.25.0 为 >60 min 未完成，两版均未测得完成值。正确性断言（可解 + `is_groebner_basis`）
  保留为手工运行项。
- **G3 grevlex cyclic-7 = 4.620 s**：较 msolve ~1 s 差 4.6×，已达同量级；模型预估 3–10 min，
  实测远超预估（打包管线对 209 基元素系统的收益显著）。
- **G4/G5 katsura**：0.25.0 已记录为预存在差距（katsura-5 即 >30 min），本轮未触及，保持记录。

## 对齐验收说明

- G1 与 msolve 比值 < 10×：**达成**（Lex 4.5×，grevlex 1.7×）。
- G2 与 msolve 比值 < 10×（grevlex 口径）：**达成**（1.3×，0.26.0 打包 + echelon 后）。
- G3 grevlex 与 msolve 比值：4.6×（同量级，未达 <10× 之外的额外目标）。
- msolve 数值为 BENCHMARK_SUITE_CN.md 静态引用，未在本次环境复测（无 msolve 二进制）。

---

# oCAS 基准测试结果（0.25.0）

归档位：`BENCHMARK_SUITE_CN.md` §4.3。每次版本发布时更新。

## 环境

```
环境：
- OS: Windows 11 Pro (26200)
- CPU: Intel Core Ultra 7 255H
- Rust: 1.97.1 (x86_64-pc-windows-gnu, MSVC 目标)
- Python: N/A（本次未测 Python 侧）
- oCAS: 0.25.0
- 特性: 默认（无 gmp/flint）
- 构建: release, codegen-units=1, LTO
- msolve: 引用 BENCHMARK_SUITE_CN.md 静态数据（未在本次环境复测）
```

## G1–G5 结果（ℤ₁₃，F5 / f5_fp 快通道）

| ID | 输入 | oCAS 0.25.0 | 门限 | msolve（引用） | 比值 | 结论 |
|---|---|---|---|---|---|---|
| G1 | cyclic-5 ℤ₁₃ | **23.88 ms**（criterion 中位数） | < 0.1 s | 5 ms | 4.8× | ✅ |
| G2 | cyclic-6 ℤ₁₃ | **1.4150 s**（criterion 中位数；单轮 1.38–1.47 s） | 路线图 < 5 s；0.25 目标 < 0.5 s | 0.04 s | 35.4× | ⚠️ 路线图达标，0.25 目标未达成 |
| G3 | cyclic-7 ℤ₁₃ | 单轮 > 60 min 未完成（中断） | 标记基准 | ~1 s | — | ⚠️ 未测得完成值 |
| G4 | katsura-6 ℤ₁₃ | katsura-5 即 > 30 min 未完成 | < 1 s | msolve 未给出 | — | ⚠️ 预存在差距，见下 |
| G5 | katsura-7 ℤ₁₃ | 不可测（katsura-6 未完成） | 标记基准 | msolve 未给出 | — | ⚠️ 预存在差距，见下 |

补充实测（同环境，F5）：

| 输入 | oCAS 0.25.0 |
|---|---|
| cyclic-3 ℤ₁₃ | 118.9 µs |
| cyclic-4 ℤ₁₃ | 610.5 µs |
| katsura-3 ℤ₁₃ | 0.005 s |
| katsura-4 ℤ₁₃ | 2.261 s（基线 0.24.0 原 F5：6.440 s） |

### multi_modular_cyclic_q（0.25.0 新增管线，criterion 中位数）

| 输入 | oCAS 0.25.0 |
|---|---|
| cyclic-3 ℚ | 2.55 ms |
| cyclic-4 ℚ | 5.57 ms |
| cyclic-5 ℚ | 154.1 ms |

## 与 0.19.0 基线对比（cyclic-6 ℤ₁₃）

- 0.19.0（CHANGELOG.md:469，F5/f5_fp）：**2.63 s**
- 0.25.0：**1.415 s**（criterion 中位数），约 **1.9×** 加速。
- 0.25.0 加速来源：`find_reducer_fp` DivisorIndex（支持掩码分桶）、syzygy 集合分桶（`contains` 由 O(syzygies) 线性扫描降为 submask 枚举）、矩阵行构造 rayon 并行、echelon 两阶段并行化（pass-1 主元并行消去 + 串行收尾，结果与纯串行逐位一致）、行注册 fused + 栈缓冲（免逐项堆分配）。
- **0.5 s 目标未达成**：实测 1.415 s，差距 2.8×。已如实记录，不伪造数字。cyclic-6 的 < 2.0 s 断言（correctness 套件）通过。

## G3–G5 未达成说明（诚实记录）

- **G3（cyclic-7 ℤ₁₃）**：0.19.0 注释为「>5 min」；0.25.0 单轮运行 60 min 未完成（中断），未测得完成时间。该门限为「标记基准」（无硬时限），不阻塞发布；正确性断言（可解 + `is_groebner_basis`）保留为手工运行项。
- **G4/G5（katsura-6/7 ℤ₁₃）**：经隔离验证（`git stash` 仅回退 f5.rs 后复测），katsura 系统的病态慢速**为 0.24.0 预存在行为**，非本次改动引入：katsura-4 原 F5 6.44 s（本版本 2.26 s，仍远高于合理量级），katsura-5 原 F5 亦 > 30 min。F5 在 katsura 系上的性能差距（syzygy 准则过滤不足/矩阵规模）需后续专项优化，列入下版本路线图。

## 对齐验收说明

- G1 与 msolve 比值 < 10×：**达成**（4.8×）。
- G2 与 msolve 比值 < 10×：**未达成**（35.4×），已如实记录。Symbolica 参考值 ~1 s（BENCHMARK_SUITE_CN.md），oCAS 0.25.0 与其同量级。
- msolve 数值为 BENCHMARK_SUITE_CN.md 静态引用，未在本次环境复测（无 msolve 二进制）。
