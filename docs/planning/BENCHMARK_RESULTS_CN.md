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
