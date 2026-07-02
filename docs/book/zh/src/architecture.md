# 架构

oCAS 由 12 个 crate 组成的 Cargo workspace，层间依赖严格向下。每个 crate 只能依赖更底层的 crate，不允许反向或循环依赖。

| 层级 | crate | 职责 |
|---|---|---|
| 5 绑定 | `ocas`、`ocas-py`、`ocas-c` | Rust、Python、C/C++ 公共 API |
| 4 应用 | `ocas-calc`、`ocas-eval`、`ocas-parse` | 微积分、求值、解析 |
| 3 符号引擎 | `ocas-atom`、`ocas-rewrite` | Atom、转换器、模式匹配、e-graph |
| 2 代数核 | `ocas-domain`、`ocas-poly` | 域、多项式、数论 |
| 1 数值后端 | `ocas-core` | GMP/FLINT 封装 |
| 0 运行时 | `ocas-core` | arena、错误、线程池、FFI |

完整设计文档见仓库中的 [ARCHITECTURE_CN.md](https://github.com/charleshzh/ocas/blob/main/ARCHITECTURE_CN.md)。
