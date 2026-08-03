# 外部参考

- [Rust API 文档（docs.rs）](https://docs.rs/ocas)：公共 Rust API 的自动生成
  参考文档，由已发布 crate 以**可移植特性**构建（不含 `gmp`/`mpfr`/`flint`
  等需要系统 C 库的后端 feature）；如需完整后端 API 文档，可本地运行
  `cargo doc -p ocas --features gmp,mpfr,flint --no-deps`。
- [crates.io 注册页](https://crates.io/crates/ocas)：已发布版本、feature 标签
  与依赖图。
- [GitHub 仓库](https://github.com/charleshzh/ocas)：源代码、Issue 跟踪与贡献
  指南。
- [文档站点](https://charleshzh.github.io/ocas/latest/)：本书，由 GitHub Pages
  提供。每个发布 tag 下 `/ocas/v<版本>/` 存有对应的版本快照。
