//! Rubi 1892-case integration coverage harness.
//!
//! Single-round report-style benchmark (no criterion timing iterations; see
//! `tests/groebner_timing.rs` for the same style). Reads the deterministic
//! 1892-case sample produced by `ocas-tests/scripts/fetch_rubi_corpus.py`
//! and classifies every problem as solved or fallback (`Integral(...)`),
//! bucketed by integrand head so the largest fallback families are easy to
//! read off the report.
//!
//! Run after fetching the corpus:
//!
//! ```text
//! uv run python ocas-tests/scripts/fetch_rubi_corpus.py
//! cargo bench -p ocas-tests --bench integrate_1892
//! ```
//!
//! The report is written to `ocas-tests/data/integrate_1892_report.json` so
//! two runs (rules on/off, before/after) can be diffed. When the corpus is
//! missing the bench prints the fetch instructions and exits 0 (same
//! skip-by-convention as the uv-absent correctness tests).

use ocas::prelude::*;
use ocas_atom::{Atom, AtomArena, AtomNode};
use ocas_core::arena::Arena;
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

const TSV_PATH: &str = "data/rubi_1892.tsv";
const REPORT_PATH: &str = "data/integrate_1892_report.json";

/// Bucket names from the 0.27.0 plan (S1).
const BUCKETS: [&str; 9] = [
    "power-binomial",
    "exp-log",
    "trig",
    "hyperbolic",
    "inverse-trig-hyper",
    "radical",
    "rational",
    "special",
    "mixed-other",
];

fn contains_fractional_pow(atom: Atom<'_>) -> bool {
    match atom.node() {
        AtomNode::Pow(_, e) => !matches!(e.node(), AtomNode::Num(_)),
        AtomNode::Add(terms) | AtomNode::Mul(terms) => {
            terms.iter().any(|t| contains_fractional_pow(*t))
        }
        AtomNode::Fun(_, args) => args.iter().any(|a| contains_fractional_pow(*a)),
        AtomNode::Num(_) | AtomNode::Var(_) => false,
    }
}

/// Collect the distinct function heads appearing anywhere in the atom.
fn collect_heads(atom: Atom<'_>, out: &mut Vec<String>) {
    match atom.node() {
        AtomNode::Fun(name, args) => {
            out.push(name.as_str().to_string());
            for a in *args {
                collect_heads(*a, out);
            }
        }
        AtomNode::Add(terms) | AtomNode::Mul(terms) => {
            for t in *terms {
                collect_heads(*t, out);
            }
        }
        AtomNode::Pow(b, e) => {
            collect_heads(*b, out);
            collect_heads(*e, out);
        }
        AtomNode::Num(_) | AtomNode::Var(_) => {}
    }
}

/// Classify an integrand into one of the plan's nine buckets.
fn bucket(expr: Atom<'_>) -> &'static str {
    let mut heads = Vec::new();
    collect_heads(expr, &mut heads);
    heads.sort();
    heads.dedup();
    let head = heads.as_slice();
    if head.is_empty() {
        // Pure algebraic. Radicals win over binomial/power shapes.
        if contains_fractional_pow(expr) {
            return "radical";
        }
        match expr.node() {
            AtomNode::Pow(_, _) | AtomNode::Mul(_) => {
                // `x^n` / `(a+b*x)^n` / products of powers.
                "power-binomial"
            }
            AtomNode::Add(_) | AtomNode::Num(_) | AtomNode::Var(_) => "rational",
            // Unreachable: any Fun head would have made `heads` non-empty.
            AtomNode::Fun(_, _) => "mixed-other",
        }
    } else if head.len() == 1 {
        let name = head[0].as_str();
        if matches!(name, "sin" | "cos" | "tan" | "cot" | "sec" | "csc") {
            "trig"
        } else if matches!(name, "sinh" | "cosh" | "tanh" | "coth" | "sech" | "csch") {
            "hyperbolic"
        } else if matches!(
            name,
            "asin"
                | "acos"
                | "atan"
                | "acot"
                | "asec"
                | "acsc"
                | "asinh"
                | "acosh"
                | "atanh"
                | "acoth"
                | "asech"
                | "acsch"
        ) {
            "inverse-trig-hyper"
        } else if matches!(name, "exp" | "log") {
            "exp-log"
        } else if matches!(name, "sqrt") {
            "radical"
        } else if matches!(
            name,
            "erf" | "erfc" | "erfi" | "Ei" | "Si" | "Ci" | "Shi" | "Chi" | "fresnels" | "fresnelc"
        ) {
            "special"
        } else {
            "mixed-other"
        }
    } else {
        "mixed-other"
    }
}

/// Parse a TSV row `id\tintegrand\tvar`.
fn parse_row(line: &str) -> Option<(&str, &str, &str)> {
    let mut it = line.splitn(3, '\t');
    let id = it.next()?;
    let integrand = it.next()?;
    let var = it.next()?;
    Some((id, integrand, var))
}

/// Outcome of a single corpus problem.
enum CaseOutcome {
    /// Antiderivative found; no `Integral(...)` residue.
    Solved(&'static str),
    /// `Integral(...)` residue (or abandoned at the budget / crashed).
    Fallback(&'static str),
    /// The integrand did not parse.
    ParseErr,
    /// The child exited abnormally (panic, stack overflow, ...).
    Crashed,
    /// The child exceeded the per-case budget.
    TimedOut,
}

/// Run one corpus problem in a child process so stack overflows, infinite
/// loops, and panics cannot take the whole report down. The child writes one
/// line and exits; the parent polls `try_wait` so a hung child is abandoned
/// at the budget without ever blocking on the pipe.
fn run_case_in_child(integrand: &str, var: &str, rules: bool, budget_ms: u64) -> CaseOutcome {
    let exe = env::current_exe().expect("current exe");
    let mut child = std::process::Command::new(exe)
        .args(["--case", integrand, var])
        .env("OCAS_INTEGRATE_RULES", if rules { "1" } else { "0" })
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("spawn case child");
    let stdout = child.stdout.take().expect("child stdout");
    let deadline = Instant::now() + std::time::Duration::from_millis(budget_ms);
    loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            // Child exited: its single output line is already in the pipe.
            let mut out = String::new();
            use std::io::Read;
            let _ = std::io::BufReader::new(stdout).read_to_string(&mut out);
            if !status.success() {
                return CaseOutcome::Crashed;
            }
            return parse_child_line(out.trim_end());
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            return CaseOutcome::TimedOut;
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
}

/// Parse the single output line of a case child (`OK\t<bucket>\t<result>`
/// or `PARSE_ERR`).
fn parse_child_line(line: &str) -> CaseOutcome {
    if line == "PARSE_ERR" {
        return CaseOutcome::ParseErr;
    }
    let mut fields = line.splitn(3, '\t');
    match (fields.next(), fields.next(), fields.next()) {
        (Some("OK"), Some(b), Some(result)) => {
            let bucket = BUCKETS
                .iter()
                .copied()
                .find(|known| *known == b)
                .unwrap_or("mixed-other");
            if result.contains("Integral(") {
                CaseOutcome::Fallback(bucket)
            } else {
                CaseOutcome::Solved(bucket)
            }
        }
        _ => CaseOutcome::Crashed,
    }
}

fn main() {
    // Child mode: one problem per invocation, one line on stdout.
    let mut args = env::args().skip(1);
    if args.next().as_deref() == Some("--case") {
        let integrand = args.next().expect("case integrand");
        let var = args.next().expect("case var");
        // Rules path toggle: `OCAS_INTEGRATE_RULES=0` uses the pre-0.27
        // entry point for baseline comparison.
        let rules = env::var("OCAS_INTEGRATE_RULES").map_or(true, |v| v != "0");
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let expr = match parse(&ctx, &integrand) {
            Ok(e) => e,
            Err(_) => {
                println!("PARSE_ERR");
                return;
            }
        };
        let var_sym = Symbol::new(&var);
        let result = if rules {
            integrate(&ctx, expr, var_sym)
        } else {
            integrate_with_options(&ctx, expr, var_sym, IntegrateOptions { rules: false })
        };
        let b = bucket(expr);
        println!("OK\t{b}\t{result}");
        return;
    }

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let tsv = manifest_dir.join(TSV_PATH);
    if !tsv.exists() {
        println!(
            "integrate_1892: corpus not found at {}; fetch it first:\n\
             \n    uv run python ocas-tests/scripts/fetch_rubi_corpus.py\n",
            tsv.display()
        );
        return;
    }

    let rules_enabled = env::var("OCAS_INTEGRATE_RULES").map_or(true, |v| v != "0");
    // Per-problem wall-clock budget: a handful of corpus problems can take
    // minutes in the rational/Risch backends; they are counted as fallback
    // instead of stalling the report. Override with
    // `OCAS_INTEGRATE_CASE_TIMEOUT_MS`.
    let case_timeout_ms: u64 = env::var("OCAS_INTEGRATE_CASE_TIMEOUT_MS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(10_000);
    let corpus = fs::read_to_string(&tsv).expect("read corpus tsv");
    let mut solved = 0usize;
    let mut fallback = 0usize;
    let mut timed_out = 0usize;
    let mut crashed = 0usize;
    let mut parse_errors = 0usize;
    let mut counts: BTreeMap<&'static str, (usize, usize)> =
        BUCKETS.iter().map(|b| (*b, (0, 0))).collect();
    let start = Instant::now();

    for (row_idx, line) in corpus.lines().enumerate() {
        if row_idx % 200 == 0 {
            eprintln!("integrate_1892: {row_idx}/1892 cases");
        }
        let Some((_id, integrand, var)) = parse_row(line) else {
            eprintln!("integrate_1892: malformed row {}", row_idx + 1);
            continue;
        };
        // Each problem runs in a child process: a pathological case can
        // overflow the stack, loop forever, or panic without taking the
        // whole report down. The child prints one line and exits; the parent
        // abandons it at the budget and counts it as fallback.
        let outcome = run_case_in_child(integrand, var, rules_enabled, case_timeout_ms);
        match outcome {
            CaseOutcome::Solved(b) => {
                solved += 1;
                counts.get_mut(b).unwrap().0 += 1;
            }
            CaseOutcome::Fallback(b) => {
                fallback += 1;
                counts.get_mut(b).unwrap().1 += 1;
            }
            CaseOutcome::ParseErr => {
                parse_errors += 1;
            }
            CaseOutcome::Crashed => {
                crashed += 1;
                fallback += 1;
                counts.get_mut("mixed-other").unwrap().1 += 1;
            }
            CaseOutcome::TimedOut => {
                timed_out += 1;
                fallback += 1;
                counts.get_mut("mixed-other").unwrap().1 += 1;
            }
        }
    }

    let total = solved + fallback;
    let total_ms = start.elapsed().as_millis();
    let coverage = if total > 0 {
        100.0 * solved as f64 / total as f64
    } else {
        0.0
    };

    // Load corpus provenance from the sibling meta file when present.
    let meta_path = manifest_dir.join("data/rubi_1892.meta.json");
    let meta = fs::read_to_string(&meta_path)
        .ok()
        .and_then(|s| MiniJson::parse(&s).ok())
        .unwrap_or(MiniJson::Null);

    println!("integrate_1892 report");
    println!("  corpus:      {}", tsv.display());
    println!(
        "  rules:       {}",
        if rules_enabled { "on" } else { "off" }
    );
    println!("  solved:      {solved}");
    println!("  fallback:    {fallback}");
    println!("  timed out:   {timed_out}");
    println!("  crashed:     {crashed}");
    println!("  parse errs:  {parse_errors}");
    println!("  coverage:    {coverage:.2}% ({solved}/{total})");
    println!("  wall time:   {total_ms} ms");
    println!("  buckets:");
    for (b, (s, f)) in &counts {
        println!("    {b:>18}: solved {s:5} fallback {f:5}");
    }

    let buckets_obj = counts
        .iter()
        .map(|(k, (s, f))| {
            format!(
                "\"{k}\": {{\"solved\": {s}, \"fallback\": {f}, \"total\": {}}}",
                s + f
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let report = format!(
        "{{\n  \"seed\": {},\n  \"n\": {total},\n  \"requested\": {},\n  \"source_url\": {},\n  \
         \"source_sha256\": {},\n  \"digest_matched\": {},\n  \"timestamp\": \"{}\",\n  \
         \"rules_enabled\": {rules_enabled},\n  \"solved\": {solved},\n  \"fallback\": {fallback},\n  \
         \"timed_out\": {timed_out},\n  \"crashed\": {crashed},\n  \"parse_errors\": {parse_errors},\n  \"coverage_pct\": {coverage},\n  \"total_ms\": {total_ms},\n  \
         \"buckets\": {{{buckets_obj}}}\n}}",
        meta.get("seed"),
        meta.get("requested"),
        meta.get("source_url"),
        meta.get("source_sha256"),
        meta.get("digest_matched"),
        utc_now_rfc3339(),
    );
    let report_path = manifest_dir.join(REPORT_PATH);
    fs::write(&report_path, report).expect("write report json");
    println!("  report:      {}", report_path.display());
}

/// RFC 3339 UTC timestamp (avoids pulling a time crate into the bench).
fn utc_now_rfc3339() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("clock before epoch");
    let secs = now.as_secs() as i64;
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    // Civil-from-days (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mth = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mth <= 2 { y + 1 } else { y };
    format!("{y:04}-{mth:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

/// Minimal JSON for the small flat `rubi_1892.meta.json` file produced by
/// the fetch script (top-level object of scalars only).
///
/// Only the `Obj` variant's payload is read (`get`); the scalar variants are
/// parsed for structural validation but never inspected.
#[allow(dead_code)]
enum MiniJson {
    Null,
    Bool(bool),
    Num(String),
    Str(String),
    /// Raw text of the top-level object, scanned by [`MiniJson::get`].
    Obj(String),
}

impl MiniJson {
    fn parse(input: &str) -> std::result::Result<Self, String> {
        let mut p = MiniParser { s: input, pos: 0 };
        let v = p.value()?;
        p.ws();
        if p.pos != p.s.len() {
            return Err("trailing json".to_string());
        }
        Ok(v)
    }

    /// Render a top-level field value as JSON text; missing keys render as
    /// `null`.
    fn get(&self, key: &str) -> String {
        let MiniJson::Obj(text) = self else {
            return "null".to_string();
        };
        let needle = format!("\"{key}\":");
        let Some(idx) = text.find(&needle) else {
            return "null".to_string();
        };
        let rest = &text[idx + needle.len()..];
        let end = rest.find(['\n', ',', '}']).unwrap_or(rest.len());
        rest[..end].trim().to_string()
    }
}

struct MiniParser<'a> {
    s: &'a str,
    pos: usize,
}

impl<'a> MiniParser<'a> {
    fn ws(&mut self) {
        while let Some(c) = self.s[self.pos..].chars().next() {
            if c.is_whitespace() {
                self.pos += c.len_utf8();
            } else {
                break;
            }
        }
    }

    fn value(&mut self) -> std::result::Result<MiniJson, String> {
        self.ws();
        let rest = &self.s[self.pos..];
        if rest.starts_with("null") {
            self.pos += 4;
            Ok(MiniJson::Null)
        } else if rest.starts_with("true") {
            self.pos += 4;
            Ok(MiniJson::Bool(true))
        } else if rest.starts_with("false") {
            self.pos += 5;
            Ok(MiniJson::Bool(false))
        } else if rest.starts_with('"') {
            let mut i = 1;
            while i < rest.len() && !rest[i..].starts_with('"') {
                i += 1;
            }
            if i >= rest.len() {
                return Err("unterminated string".to_string());
            }
            let inner = &rest[1..i];
            self.pos += i + 1;
            Ok(MiniJson::Str(inner.to_string()))
        } else if rest.starts_with('{') {
            let mut depth = 0usize;
            let mut i = 0usize;
            let bytes = rest.as_bytes();
            while i < bytes.len() {
                match bytes[i] {
                    b'{' => depth += 1,
                    b'}' => {
                        depth -= 1;
                        if depth == 0 {
                            break;
                        }
                    }
                    b'"' => {
                        i += 1;
                        while i < bytes.len() && bytes[i] != b'"' {
                            if bytes[i] == b'\\' {
                                i += 1;
                            }
                            i += 1;
                        }
                    }
                    _ => {}
                }
                i += 1;
            }
            if depth != 0 {
                return Err("unbalanced object".to_string());
            }
            let raw = rest[..=i].to_string();
            self.pos += i + 1;
            Ok(MiniJson::Obj(raw))
        } else {
            // Number: digits with optional sign, dot, exponent.
            let end = rest
                .find(|c: char| {
                    !(c.is_ascii_digit()
                        || c == '-'
                        || c == '.'
                        || c == 'e'
                        || c == 'E'
                        || c == '+')
                })
                .unwrap_or(rest.len());
            if end == 0 {
                return Err("bad json value".to_string());
            }
            let num = rest[..end].to_string();
            self.pos += end;
            Ok(MiniJson::Num(num))
        }
    }
}
