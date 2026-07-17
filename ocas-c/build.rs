use std::env;
use std::path::PathBuf;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Write the generated header into OUT_DIR, never into the source tree.
    // `cargo publish` verifies the packaged tarball is not modified by the
    // build script; writing into `include/` (the source dir) aborts publish
    // with "Source directory was modified by build.rs". The committed
    // `include/ocas.h` is the canonical copy shipped in the package; OUT_DIR
    // is only for building from source.
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let output = out_dir.join("ocas.h");

    cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(cbindgen::Config::from_root_or_default(&crate_dir))
        .generate()
        .expect("unable to generate C header")
        .write_to_file(&output);

    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=cbindgen.toml");
}
