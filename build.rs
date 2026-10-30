use std::env;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn run(command: &mut Command) {
    let output = command
        .output()
        .expect("could not launch native build tool");
    if !output.status.success() {
        panic!(
            "{command:?} failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

fn main() {
    for path in ["CMakeLists.txt", "native", "dependencies/iqtree2"] {
        println!("cargo:rerun-if-changed={path}");
    }
    for name in [
        "CMAKE",
        "CC",
        "CXX",
        "AR",
        "EIGEN3_INCLUDE_DIR",
        "CMAKE_PREFIX_PATH",
    ] {
        println!("cargo:rerun-if-env-changed={name}");
    }

    assert_eq!(
        env::var("HOST").unwrap(),
        env::var("TARGET").unwrap(),
        "cross-compilation is not configured for IQ-TREE"
    );
    let os = env::var("CARGO_CFG_TARGET_OS").unwrap();
    assert!(
        os == "macos" || os == "linux",
        "native build supports macOS and Linux"
    );

    let root = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    let native = out.join("native");
    let cmake = env::var_os("CMAKE").unwrap_or_else(|| "cmake".into());
    let mut configure = Command::new(&cmake);
    configure.arg("-S").arg(&root).arg("-B").arg(&native).args([
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON",
    ]);
    for name in ["EIGEN3_INCLUDE_DIR", "CMAKE_PREFIX_PATH"] {
        if let Some(value) = env::var_os(name) {
            let mut argument = std::ffi::OsString::from(format!("-D{name}="));
            argument.push(value);
            configure.arg(argument);
        }
    }
    run(&mut configure);
    run(Command::new(&cmake)
        .arg("--build")
        .arg(&native)
        .args(["--target", "pf_alisim_bridge", "--parallel"])
        .arg(env::var("NUM_JOBS").unwrap_or_else(|_| "2".into())));

    let archives = fs::read_to_string(native.join("native-archives.txt")).unwrap();
    let bundle = out.join("libpf_alisim.a");
    if os == "macos" {
        let mut tool = Command::new("/usr/bin/libtool");
        tool.args(["-static", "-o"]).arg(&bundle);
        for archive in archives.lines() {
            tool.arg(Path::new(archive));
        }
        run(&mut tool);
    } else {
        // MRI merging preserves duplicate object basenames across archives.
        let mut script = format!("CREATE \"{}\"\n", bundle.display());
        for archive in archives.lines() {
            script.push_str(&format!("ADDLIB \"{archive}\"\n"));
        }
        script.push_str("SAVE\nEND\n");
        let mut child = Command::new(env::var_os("AR").unwrap_or_else(|| "ar".into()))
            .arg("-M")
            .stdin(Stdio::piped())
            .spawn()
            .expect("could not launch ar");
        child
            .stdin
            .take()
            .unwrap()
            .write_all(script.as_bytes())
            .unwrap();
        assert!(
            child.wait().unwrap().success(),
            "could not bundle native archives"
        );
    }

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=pf_alisim");
    println!(
        "cargo:rustc-link-lib=dylib={}",
        if os == "macos" { "c++" } else { "stdc++" }
    );
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=m");
    if os == "linux" {
        println!("cargo:rustc-link-lib=pthread");
    }
}
