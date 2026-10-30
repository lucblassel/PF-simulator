use std::path::Path;
use std::process::Command;

#[test]
fn executable_contract() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let python = std::env::var_os("PYTHON").unwrap_or_else(|| "python3".into());
    let result = Command::new(python)
        .arg(root.join("tests/test_simulation.py"))
        .arg(env!("CARGO_BIN_EXE_pf-sim"))
        .output()
        .expect("end-to-end tests require Python 3 (or set PYTHON)");
    assert!(
        result.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}
