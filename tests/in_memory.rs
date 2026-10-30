use std::fs;
use std::sync::atomic::{AtomicU64, Ordering};

use pf_simulator::alignment::{AlignmentConfig, simulate_alignment_in_memory};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

#[test]
fn native_streams_return_fasta_without_creating_files() {
    let directory = std::env::temp_dir().join(format!(
        "pf-memory-{}-{}",
        std::process::id(),
        NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&directory).unwrap();
    let tree = directory.join("input.nwk");
    fs::write(&tree, "(T1:0.3,T2:0.3,T3:0.3);\n").unwrap();
    // This parent deliberately does not exist: no native output can be written
    // here, including the temporary-data file used by AliSim's indel algorithm.
    let prefix = directory.join("missing").join("alignment");
    for with_indels in [true, false, true] {
        let config = AlignmentConfig {
            length: 100,
            model: "LG+G4{0.5}".into(),
            indel: with_indels.then(|| "0.1,0.05".into()),
            indel_size: with_indels.then(|| "GEO{2},GEO{3}".into()),
            ..Default::default()
        };
        let output = simulate_alignment_in_memory(&tree, &prefix, &config, 42, None).unwrap();
        assert!(output.aligned_fasta.starts_with(b">T"));
        assert_eq!(output.unaligned_fasta.is_some(), with_indels);
        let entries: Vec<_> = fs::read_dir(&directory)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect();
        assert_eq!(entries, ["input.nwk"]);
    }
    fs::remove_dir_all(directory).unwrap();
}
