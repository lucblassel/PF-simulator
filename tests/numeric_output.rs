use std::collections::BTreeMap;
use std::fs;

use pf_simulator::alignment::{
    AlignmentConfig, SequenceType, has_duplicate_sequences, simulate_alignment_in_memory,
    simulate_alignment_numeric,
};

fn rows(fasta: &[u8]) -> BTreeMap<Vec<u8>, Vec<u8>> {
    let mut rows = BTreeMap::new();
    let mut name = Vec::new();
    for line in fasta.split(|&byte| byte == b'\n').map(<[u8]>::trim_ascii) {
        if let Some(header) = line.strip_prefix(b">") {
            name = header.trim_ascii().to_vec();
            rows.insert(name.clone(), Vec::new());
        } else if !line.is_empty() {
            rows.get_mut(&name).unwrap().extend_from_slice(line);
        }
    }
    rows
}

#[test]
fn numeric_output_matches_reference_fasta_sequences() {
    let directory = std::env::temp_dir().join(format!("pf-raw-test-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let tree = directory.join("tree.nwk");
    fs::write(&tree, "(T1:0.5,T2:0.5,(T3:0.3,T4:0.3):0.2);\n").unwrap();
    let prefix = directory.join("missing").join("alignment");
    for seed in [42, 777] {
        for model in ["LG", "LG+G4", "LG+I{0.2}", "LG+GC{0.5}", "JC"] {
            for indels in [true, false] {
                let config = AlignmentConfig {
                    length: 100,
                    model: model.into(),
                    sequence_type: if model.starts_with("JC") {
                        SequenceType::Dna
                    } else {
                        SequenceType::AminoAcid
                    },
                    indel: indels.then(|| "0.1,0.05".into()),
                    indel_size: indels.then(|| "GEO{2},GEO{3}".into()),
                    no_unaligned: false,
                    ancestral: false,
                };
                let fasta =
                    simulate_alignment_in_memory(&tree, &prefix, &config, seed, None).unwrap();
                let raw = simulate_alignment_numeric(&tree, &prefix, &config, seed, None).unwrap();
                assert_eq!(raw.rows.len(), 4);
                assert_eq!(rows(&fasta.aligned_fasta), rows(&raw.to_fasta(false)));
                assert_eq!(
                    has_duplicate_sequences(&fasta.aligned_fasta).unwrap(),
                    raw.has_duplicates()
                );
                assert_eq!(raw.export_unaligned, fasta.unaligned_fasta.is_some());
                if let Some(unaligned) = fasta.unaligned_fasta {
                    assert_eq!(rows(&unaligned), rows(&raw.to_fasta(true)));
                }
                assert!(!prefix.parent().unwrap().exists());
            }
        }
    }
    // Exercise the larger-tree path and zero-rate/no-unaligned metadata.
    let large_tree = pf_simulator::tree::Tree::uniform(200, 0.1, 17).unwrap();
    fs::write(&tree, large_tree.newick()).unwrap();
    for (indel, no_unaligned) in [
        (None, false),
        (Some("0,0"), false),
        (Some("0.01,0.01"), false),
        (Some("0.01,0.01"), true),
    ] {
        let config = AlignmentConfig {
            length: 500,
            indel: indel.map(str::to_owned),
            no_unaligned,
            ..Default::default()
        };
        let fasta = simulate_alignment_in_memory(&tree, &prefix, &config, 123, None).unwrap();
        let raw = simulate_alignment_numeric(&tree, &prefix, &config, 123, None).unwrap();
        assert_eq!(raw.rows.len(), 200);
        assert_eq!(rows(&fasta.aligned_fasta), rows(&raw.to_fasta(false)));
        assert_eq!(
            has_duplicate_sequences(&fasta.aligned_fasta).unwrap(),
            raw.has_duplicates()
        );
        assert_eq!(raw.export_unaligned, fasta.unaligned_fasta.is_some());
        if let Some(unaligned) = fasta.unaligned_fasta {
            assert_eq!(rows(&unaligned), rows(&raw.to_fasta(true)));
        }
    }
    fs::remove_dir_all(directory).unwrap();
}
