use pf_simulator::{
    alignment::{AlignmentConfig, simulate_alignment_from_tree, simulate_alignment_in_memory},
    tree::Tree,
};
use std::{collections::BTreeMap, fs};

fn rows(bytes: &[u8]) -> BTreeMap<String, String> {
    let mut rows = BTreeMap::new();
    let mut name = String::new();
    for line in std::str::from_utf8(bytes).unwrap().lines() {
        if let Some(label) = line.strip_prefix('>') {
            name = label.trim().into();
            rows.insert(name.clone(), String::new());
        } else {
            rows.get_mut(&name).unwrap().push_str(line);
        }
    }
    rows
}

#[test]
fn ancestral_numeric_matches_reference_for_all_rows_and_preserves_tips() {
    let directory = std::env::temp_dir().join(format!("pf-ancestral-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    for rooted in [false, true] {
        let source = if rooted {
            "((T1:0.1,T2:0.2):0.3,(T3:0.4,T4:0.5):0.6);"
        } else {
            "(T1:0.1,T2:0.2,(T3:0.3,T4:0.4):0.5);"
        };
        let tree =
            Tree::from_phylogeny(phylotree::tree::Tree::from_newick(source).unwrap()).unwrap();
        let file = directory.join("tree.nwk");
        fs::write(&file, tree.newick_with_internal_names().unwrap()).unwrap();
        for indels in [false, true] {
            let config = AlignmentConfig {
                model: "LG+GC{0.5}".into(),
                length: 100,
                indel: indels.then(|| "0.02,0.01".into()),
                indel_size: indels.then(|| "GEO{2},GEO{3}".into()),
                ancestral: true,
                ..Default::default()
            };
            let reference = simulate_alignment_in_memory(
                &file,
                &directory.join("reference"),
                &config,
                42,
                None,
            )
            .unwrap();
            let raw = simulate_alignment_from_tree(
                &tree,
                &directory.join("missing/raw"),
                &config,
                42,
                None,
            )
            .unwrap();
            assert_eq!(raw.rows.len(), if rooted { 7 } else { 6 });
            assert_eq!(rows(&raw.to_fasta(false)), rows(&reference.aligned_fasta));
            if let Some(unaligned) = reference.unaligned_fasta {
                assert_eq!(rows(&raw.to_fasta(true)), rows(&unaligned));
            }
            assert!(!directory.join("missing").exists());
            let tips = simulate_alignment_from_tree(
                &tree,
                &directory.join("missing/tips"),
                &AlignmentConfig {
                    ancestral: false,
                    ..config.clone()
                },
                42,
                None,
            )
            .unwrap();
            let all = rows(&raw.to_fasta(true));
            assert_eq!(tips.rows.len(), 4);
            for (name, row) in rows(&tips.to_fasta(true)) {
                assert_eq!(all[&name], row);
            }
        }
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn unnamed_newick_internal_nodes_use_the_reference_ids() {
    use pf_simulator::alignment::simulate_alignment_numeric;
    let directory =
        std::env::temp_dir().join(format!("pf-ancestral-unnamed-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let file = directory.join("tree.nwk");
    fs::write(&file, "((T1:0.1,T2:0.2):0.3,(T3:0.4,T4:0.5):0.6);").unwrap();
    for indels in [false, true] {
        let config = AlignmentConfig {
            ancestral: true,
            length: 100,
            indel: indels.then(|| "0.02,0.01".into()),
            ..Default::default()
        };
        let reference =
            simulate_alignment_in_memory(&file, &directory.join("reference"), &config, 42, None)
                .unwrap();
        let raw =
            simulate_alignment_numeric(&file, &directory.join("missing/raw"), &config, 42, None)
                .unwrap();
        assert_eq!(rows(&raw.to_fasta(false)), rows(&reference.aligned_fasta));
        assert!(!directory.join("missing").exists());
    }
    fs::remove_dir_all(directory).unwrap();
}
