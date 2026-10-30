use std::fs;

use pf_simulator::alignment::{
    AlignmentConfig, SequenceType, simulate_alignment_from_tree, simulate_alignment_numeric,
};
use pf_simulator::tree::Tree;

#[test]
fn direct_tree_matches_file_input_without_creating_files() {
    let directory = std::env::temp_dir().join(format!("pf-tree-input-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let prefix = directory.join("missing").join("alignment");
    let filename = directory.join("reference.nwk");
    for tips in [3, 10, 200] {
        for tree_seed in [1, 42] {
            let tree = Tree::uniform(tips, 0.1, tree_seed).unwrap();
            for model in ["LG", "LG+G4{0.5}", "LG+GC{0.5}", "LG+H2", "JC"] {
                for with_indels in [false, true] {
                    let config = AlignmentConfig {
                        length: 100,
                        model: model.into(),
                        sequence_type: if model == "JC" {
                            SequenceType::Dna
                        } else {
                            SequenceType::AminoAcid
                        },
                        indel: with_indels.then(|| "0.1,0.05".into()),
                        indel_size: with_indels.then(|| "GEO{2},GEO{3}".into()),
                        ..Default::default()
                    };
                    // The entire directory is empty while direct simulation runs.
                    assert_eq!(fs::read_dir(&directory).unwrap().count(), 0);
                    let direct =
                        simulate_alignment_from_tree(&tree, &prefix, &config, 777, None).unwrap();
                    assert_eq!(fs::read_dir(&directory).unwrap().count(), 0);
                    assert_eq!(direct.rows.len(), tips);

                    fs::write(&filename, tree.newick()).unwrap();
                    let reference =
                        simulate_alignment_numeric(&filename, &prefix, &config, 777, None).unwrap();
                    fs::remove_file(&filename).unwrap();
                    assert_eq!(direct.names, reference.names, "{tips} tips, {model}");
                    assert_eq!(direct.alphabet, reference.alphabet);
                    assert_eq!(direct.rows, reference.rows, "{tips} tips, {model}");
                    assert_eq!(direct.export_unaligned, reference.export_unaligned);
                    assert_eq!(direct.has_duplicates(), reference.has_duplicates());
                }
            }
        }
    }
    fs::remove_dir_all(directory).unwrap();
}
