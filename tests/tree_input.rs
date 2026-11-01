use std::fs;

use pf_simulator::alignment::{
    AlignmentConfig, SequenceType, simulate_alignment_from_tree, simulate_alignment_numeric,
};
use pf_simulator::tree::Tree;
use pf_simulator::tree_simulation::{TreeArgs, TreePrior, TreeSimulation};

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

#[test]
fn rooted_tree_priors_match_newick_input() {
    let dir = std::env::temp_dir().join(format!("pf-rooted-input-{}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let path = dir.join("tree.nwk");
    let prefix = dir.join("missing").join("alignment");
    for prior in [
        TreePrior::UniformRooted,
        TreePrior::Ete3,
        TreePrior::Yule,
        TreePrior::BirthDeath,
    ] {
        let generator = TreeSimulation::new(&TreeArgs {
            tree_prior: Some(prior),
            ..Default::default()
        })
        .unwrap();
        for tips in [3, 20, 200] {
            let tree = generator.generate(tips, 773).unwrap();
            assert_eq!(tree.edges().len(), 2 * tips - 2);
            for model in ["LG", "LG+H2"] {
                for indels in [false, true] {
                    let config = AlignmentConfig {
                        length: 100,
                        model: model.into(),
                        indel: indels.then(|| "0.01,0.01".into()),
                        ..Default::default()
                    };
                    let direct =
                        simulate_alignment_from_tree(&tree, &prefix, &config, 42, None).unwrap();
                    fs::write(&path, tree.newick()).unwrap();
                    let file =
                        simulate_alignment_numeric(&path, &prefix, &config, 42, None).unwrap();
                    assert_eq!(direct.rows.len(), tips);
                    assert_eq!(direct.names, file.names);
                    assert_eq!(direct.rows, file.rows);
                    assert_eq!(direct.alphabet, file.alphabet);
                    fs::remove_file(&path).unwrap();
                }
            }
        }
    }
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 0);
    fs::remove_dir(dir).unwrap();
}
