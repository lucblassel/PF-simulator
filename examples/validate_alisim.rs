//! Emit unconditional alignments for an external standalone-AliSim comparison.
//! TSV rows: id, tree, model, alphabet, root length, seed, indel rates, indel sizes.
use std::error::Error;
use std::fs;
use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use pf_simulator::alignment::{
    AlignmentConfig, SequenceType, simulate_alignment_from_tree, simulate_alignment_in_memory,
    simulate_alignment_numeric,
};
use pf_simulator::tree::Tree;

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    Direct,
    File,
    Fasta,
}

#[derive(Parser)]
struct Args {
    #[arg(long)]
    manifest: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, value_enum, default_value = "direct")]
    mode: Mode,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    fs::create_dir(&args.output)?;
    let manifest = fs::read_to_string(&args.manifest)?;
    let mut count = 0;
    for line in manifest
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
    {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 8 {
            return Err("Manifest must have exactly eight TSV fields".into());
        }
        let [id, path, model, alphabet, length, seed, indel, size] = fields.as_slice() else {
            unreachable!()
        };
        if !id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
            || id.is_empty()
        {
            return Err("Unsafe fixture identifier".into());
        }
        let tree_file = PathBuf::from(path);
        let prefix = args.output.join(id);
        let config = AlignmentConfig {
            model: (*model).into(),
            length: length.parse()?,
            sequence_type: match *alphabet {
                "AA" => SequenceType::AminoAcid,
                "DNA" => SequenceType::Dna,
                _ => return Err("Unknown alphabet".into()),
            },
            indel: (!indel.is_empty()).then(|| (*indel).to_owned()),
            indel_size: (!size.is_empty()).then(|| (*size).to_owned()),
            no_unaligned: false,
            ancestral: false,
        };
        let seed = seed.parse()?;
        let output = match args.mode {
            Mode::Direct => {
                let tree = Tree::from_phylogeny(phylotree::tree::Tree::from_file(&tree_file)?)?;
                simulate_alignment_from_tree(&tree, &prefix, &config, seed, None)?.into_fasta()
            }
            Mode::File => {
                simulate_alignment_numeric(&tree_file, &prefix, &config, seed, None)?.into_fasta()
            }
            Mode::Fasta => simulate_alignment_in_memory(&tree_file, &prefix, &config, seed, None)?,
        };
        fs::write(args.output.join(format!("{id}.fa")), output.aligned_fasta)?;
        if let Some(unaligned) = output.unaligned_fasta {
            fs::write(args.output.join(format!("{id}.unaligned.fa")), unaligned)?;
        }
        count += 1;
    }
    println!("Wrote {count} unconditional alignment fixtures");
    Ok(())
}
