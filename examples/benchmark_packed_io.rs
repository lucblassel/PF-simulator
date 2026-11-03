//! Exercise the production packed writer with prebuilt payloads, without AliSim.
use std::fs;
use std::path::PathBuf;
use std::thread;
use std::time::Instant;

use clap::{Parser, ValueEnum};
use serde_json::json;

#[path = "../src/packed.rs"]
mod packed;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Format {
    Packed,
    Files,
}

#[derive(Parser)]
struct Args {
    #[arg(long)]
    out: PathBuf,
    #[arg(long, value_enum, default_value = "packed")]
    format: Format,
    #[arg(long, default_value_t = 2048)]
    samples: u64,
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    workers: u16,
    /// Approximate FASTA payload bytes; payload construction is outside timing.
    #[arg(long, default_value_t = 65_536)]
    payload_bytes: usize,
    #[arg(long, default_value_t = 1_073_741_824)]
    shard_size_bytes: u64,
    #[arg(long, default_value_t = packed::DEFAULT_SHARD_BUFFER_BYTES, value_parser = clap::value_parser!(u32).range(1024..))]
    shard_buffer_bytes: u32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    if args.samples < u64::from(args.workers) || args.payload_bytes < 16 {
        return Err("samples must cover all workers and payload must be at least 16 bytes".into());
    }
    fs::create_dir(&args.out)?;
    // Equal payloads across formats/variants isolate output handling, not biology
    // or FASTA serialization. Threads own separate SQLite connections and shards.
    let mut sequence = Vec::from(&b">T1\n"[..]);
    sequence.resize(args.payload_bytes - 1, b'A');
    sequence.push(b'\n');
    let tree = b"(T1:0.1,T2:0.1);\n";
    let started = Instant::now();
    match args.format {
        Format::Packed => fs::create_dir(args.out.join("shards"))?,
        Format::Files => {
            for name in ["trees", "msa_untrimmed", "metadata"] {
                fs::create_dir(args.out.join(name))?;
            }
        }
    }
    let timings = thread::scope(|scope| {
        let handles = (0..args.workers)
            .map(|worker| {
                let args = &args;
                let sequence = &sequence;
                scope.spawn(move || -> std::io::Result<serde_json::Value> {
                    let initialized = Instant::now();
                    let mut writer = match args.format {
                        Format::Packed => Some(packed::PackedWriter::new(
                            args.out.join("shards"),
                            args.shard_size_bytes,
                            worker,
                            String::new(),
                            args.shard_buffer_bytes,
                        )?),
                        Format::Files => None,
                    };
                    let init_s = initialized.elapsed().as_secs_f64();
                    let writing = Instant::now();
                    let mut count = 0;
                    for candidate in
                        (u64::from(worker)..args.samples).step_by(usize::from(args.workers))
                    {
                        let name = format!("sample_{candidate:06}");
                        let metadata = json!({"sample": name, "candidate_id": candidate,
                            "benchmark": true, "padding": "x".repeat(1500)});
                        if let Some(writer) = &mut writer {
                            writer.write(candidate, &name, tree, sequence, &metadata)?;
                        } else {
                            fs::write(args.out.join("trees").join(format!("{name}.nwk")), tree)?;
                            fs::write(
                                args.out.join("msa_untrimmed").join(format!("{name}.fa")),
                                sequence,
                            )?;
                            fs::write(
                                args.out.join("metadata").join(format!("{name}.json")),
                                serde_json::to_vec_pretty(&metadata)?,
                            )?;
                        }
                        count += 1;
                    }
                    let write_s = writing.elapsed().as_secs_f64();
                    let finishing = Instant::now();
                    if let Some(writer) = writer {
                        writer.finish()?;
                    }
                    Ok(json!({"worker": worker, "samples": count, "init_s": init_s,
                        "write_s": write_s, "finish_s": finishing.elapsed().as_secs_f64()}))
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .map_err(|_| std::io::Error::other("writer thread panicked"))?
            })
            .collect::<std::io::Result<Vec<_>>>()
    })?;
    let merging = Instant::now();
    if matches!(args.format, Format::Packed) {
        packed::merge_workers(&args.out, u64::from(args.workers), args.samples)?;
    }
    let merge_s = merging.elapsed().as_secs_f64();
    let total_s = started.elapsed().as_secs_f64();
    // Validation is outside the measured interval, including database reads.
    if matches!(args.format, Format::Packed) {
        assert!(!args.out.join(packed::PARTIAL_INDEX).exists());
        let db = rusqlite::Connection::open(args.out.join(packed::INDEX))?;
        let count: i64 = db.query_row("SELECT COUNT(*) FROM pairs", [], |row| row.get(0))?;
        assert_eq!(count, i64::try_from(args.samples)?);
        let integrity: String = db.query_row("PRAGMA integrity_check", [], |row| row.get(0))?;
        assert_eq!(integrity, "ok");
    }
    println!(
        "{}",
        json!({"format": format!("{:?}", args.format).to_lowercase(),
        "samples": args.samples, "workers": args.workers, "payload_bytes": sequence.len(),
        "shard_size_bytes": args.shard_size_bytes, "shard_buffer_bytes": args.shard_buffer_bytes, "total_s": total_s,
        "merge_s": merge_s, "worker_timings": timings})
    );
    Ok(())
}
