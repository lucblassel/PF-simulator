//! Emit empirical alpha draws for prior validation or individual-attempt replay.
use clap::Parser;
use pf_simulator::site_rates::sample_alpha;

#[derive(Parser)]
struct Args {
    #[arg(long, default_value_t = 42)]
    seed: u64,
    #[arg(long, default_value_t = 20000)]
    count: u64,
    #[arg(long, default_value_t = 0)]
    candidate: u64,
    #[arg(long, default_value_t = 0)]
    attempt: u32,
}

fn main() {
    let args = Args::parse();
    for offset in 0..args.count {
        println!(
            "{}",
            sample_alpha(args.seed, args.candidate + offset, args.attempt)
        );
    }
}
