/// SplitMix64-derived deterministic seeds, independent across samples and MSA retries.
pub fn sample_seed(seed: u64, sample: u64) -> u64 {
    let mut value = seed.wrapping_add(0x9e3779b97f4a7c15_u64.wrapping_mul(sample.wrapping_add(1)));
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_vectors() {
        assert_eq!(sample_seed(0, 0), 0xe220a8397b1dcdaf);
        assert_eq!(sample_seed(0, 1), 0x6e789e6aa1b965f4);
    }
}
