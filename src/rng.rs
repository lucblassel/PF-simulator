/// Small, deterministic xoshiro256** generator seeded with SplitMix64.
/// Used for topology and branch lengths; AliSim has its own RNG.
pub(crate) struct Rng([u64; 4]);

pub fn sample_seed(seed: u64, sample: u64) -> u64 {
    let mut value = seed.wrapping_add(0x9e3779b97f4a7c15_u64.wrapping_mul(sample.wrapping_add(1)));
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(std::array::from_fn(|i| sample_seed(seed, i as u64)))
    }

    fn next(&mut self) -> u64 {
        let result = self.0[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let temporary = self.0[1] << 17;
        self.0[2] ^= self.0[0];
        self.0[3] ^= self.0[1];
        self.0[1] ^= self.0[2];
        self.0[0] ^= self.0[3];
        self.0[2] ^= temporary;
        self.0[3] = self.0[3].rotate_left(45);
        result
    }

    pub fn index(&mut self, upper: usize) -> usize {
        let upper = upper as u64;
        let threshold = upper.wrapping_neg() % upper;
        loop {
            let value = self.next();
            if value >= threshold {
                return (value % upper) as usize;
            }
        }
    }

    pub fn exponential(&mut self, mean: f64) -> f64 {
        // 52 bits plus an exactly representable midpoint gives 0 < U < 1.
        let uniform = ((self.next() >> 12) as f64 + 0.5) / 4_503_599_627_370_496.0;
        -mean * uniform.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_vectors() {
        assert_eq!(sample_seed(0, 0), 0xe220a8397b1dcdaf);
        assert_eq!(sample_seed(0, 1), 0x6e789e6aa1b965f4);
        let mut rng = Rng([1, 2, 3, 4]);
        assert_eq!(rng.next(), 11520);
        assert_eq!(rng.next(), 0);
        assert_eq!(rng.next(), 1509978240);
    }
}
