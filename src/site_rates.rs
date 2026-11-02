//! Empirical Gamma shapes from the Python alignment pipeline.
//!
//! Values are bundled at compile time; Python and the original pickle are not
//! needed at runtime. Reproducibility uses a separate RNG stream per attempt.
use std::fmt;
use std::str::FromStr;
use std::sync::OnceLock;

use rand::{Rng, SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, Normal};

use crate::rng::sample_seed;

/// Embedded empirical observations and their SHA-256 provenance headers.
pub const EMPIRICAL_GAMMA: &str = include_str!("../data/empirical_gamma.tsv");

/// Empirical site-rate mode, independent of branch-length Gamma priors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GammaMode {
    /// Keep the supplied model unchanged.
    None,
    /// Continuous Gamma site rates.
    Continuous,
    /// Discrete Gamma with a positive number of categories.
    Discrete(u32),
}

impl FromStr for GammaMode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_uppercase().as_str() {
            "NONE" => Ok(Self::None),
            "GC" => Ok(Self::Continuous),
            text => text
                .strip_prefix('G')
                .and_then(|n| n.parse::<u32>().ok())
                .filter(|&n| n > 0 && n <= i32::MAX as u32)
                .map(Self::Discrete)
                .ok_or_else(|| "expected none, GC, or G[n] with n positive (e.g. G8)".into()),
        }
    }
}

impl fmt::Display for GammaMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Continuous => write!(f, "GC"),
            Self::Discrete(n) => write!(f, "G{n}"),
        }
    }
}

fn values() -> &'static [f64] {
    static VALUES: OnceLock<Vec<f64>> = OnceLock::new();
    VALUES.get_or_init(|| {
        EMPIRICAL_GAMMA
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| line.parse().expect("bundled positive alpha"))
            .collect()
    })
}

/// Draw a shape for a zero-based candidate and attempt using a dedicated stream.
///
/// Matches Python's law: select an observation uniformly (including repeated
/// observations), draw Normal(mean, mean/10), then floor at 0.05. RNG algorithms
/// differ from Python/NumPy; equal numeric seeds do not imply equal draws there.
///
/// ```
/// use pf_simulator::site_rates::sample_alpha;
/// let alpha = sample_alpha(42, 3, 0);
/// assert!(alpha >= 0.05);
/// assert_eq!(alpha, sample_alpha(42, 3, 0));
/// ```
pub fn sample_alpha(master_seed: u64, candidate: u64, attempt: u32) -> f64 {
    // A domain distinct from the tree and AliSim streams, versioned in summaries.
    let seed = sample_seed(
        sample_seed(master_seed, candidate),
        0x47414d4d41000000 ^ u64::from(attempt),
    );
    let mut rng = StdRng::seed_from_u64(seed);
    let values = values();
    let mean = values[rng.gen_range(0..values.len())];
    Normal::new(mean, mean / 10.0)
        .expect("positive empirical mean")
        .sample(&mut rng)
        .max(0.05)
}

impl GammaMode {
    /// Resolve a model for an attempt; `none` preserves it verbatim.
    ///
    /// ```
    /// use pf_simulator::site_rates::GammaMode;
    /// let model = GammaMode::Discrete(8).model("LG", 42, 0, 0);
    /// assert!(model.starts_with("LG+G8{"));
    /// ```
    pub fn model(self, base: &str, seed: u64, candidate: u64, attempt: u32) -> String {
        if self == Self::None {
            base.to_owned()
        } else {
            format!(
                "{base}+{self}{{{}}}",
                sample_alpha(seed, candidate, attempt)
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundle_and_modes() {
        assert_eq!(values().len(), 12408);
        assert!(values().iter().all(|v| v.is_finite() && *v > 0.0));
        assert_eq!(values().iter().copied().fold(f64::INFINITY, f64::min), 0.02);
        for mode in ["none", "GC", "G1", "G8", "G32"] {
            assert_eq!(mode.parse::<GammaMode>().unwrap().to_string(), mode);
        }
        for mode in ["G0", "G-1", "G", "G8{0.5}", "R4"] {
            assert!(mode.parse::<GammaMode>().is_err());
        }
    }

    #[test]
    fn attempts_have_independent_reproducible_draws() {
        let a = sample_alpha(123, 0, 0);
        assert_eq!(a, sample_alpha(123, 0, 0));
        assert_eq!(a, 0.7767144924613649); // Version-1 sampling reference vector.
        assert_ne!(a, sample_alpha(123, 0, 1));
        assert_ne!(a, sample_alpha(123, 1, 0));
        assert_ne!(a, sample_alpha(124, 0, 0));
    }
}
