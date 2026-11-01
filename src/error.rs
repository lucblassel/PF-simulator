//! Errors from the simulation workflow.
use std::{error::Error, fmt, io};

/// Duplicate exhaustion is recoverable at the sample level; other errors are fatal.
#[derive(Debug)]
pub enum SimulationError {
    /// Every allowed MSA attempt for this sample contained duplicate tip sequences.
    DuplicateSequences {
        /// Stable sample identifier, including skipped samples.
        sample: String,
        /// Number of attempts exhausted on the fixed tree.
        attempts: u32,
    },
    /// An underlying simulation or I/O error that must stop the run.
    Other(Box<dyn Error + Send + Sync>),
}

impl fmt::Display for SimulationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSequences { sample, attempts } => write!(
                formatter,
                "No duplicate-free MSA for {sample} after {attempts} attempts"
            ),
            Self::Other(source) => fmt::Display::fmt(source, formatter),
        }
    }
}

impl Error for SimulationError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::DuplicateSequences { .. } => None,
            Self::Other(source) => Some(source.as_ref()),
        }
    }
}

impl From<io::Error> for SimulationError {
    fn from(error: io::Error) -> Self {
        Self::Other(Box::new(error))
    }
}

impl From<String> for SimulationError {
    fn from(error: String) -> Self {
        // The native bridge currently returns error messages as strings.
        Self::Other(error.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrapper_preserves_the_underlying_io_error() {
        let error = SimulationError::from(io::Error::new(io::ErrorKind::PermissionDenied, "log"));
        let source = error.source().unwrap().downcast_ref::<io::Error>().unwrap();
        assert_eq!(source.kind(), io::ErrorKind::PermissionDenied);
    }
}
