use std::collections::HashSet;
use std::ffi::{CString, c_char, c_int, c_void};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::sync::Mutex;

use crate::tree::Tree;

#[repr(C)]
struct TreeEdge {
    a: u32,
    b: u32,
    length: f64,
}

enum TreeInput<'a> {
    File(&'a Path),
    Structure(&'a Tree),
}

static ALISIM_LOCK: Mutex<()> = Mutex::new(());

unsafe extern "C" {
    fn pf_alisim_simulate(
        tree_file: *const c_char,
        edges: *const TreeEdge,
        edge_count: usize,
        tips: u32,
        output_prefix: *const c_char,
        model: *const c_char,
        sequence_type: *const c_char,
        length: i32,
        seed: i32,
        indel: *const c_char,
        indel_size: *const c_char,
        no_unaligned: c_int,
        ancestral: c_int,
        log_path: *const c_char,
        output_mode: c_int,
        callback: unsafe extern "C" fn(*const u8, usize, c_int, *mut c_void) -> c_int,
        context: *mut c_void,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
}

#[derive(Clone, Copy, Debug)]
pub enum SequenceType {
    AminoAcid,
    Dna,
}

impl SequenceType {
    fn as_str(self) -> &'static str {
        match self {
            Self::AminoAcid => "AA",
            Self::Dna => "DNA",
        }
    }
}

#[derive(Clone, Debug)]
pub struct AlignmentConfig {
    pub length: u32,
    pub model: String,
    pub sequence_type: SequenceType,
    /// AliSim's insertion,deletion rate string, passed through unchanged.
    pub indel: Option<String>,
    /// AliSim's insertion/deletion size distribution string, passed through unchanged.
    pub indel_size: Option<String>,
    /// Suppress AliSim's additional `.unaligned.fa` output.
    pub no_unaligned: bool,
    /// Include simulated ancestral sequences in the shared alignment coordinates.
    pub ancestral: bool,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            length: 500,
            model: "LG".into(),
            sequence_type: SequenceType::AminoAcid,
            indel: None,
            indel_size: None,
            no_unaligned: false,
            ancestral: false,
        }
    }
}

/// Simulate tip sequences through linked AliSim. Calls from all Rust threads
/// are serialized. Some upstream IQ-TREE errors can terminate the process.
pub fn simulate_alignment(
    tree_file: &Path,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
) -> Result<(), String> {
    simulate_alignment_with_log(tree_file, output_prefix, config, seed, None)
}

/// Simulate quietly, optionally appending AliSim's informational output to a
/// file. Native errors remain on stderr. Calls share the same process-wide lock.
pub fn simulate_alignment_with_log(
    tree_file: &Path,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
    log_path: Option<&Path>,
) -> Result<(), String> {
    let data =
        simulate_alignment_numeric(tree_file, output_prefix, config, seed, log_path)?.into_fasta();
    let path = |suffix: &str| {
        let mut path = output_prefix.as_os_str().to_owned();
        path.push(suffix);
        std::path::PathBuf::from(path)
    };
    std::fs::write(path(".fa"), &data.aligned_fasta).map_err(|error| error.to_string())?;
    if let Some(unaligned) = data.unaligned_fasta {
        std::fs::write(path(".unaligned.fa"), unaligned).map_err(|error| error.to_string())?;
    }
    Ok(())
}

/// FASTA bytes owned by Rust. No alignment or indel temporary files are written
/// by the in-memory backend; the input Newick file is still read by AliSim.
#[derive(Default, Debug)]
pub struct AlignmentOutput {
    pub aligned_fasta: Vec<u8>,
    pub unaligned_fasta: Option<Vec<u8>>,
}

#[derive(Default)]
struct OutputCollector {
    fasta: AlignmentOutput,
    raw: Option<RawAlignment>,
}

/// Rust-owned numeric rows. Duplicate checks avoid FASTA formatting and parsing.
#[derive(Default, Debug)]
pub struct RawAlignment {
    pub names: Vec<String>,
    pub rows: Vec<Vec<u8>>,
    pub alphabet: Vec<u8>,
    pub export_unaligned: bool,
}

impl RawAlignment {
    /// Format an accepted attempt once, including optional gap-free output.
    pub fn into_fasta(self) -> AlignmentOutput {
        AlignmentOutput {
            aligned_fasta: self.to_fasta(false),
            unaligned_fasta: self.export_unaligned.then(|| self.to_fasta(true)),
        }
    }

    pub fn has_duplicates(&self) -> bool {
        let mut seen = HashSet::new();
        self.rows.iter().any(|row| !seen.insert(row.as_slice()))
    }

    /// Format a column interval without changing full-alignment coordinates.
    pub fn window_fasta(&self, start: usize, end: usize) -> Vec<u8> {
        self.format_fasta(false, Some(start..end))
    }

    pub fn to_fasta(&self, gap_free: bool) -> Vec<u8> {
        self.format_fasta(gap_free, None)
    }

    fn format_fasta(&self, gap_free: bool, window: Option<std::ops::Range<usize>>) -> Vec<u8> {
        let mut result = Vec::new();
        for (name, row) in self.names.iter().zip(&self.rows) {
            result.push(b'>');
            result.extend_from_slice(name.as_bytes());
            result.push(b'\n');
            let row = match &window {
                Some(range) => &row[range.clone()],
                None => row.as_slice(),
            };
            result.extend(
                row.iter()
                    .map(|&code| self.alphabet[code as usize])
                    .filter(|&symbol| !gap_free || symbol != b'-'),
            );
            result.push(b'\n');
        }
        result
    }
}

unsafe extern "C" fn collect_output(
    data: *const u8,
    size: usize,
    kind: c_int,
    context: *mut c_void,
) -> c_int {
    if (data.is_null() && size != 0) || context.is_null() || size > isize::MAX as usize {
        return 1;
    }
    // SAFETY: the bridge invokes this synchronously with a live byte buffer and
    // the AlignmentOutput passed below, with exclusive access for this call.
    let output = unsafe { &mut *context.cast::<OutputCollector>() };
    let mut bytes = Vec::new();
    if bytes.try_reserve_exact(size).is_err() {
        return 1;
    }
    // The allocation is already reserved, so copying cannot unwind into C++.
    if size != 0 {
        bytes.extend_from_slice(unsafe { std::slice::from_raw_parts(data, size) });
    }
    if let Some(raw) = &mut output.raw {
        match kind {
            0 => raw.rows.push(bytes),
            1 => match String::from_utf8(bytes) {
                Ok(name) => raw.names.push(name),
                Err(_) => return 1,
            },
            2 => raw.alphabet = bytes,
            3 if bytes.len() == 1 => raw.export_unaligned = bytes[0] != 0,
            _ => return 1,
        }
        return 0;
    }
    match kind {
        0 => output.fasta.aligned_fasta = bytes,
        1 => output.fasta.unaligned_fasta = Some(bytes),
        _ => return 1,
    }
    0
}

/// Simulate one MSA into Rust-owned memory, without creating sequence files.
/// `output_prefix` identifies virtual AliSim streams; its directory need not exist.
pub fn simulate_alignment_in_memory(
    tree_file: &Path,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
    log_path: Option<&Path>,
) -> Result<AlignmentOutput, String> {
    simulate_output(
        TreeInput::File(tree_file),
        output_prefix,
        config,
        seed,
        log_path,
        0,
    )
    .map(|output| output.fasta)
}

/// Simulate numeric tip rows without FASTA or temporary alignment files.
pub fn simulate_alignment_numeric(
    tree_file: &Path,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
    log_path: Option<&Path>,
) -> Result<RawAlignment, String> {
    numeric_output(simulate_output(
        TreeInput::File(tree_file),
        output_prefix,
        config,
        seed,
        log_path,
        1,
    )?)
}

/// Simulate directly from Rust's edge list, without tree I/O or Newick conversion.
/// IQ-TREE owns its native nodes; the Rust tree is borrowed only for this call.
pub fn simulate_alignment_from_tree(
    tree: &Tree,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
    log_path: Option<&Path>,
) -> Result<RawAlignment, String> {
    numeric_output(simulate_output(
        TreeInput::Structure(tree),
        output_prefix,
        config,
        seed,
        log_path,
        1,
    )?)
}

fn numeric_output(data: OutputCollector) -> Result<RawAlignment, String> {
    let raw = data.raw.ok_or("No numeric alignment returned")?;
    if raw.rows.is_empty()
        || raw.names.len() != raw.rows.len()
        || raw.alphabet.is_empty()
        || raw.rows.iter().any(|row| {
            row.len() != raw.rows[0].len()
                || row
                    .iter()
                    .any(|&code| usize::from(code) >= raw.alphabet.len())
        })
        || raw.names.iter().collect::<HashSet<_>>().len() != raw.names.len()
    {
        return Err("AliSim returned an invalid numeric alignment".into());
    }
    Ok(raw)
}

/// Compatibility entry point for the transport benchmark.
#[cfg(feature = "raw-benchmark")]
pub fn simulate_alignment_raw_benchmark(
    tree_file: &Path,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
) -> Result<RawAlignment, String> {
    simulate_alignment_numeric(tree_file, output_prefix, config, seed, None)
}

fn simulate_output(
    tree_input: TreeInput<'_>,
    output_prefix: &Path,
    config: &AlignmentConfig,
    seed: i32,
    log_path: Option<&Path>,
    output_mode: c_int,
) -> Result<OutputCollector, String> {
    if config.length == 0 || config.length > i32::MAX as u32 || !(1..i32::MAX).contains(&seed) {
        return Err("Invalid alignment length or seed".into());
    }
    if config.model.is_empty() {
        return Err("Model must not be empty".into());
    }
    let (tree_file, edges, tips) = match tree_input {
        TreeInput::File(path) => {
            if !path.is_file() {
                return Err(format!("Tree does not exist: {}", path.display()));
            }
            (path, Vec::new(), 0)
        }
        TreeInput::Structure(tree) => {
            let edges = tree
                .edges()
                .iter()
                .map(|edge| TreeEdge {
                    a: edge.a as u32,
                    b: edge.b as u32,
                    length: edge.length,
                })
                .collect();
            (Path::new("<in-memory-tree>"), edges, tree.tips() as u32)
        }
    };

    let cstring = |bytes: &[u8]| {
        CString::new(bytes).map_err(|_| "FFI strings cannot contain NUL bytes".to_string())
    };
    let tree = cstring(tree_file.as_os_str().as_bytes())?;
    let output = cstring(output_prefix.as_os_str().as_bytes())?;
    let model = cstring(config.model.as_bytes())?;
    let sequence_type = cstring(config.sequence_type.as_str().as_bytes())?;
    let indel = config
        .indel
        .as_deref()
        .map(|value| cstring(value.as_bytes()))
        .transpose()?;
    let indel_size = config
        .indel_size
        .as_deref()
        .map(|value| cstring(value.as_bytes()))
        .transpose()?;
    let log = log_path
        .map(|path| cstring(path.as_os_str().as_bytes()))
        .transpose()?;
    let mut error = [0_u8; 4096];
    let mut data = OutputCollector::default();
    if output_mode == 1 {
        data.raw = Some(RawAlignment::default());
    }
    let _guard = ALISIM_LOCK
        .lock()
        .map_err(|_| "AliSim lock was poisoned".to_string())?;

    // SAFETY: all input strings live through the call and are NUL-terminated;
    // the C-compatible edge array also lives through the call and is read-only.
    // error is writable for the supplied capacity. The mutex excludes concurrent
    // use of IQ-TREE globals, and the C++ bridge catches foreign exceptions.
    let status = unsafe {
        pf_alisim_simulate(
            tree.as_ptr(),
            if edges.is_empty() {
                std::ptr::null()
            } else {
                edges.as_ptr()
            },
            edges.len(),
            tips,
            output.as_ptr(),
            model.as_ptr(),
            sequence_type.as_ptr(),
            config.length as i32,
            seed,
            indel
                .as_ref()
                .map_or(std::ptr::null(), |value| value.as_ptr()),
            indel_size
                .as_ref()
                .map_or(std::ptr::null(), |value| value.as_ptr()),
            c_int::from(config.no_unaligned),
            c_int::from(config.ancestral),
            log.as_ref().map_or(std::ptr::null(), |path| path.as_ptr()),
            output_mode,
            collect_output,
            (&mut data as *mut OutputCollector).cast(),
            error.as_mut_ptr().cast(),
            error.len(),
        )
    };
    if status == 0 {
        return Ok(data);
    }
    let end = error
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(error.len());
    Err(String::from_utf8_lossy(&error[..end]).into_owned())
}

/// Exact equality of complete aligned tip sequences, including gap positions.
/// FASTA IDs and line wrapping do not affect duplicate detection.
pub fn has_duplicate_sequences(fasta: &[u8]) -> Result<bool, String> {
    let mut seen = HashSet::new();
    let mut sequence = Vec::new();
    let mut in_record = false;
    for line in fasta.split(|&byte| byte == b'\n').map(<[u8]>::trim_ascii) {
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b">") {
            if in_record && !seen.insert(std::mem::take(&mut sequence)) {
                return Ok(true);
            }
            in_record = true;
        } else if in_record {
            sequence.extend_from_slice(line);
        } else {
            return Err("Expected a FASTA header before sequence data".into());
        }
    }
    if !in_record {
        return Err("AliSim returned no FASTA records".into());
    }
    Ok(!seen.insert(sequence))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn duplicate_detection_uses_full_sequences_and_ignores_wrapping() {
        assert!(has_duplicate_sequences(b">T1\r\nAC-\r\nGT\r\n>T2\nAC-GT\n").unwrap());
        assert!(!has_duplicate_sequences(b">T1\nAC-GT\n>T2\nA-CGT\n").unwrap());
        assert!(has_duplicate_sequences(b">T1\n\n>T2\n\n").unwrap());
        assert!(has_duplicate_sequences(b"ACGT\n").is_err());
        assert!(has_duplicate_sequences(b"").is_err());
    }

    #[test]
    fn rejects_invalid_inputs_before_ffi() {
        let config = AlignmentConfig {
            length: 0,
            ..Default::default()
        };
        assert!(simulate_alignment(Path::new("missing"), Path::new("out"), &config, 1).is_err());
        assert!(
            simulate_alignment(
                Path::new("missing"),
                Path::new("out"),
                &Default::default(),
                0
            )
            .is_err()
        );
    }
}
