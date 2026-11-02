//! Exact tip-only selection on numeric rows, before FASTA formatting.
use std::collections::{HashMap, HashSet};

use clap::ValueEnum;

use crate::alignment::RawAlignment;

/// Which sliding-window endpoints to consider.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum WindowPolicy {
    /// Include every valid start, including the final window.
    #[default]
    Inclusive,
    /// Reproduce Python's excluded final start when M > L; M == L uses zero.
    Python,
}

/// Independently enabled acceptance checks; callers apply the allow-duplicates override.
#[derive(Clone, Copy, Debug)]
pub struct TipChecks {
    pub gapped: bool,
    pub ungapped: bool,
    pub window: bool,
    pub window_policy: WindowPolicy,
}

/// Why an attempt failed the first enabled check (checks are evaluated in this order).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rejection {
    Gapped,
    Ungapped,
    Window,
}

impl Rejection {
    pub fn name(self) -> &'static str {
        match self {
            Self::Gapped => "duplicate_gapped_tips",
            Self::Ungapped => "duplicate_ungapped_tips",
            Self::Window => "no_unique_tip_window",
        }
    }
}

/// A zero-based, half-open interval in the full alignment's column coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Window {
    pub start: usize,
    pub end: usize,
    pub checked: bool,
}

/// Select a window, or a prefix when the window check is disabled.
///
/// Tip membership uses explicit T1..Tn identities from the simulator's tree ABI,
/// never internal rows or a broad prefix match. Structural errors remain fatal.
/// ```
/// use pf_simulator::{alignment::RawAlignment, selection::{select, TipChecks, WindowPolicy}};
/// let raw = RawAlignment { names: vec!["T1".into(), "T2".into()],
///     rows: vec![vec![0, 1], vec![0, 2]], alphabet: b"ACG-".to_vec(), ..Default::default() };
/// let checks = TipChecks { gapped: true, ungapped: false, window: true,
///     window_policy: WindowPolicy::Inclusive };
/// let window = select(&raw, 2, 1, checks)?.unwrap();
/// assert_eq!(window.start, 1);
/// # Ok::<(), String>(())
/// ```
pub fn select(
    raw: &RawAlignment,
    tips: usize,
    length: usize,
    checks: TipChecks,
) -> Result<Result<Window, Rejection>, String> {
    if raw.names.len() != raw.rows.len() || length == 0 {
        return Err("Invalid numeric alignment dimensions or window length".into());
    }
    let columns = raw.rows.first().ok_or("Empty numeric alignment")?.len();
    if raw.rows.iter().any(|row| row.len() != columns) {
        return Err("Numeric alignment has unequal row lengths".into());
    }
    if tips == 0
        || raw
            .rows
            .iter()
            .flatten()
            .any(|&code| usize::from(code) >= raw.alphabet.len())
    {
        return Err("Invalid tip count or numeric state outside alphabet".into());
    }
    let mut indices = HashMap::with_capacity(raw.names.len());
    for (i, name) in raw.names.iter().enumerate() {
        if indices.insert(name.as_str(), i).is_some() {
            return Err(format!("Repeated sequence name: {name}"));
        }
    }
    let rows: Vec<&[u8]> = (1..=tips)
        .map(|tip| {
            indices
                .get(format!("T{tip}").as_str())
                .map(|&i| raw.rows[i].as_slice())
                .ok_or_else(|| format!("Missing tip T{tip}"))
        })
        .collect::<Result<_, _>>()?;
    if checks.gapped && duplicates(&rows) {
        return Ok(Err(Rejection::Gapped));
    }
    if checks.ungapped {
        let mut seen = HashSet::with_capacity(tips);
        for row in &rows {
            let ungapped: Vec<u8> = row
                .iter()
                .copied()
                .filter(|&code| raw.alphabet[usize::from(code)] != b'-')
                .collect();
            if !seen.insert(ungapped) {
                return Ok(Err(Rejection::Ungapped));
            }
        }
    }
    if !checks.window {
        return Ok(Ok(Window {
            start: 0,
            end: columns.min(length),
            checked: false,
        }));
    }
    if columns < length {
        return Ok(Err(Rejection::Window));
    }
    let starts = if checks.window_policy == WindowPolicy::Python && columns > length {
        columns - length
    } else {
        columns - length + 1
    };
    match unique_window(&rows, length, starts) {
        Some(start) => Ok(Ok(Window {
            start,
            end: start + length,
            checked: true,
        })),
        None => Ok(Err(Rejection::Window)),
    }
}

fn duplicates(rows: &[&[u8]]) -> bool {
    let mut seen = HashSet::with_capacity(rows.len());
    rows.iter().any(|row| !seen.insert(*row))
}

// Rolling hashes only group candidates: exact byte comparisons resolve collisions.
// O(tips * (L + number_of_starts)) hashing, without serializing numeric rows.
fn unique_window(rows: &[&[u8]], length: usize, starts: usize) -> Option<usize> {
    let base = 257_u64;
    let power = base.wrapping_pow((length - 1) as u32);
    let mut hashes: Vec<u64> = rows
        .iter()
        .map(|row| {
            row[..length].iter().fold(0_u64, |hash, &code| {
                hash.wrapping_mul(base).wrapping_add(u64::from(code) + 1)
            })
        })
        .collect();
    let mut heads = HashMap::with_capacity(rows.len());
    let mut previous = vec![None; rows.len()];
    for start in 0..starts {
        heads.clear();
        let mut unique = true;
        for (i, &hash) in hashes.iter().enumerate() {
            let head = heads.insert(hash, i);
            previous[i] = head;
            let mut other = head;
            while let Some(j) = other {
                if rows[i][start..start + length] == rows[j][start..start + length] {
                    unique = false;
                    break;
                }
                other = previous[j];
            }
            if !unique {
                break;
            }
        }
        if unique {
            return Some(start);
        }
        if start + 1 < starts {
            for (hash, row) in hashes.iter_mut().zip(rows) {
                *hash = hash
                    .wrapping_sub((u64::from(row[start]) + 1).wrapping_mul(power))
                    .wrapping_mul(base)
                    .wrapping_add(u64::from(row[start + length]) + 1);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    fn raw(rows: &[(&str, &str)]) -> RawAlignment {
        RawAlignment {
            names: rows.iter().map(|r| r.0.into()).collect(),
            rows: rows.iter().map(|r| r.1.as_bytes().to_vec()).collect(),
            alphabet: (0..=255).collect(),
            export_unaligned: false,
        }
    }
    fn checks(gapped: bool, ungapped: bool, window: bool) -> TipChecks {
        TipChecks {
            gapped,
            ungapped,
            window,
            window_policy: WindowPolicy::Inclusive,
        }
    }
    #[test]
    fn gapped_and_ungapped_are_independent_and_ancestors_are_ignored() {
        let data = raw(&[("T1", "A-C"), ("T2", "AC-"), ("I1", "A-C")]);
        assert!(
            select(&data, 2, 2, checks(true, false, false))
                .unwrap()
                .is_ok()
        );
        assert_eq!(
            select(&data, 2, 2, checks(false, true, false)).unwrap(),
            Err(Rejection::Ungapped)
        );
        let data = raw(&[("T1", "AC"), ("T2", "AC")]);
        assert_eq!(
            select(&data, 2, 2, checks(true, false, false)).unwrap(),
            Err(Rejection::Gapped)
        );
    }
    #[test]
    fn window_edges_prefix_short_and_missing_tips() {
        let data = raw(&[("T1", "AAA"), ("T2", "AAC")]);
        assert_eq!(
            select(&data, 2, 2, checks(false, false, true))
                .unwrap()
                .unwrap()
                .start,
            1
        );
        let mut python = checks(false, false, true);
        python.window_policy = WindowPolicy::Python;
        assert_eq!(select(&data, 2, 2, python).unwrap(), Err(Rejection::Window));
        assert!(select(&data, 2, 3, python).unwrap().is_ok());
        assert_eq!(select(&data, 2, 4, python).unwrap(), Err(Rejection::Window));
        assert_eq!(
            select(&data, 2, 4, checks(false, false, false))
                .unwrap()
                .unwrap(),
            Window {
                start: 0,
                end: 3,
                checked: false
            }
        );
        assert!(select(&data, 3, 1, python).is_err());
    }
    #[test]
    fn rolling_windows_match_exhaustive_exact_search() {
        use rand::{Rng, SeedableRng, rngs::StdRng};
        let mut rng = StdRng::seed_from_u64(523);
        for _ in 0..300 {
            let tips = rng.gen_range(2..10);
            let columns = rng.gen_range(1..40);
            let data: Vec<Vec<u8>> = (0..tips)
                .map(|_| (0..columns).map(|_| rng.gen_range(0..4)).collect())
                .collect();
            let rows: Vec<&[u8]> = data.iter().map(Vec::as_slice).collect();
            for length in 1..=columns {
                let expected = (0..=columns - length).find(|&start| {
                    let windows: Vec<&[u8]> =
                        rows.iter().map(|row| &row[start..start + length]).collect();
                    !duplicates(&windows)
                });
                assert_eq!(unique_window(&rows, length, columns - length + 1), expected);
            }
        }
    }
}
