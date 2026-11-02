//! Worker-owned, uncompressed tar shards with a packed_data.py-compatible index.
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use rusqlite::{Connection, OpenFlags, params};
use serde_json::Value;

pub(crate) const INDEX: &str = "index.sqlite3";
pub(crate) const PARTIAL_INDEX: &str = "index.sqlite3.partial";
const SCHEMA: &str = include_str!("packed_schema.sql");

fn integer(value: u64) -> io::Result<i64> {
    i64::try_from(value).map_err(io::Error::other)
}

fn worker_index_name(worker: u64) -> String {
    format!("worker-{worker:03}.index.sqlite3")
}

fn database(index: &Path) -> io::Result<Connection> {
    let path = index.with_extension("sqlite3.partial");
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)?;
    let db = Connection::open(path).map_err(io::Error::other)?;
    db.execute_batch(SCHEMA).map_err(io::Error::other)?;
    Ok(db)
}

fn finish_database(db: Connection, index: &Path) -> io::Result<()> {
    db.execute_batch("COMMIT;").map_err(io::Error::other)?;
    db.close().map_err(|(_, error)| io::Error::other(error))?;
    fs::rename(index.with_extension("sqlite3.partial"), index)
}

struct CountingWriter {
    inner: BufWriter<File>,
    bytes: u64,
}
impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let count = self.inner.write(bytes)?;
        self.bytes += count as u64;
        Ok(count)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

struct Shard {
    name: String,
    archive: tar::Builder<CountingWriter>,
}
impl Shard {
    fn append(&mut self, name: &str, bytes: &[u8]) -> io::Result<u64> {
        // Our generated member names fit a single ustar header: no hidden PAX
        // headers may intervene between this position and the payload.
        let mut header = tar::Header::new_ustar();
        header.set_path(name)?;
        header.set_size(bytes.len() as u64);
        header.set_mode(0o644);
        header.set_mtime(0);
        header.set_cksum();
        let offset = self.archive.get_ref().bytes + 512;
        self.archive.append(&header, bytes)?;
        Ok(offset)
    }

    fn finish(self, root: &Path) -> io::Result<()> {
        let mut writer = self.archive.into_inner()?;
        writer.flush()?;
        drop(writer);
        fs::rename(
            root.join(format!("{}.partial", self.name)),
            root.join(self.name),
        )
    }
}

pub(crate) struct PackedWriter {
    root: PathBuf,
    index: PathBuf,
    worker: u16,
    prefix: String,
    db: Connection,
    shard: Option<Shard>,
    next_shard: u64,
    count: u64,
    target: u64,
}
impl PackedWriter {
    pub(crate) fn new(root: PathBuf, target: u64, worker: u16, prefix: String) -> io::Result<Self> {
        let index = root.join(worker_index_name(u64::from(worker)));
        let db = database(&index)?;
        Ok(Self {
            root,
            index,
            worker,
            prefix,
            db,
            shard: None,
            next_shard: 0,
            count: 0,
            target,
        })
    }

    pub(crate) fn write(
        &mut self,
        candidate: u64,
        name: &str,
        tree: &[u8],
        sequence: &[u8],
        metadata: &Value,
    ) -> io::Result<()> {
        let footprint = 1024
            + (tree.len() as u64).div_ceil(512) * 512
            + (sequence.len() as u64).div_ceil(512) * 512;
        if self
            .shard
            .as_ref()
            .is_some_and(|s| s.archive.get_ref().bytes + footprint + 1024 > self.target)
        {
            self.shard.take().expect("open shard").finish(&self.root)?;
        }
        if self.shard.is_none() {
            let name = format!(
                "{}worker-{:03}-shard-{:06}.tar",
                self.prefix, self.worker, self.next_shard
            );
            self.next_shard += 1;
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(self.root.join(format!("{name}.partial")))?;
            self.shard = Some(Shard {
                name,
                archive: tar::Builder::new(CountingWriter {
                    inner: BufWriter::new(file),
                    bytes: 0,
                }),
            });
        }
        let shard = self.shard.as_mut().expect("created shard");
        // Physical member names include the writer identity, while index names
        // remain worker-independent logical sample IDs for reproducibility.
        let member = format!(
            "{}worker-{:03}-sample_{candidate:06}",
            self.prefix, self.worker
        );
        let tree_name = format!("trees/{member}.nwk");
        let sequence_name = format!("sequences/{member}.fa");
        let tree_offset = shard.append(&tree_name, tree)?;
        let sequence_offset = shard.append(&sequence_name, sequence)?;
        self.db
            .prepare_cached("INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")
            .map_err(io::Error::other)?
            .execute(params![
                integer(self.count)?,
                name,
                shard.name,
                tree_name,
                integer(tree_offset)?,
                integer(tree.len() as u64)?,
                sequence_name,
                integer(sequence_offset)?,
                integer(sequence.len() as u64)?
            ])
            .map_err(io::Error::other)?;
        self.db
            .prepare_cached("INSERT INTO sample_metadata VALUES (?, ?, ?)")
            .map_err(io::Error::other)?
            .execute(params![name, integer(candidate)?, metadata.to_string()])
            .map_err(io::Error::other)?;
        self.count += 1;
        if self.count.is_multiple_of(1000) {
            self.db
                .execute_batch("COMMIT; BEGIN IMMEDIATE;")
                .map_err(io::Error::other)?;
        }
        Ok(())
    }

    pub(crate) fn finish(mut self) -> io::Result<()> {
        if let Some(shard) = self.shard.take() {
            shard.finish(&self.root)?;
        }
        finish_database(self.db, &self.index)
    }
}

/// Merge only index rows; worker-prefixed tar files stay in the shared shards directory.
pub(crate) fn merge_workers(root: &Path, workers: u64, expected: u64) -> io::Result<()> {
    let db = database(&root.join(INDEX))?;
    db.execute_batch("CREATE TEMP TABLE incoming AS SELECT pairs.*, sample_metadata.candidate_id, sample_metadata.metadata_json FROM pairs JOIN sample_metadata USING(name);")
        .map_err(io::Error::other)?;
    for worker in 0..workers {
        let source = Connection::open_with_flags(
            root.join("shards").join(worker_index_name(worker)),
            OpenFlags::SQLITE_OPEN_READ_ONLY,
        )
        .map_err(io::Error::other)?;
        let mut statement = source.prepare("SELECT p.*, m.candidate_id, m.metadata_json FROM pairs p JOIN sample_metadata m USING(name) ORDER BY p.id")
            .map_err(io::Error::other)?;
        let mut rows = statement.query([]).map_err(io::Error::other)?;
        while let Some(row) = rows.next().map_err(io::Error::other)? {
            let shard: String = row.get(2).map_err(io::Error::other)?;
            let mut values = (0..11)
                .map(|column| row.get::<_, rusqlite::types::Value>(column))
                .collect::<rusqlite::Result<Vec<_>>>()
                .map_err(io::Error::other)?;
            values[2] = format!("shards/{shard}").into();
            db.prepare_cached("INSERT INTO incoming VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
                .map_err(io::Error::other)?
                .execute(rusqlite::params_from_iter(values))
                .map_err(io::Error::other)?;
        }
    }
    db.execute_batch("INSERT INTO pairs SELECT ROW_NUMBER() OVER (ORDER BY candidate_id) - 1, name, shard, tree_name, tree_offset, tree_size, sequence_name, sequence_offset, sequence_size FROM incoming;
        INSERT INTO sample_metadata SELECT name, candidate_id, metadata_json FROM incoming;
        DROP TABLE incoming;") .map_err(io::Error::other)?;
    let count: i64 = db
        .query_row("SELECT COUNT(*) FROM pairs", [], |row| row.get(0))
        .map_err(io::Error::other)?;
    if count != integer(expected)? {
        return Err(io::Error::other(format!(
            "Packed index has {count} pairs, expected {expected}"
        )));
    }
    finish_database(db, &root.join(INDEX))?;
    for worker in 0..workers {
        fs::remove_file(root.join("shards").join(worker_index_name(worker)))?;
    }
    Ok(())
}
