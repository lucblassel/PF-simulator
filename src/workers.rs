//! Process isolation keeps IQ-TREE's global state private to each worker.
use std::error::Error;
use std::ffi::OsString;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, mpsc};
use std::thread::{self, JoinHandle};

use indicatif::ProgressBar;

use crate::Config;
use pf_simulator::selection::Rejection;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Outcome {
    Written,
    Skipped(Rejection),
}

#[derive(Default)]
pub(crate) struct Counts {
    pub written: u64,
    pub failed: u64,
    pub failed_gapped: u64,
    pub failed_ungapped: u64,
    pub failed_window: u64,
}

impl Counts {
    pub fn record(&mut self, outcome: Outcome) {
        match outcome {
            Outcome::Written => self.written += 1,
            Outcome::Skipped(reason) => {
                self.failed += 1;
                match reason {
                    Rejection::Gapped => self.failed_gapped += 1,
                    Rejection::Ungapped => self.failed_ungapped += 1,
                    Rejection::Window => self.failed_window += 1,
                }
            }
        }
    }

    pub fn completed(&self) -> u64 {
        self.written + self.failed
    }

    pub fn finished(&self, config: &Config) -> bool {
        self.progress_position(config) >= config.samples
    }

    fn progress_position(&self, config: &Config) -> u64 {
        if config.ensure_samples {
            self.written
        } else {
            self.completed()
        }
    }

    pub fn update_progress(&self, config: &Config, progress: &ProgressBar) {
        progress.set_message(format!("{} written, {} failed", self.written, self.failed));
        progress.set_position(self.progress_position(config));
    }
}

struct Worker {
    child: Child,
    dispatch: Option<mpsc::Sender<u64>>,
    thread: Option<JoinHandle<()>>,
}

impl Drop for Worker {
    fn drop(&mut self) {
        self.dispatch.take();
        // Also runs on startup failures and early returns. Kill before joining
        // so a reader blocked on a failed worker cannot hang the coordinator.
        let _ = self.child.kill();
        let _ = self.child.wait();
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

struct Event {
    worker: usize,
    sample: u64,
    result: Result<Outcome, WorkerFailure>,
}

enum WorkerFailure {
    OutputClosed,
    Message(String),
}

impl From<String> for WorkerFailure {
    fn from(message: String) -> Self {
        Self::Message(message)
    }
}

fn describe_failure(child: &mut Child, failure: WorkerFailure) -> String {
    let (message, status) = match failure {
        // EOF is delivered during process teardown. Wait for the actual status
        // before cleanup can kill the child and obscure its original signal.
        WorkerFailure::OutputClosed => (
            "Worker exited without reporting its result".to_owned(),
            child.wait().map(Some),
        ),
        WorkerFailure::Message(message) => (message, child.try_wait()),
    };
    match status {
        Ok(Some(status)) => format!("{message}; worker process exited with {status}"),
        Ok(None) => message,
        Err(error) => format!("{message}; could not obtain worker exit status: {error}"),
    }
}

/// Workers inherit simulation options, but never append to the shared log.
fn worker_arguments(seed: u64) -> Vec<OsString> {
    let mut result = Vec::new();
    let mut arguments = std::env::args_os().skip(1);
    while let Some(argument) = arguments.next() {
        if argument == "--alisim-log" || argument == "--seed" {
            arguments.next();
        } else if !argument.as_encoded_bytes().starts_with(b"--alisim-log=")
            && !argument.as_encoded_bytes().starts_with(b"--seed=")
        {
            result.push(argument);
        }
    }
    result.push("--seed".into());
    result.push(seed.to_string().into());
    result.push("--worker".into());
    result
}

fn merge_logs(destination: Option<&Path>, logs: &[PathBuf]) -> std::io::Result<()> {
    if let Some(destination) = destination {
        let mut output = OpenOptions::new().append(true).open(destination)?;
        for (id, path) in logs.iter().enumerate() {
            writeln!(output, "\n=== Worker {id} ===")?;
            std::io::copy(&mut fs::File::open(path)?, &mut output)?;
        }
        output.flush()?;
        // Keep all fragments if merging fails, rather than losing diagnostics.
        for path in logs {
            fs::remove_file(path)?;
        }
    }
    Ok(())
}

pub(crate) fn run(
    config: &Config,
    progress: &ProgressBar,
    counts: &mut Counts,
) -> Result<(), Box<dyn Error>> {
    let jobs = u64::from(config.jobs).min(config.samples) as usize;
    let executable = std::env::current_exe()?;
    let arguments = worker_arguments(config.seed.expect("parent resolved master seed"));
    let mut next = 0_u64;
    let mut limit = config.samples;
    let cancelled = Arc::new(AtomicBool::new(false));
    let (send, receive) = mpsc::channel::<Event>();
    let mut workers = Vec::with_capacity(jobs);
    let mut logs = Vec::with_capacity(jobs);

    let result = (|| -> Result<(), Box<dyn Error>> {
        for id in 0..jobs {
            if cancelled.load(Ordering::Acquire) {
                break;
            }
            let mut command = Command::new(&executable);
            command
                .args(&arguments)
                .arg("--worker-id")
                .arg(id.to_string())
                .stdin(Stdio::piped())
                .stdout(Stdio::piped());
            // Native errors stay visible; normal AliSim output is redirected.
            command.stderr(Stdio::inherit());
            if let Some(path) = &config.alisim_log {
                let mut name = path.as_os_str().to_owned();
                name.push(format!(".worker-{id:03}.log"));
                let path = PathBuf::from(name);
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(&path)?;
                command.arg("--alisim-log").arg(&path);
                logs.push(path);
            }
            let child = command.spawn()?;
            workers.push(Worker {
                child,
                dispatch: None,
                thread: None,
            });
            let worker = workers.last_mut().unwrap();
            let mut input = worker.child.stdin.take().unwrap();
            let mut output = BufReader::new(worker.child.stdout.take().unwrap());
            let cancelled = Arc::clone(&cancelled);
            let send = send.clone();
            let (dispatch, tasks) = mpsc::channel::<u64>();
            worker.dispatch = Some(dispatch);
            worker.thread = Some(thread::Builder::new().name(format!("worker-{id}")).spawn(
                move || {
                    let mut line = String::new();
                    for sample in tasks {
                        if cancelled.load(Ordering::Acquire) {
                            break;
                        }
                        let result = (|| -> Result<Outcome, WorkerFailure> {
                            writeln!(input, "{sample}")
                                .and_then(|_| input.flush())
                                .map_err(|e| format!("Cannot dispatch sample: {e}"))?;
                            line.clear();
                            if output.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
                                return Err(WorkerFailure::OutputClosed);
                            }
                            match line.trim_end().split_once('\t') {
                                Some((reported, status)) if reported == sample.to_string() => {
                                    match status {
                                        "written" => Ok(Outcome::Written),
                                        "duplicate_gapped_tips" => {
                                            Ok(Outcome::Skipped(Rejection::Gapped))
                                        }
                                        "duplicate_ungapped_tips" => {
                                            Ok(Outcome::Skipped(Rejection::Ungapped))
                                        }
                                        "no_unique_tip_window" => {
                                            Ok(Outcome::Skipped(Rejection::Window))
                                        }
                                        _ => Err(format!("Invalid worker status: {status}").into()),
                                    }
                                }
                                _ => Err(format!("Invalid worker response: {line:?}").into()),
                            }
                        })();
                        if result.is_err() {
                            cancelled.store(true, Ordering::Release);
                        }
                        if send
                            .send(Event {
                                worker: id,
                                sample,
                                result,
                            })
                            .is_err()
                        {
                            break;
                        }
                    }
                    // Closing stdin tells the persistent worker to exit normally.
                },
            )?);
            worker.dispatch.as_ref().unwrap().send(next)?;
            next += 1; // At most min(jobs, requested samples) initial tasks.
        }
        drop(send);
        while !counts.finished(config) {
            let event = receive
                .recv()
                .map_err(|_| "Workers stopped before completing the batch")?;
            let outcome = event.result.map_err(|failure| {
                let error = describe_failure(&mut workers[event.worker].child, failure);
                format!(
                    "Worker {} failed on sample_{:06}: {error}",
                    event.worker, event.sample
                )
            })?;
            counts.record(outcome);
            counts.update_progress(config, progress);
            // Each exhausted candidate earns exactly one replacement. Therefore
            // written + in-flight + unissued slots never exceeds the target.
            // Consecutive candidate IDs also preserve worker-count determinism.
            if config.ensure_samples && matches!(outcome, Outcome::Skipped(_)) {
                limit = limit
                    .checked_add(1)
                    .ok_or("Candidate sample IDs exhausted")?;
            }
            if next < limit && !cancelled.load(Ordering::Acquire) {
                workers[event.worker]
                    .dispatch
                    .as_ref()
                    .unwrap()
                    .send(next)?;
                next += 1;
            }
        }
        for worker in &mut workers {
            worker.dispatch.take();
        }
        for (id, worker) in workers.iter_mut().enumerate() {
            if let Some(thread) = worker.thread.take() {
                thread
                    .join()
                    .map_err(|_| format!("Worker {id} reader panicked"))?;
            }
            let status = worker.child.wait()?;
            if !status.success() {
                return Err(format!("Worker {id} exited with {status}").into());
            }
        }
        Ok(())
    })();
    cancelled.store(true, Ordering::Release);
    // Stop every worker before joining any of them on fatal errors.
    if result.is_err() {
        for worker in &mut workers {
            worker.dispatch.take();
            let _ = worker.child.kill();
        }
    }
    drop(workers);
    let merged = merge_logs(config.alisim_log.as_deref(), &logs);
    result?;
    merged?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn closed_worker_output_reports_exit_code() {
        let mut child = Command::new("sh").args(["-c", "exit 7"]).spawn().unwrap();
        let message = describe_failure(&mut child, WorkerFailure::OutputClosed);
        assert!(message.contains("exit status: 7"), "{message}");
    }

    #[cfg(unix)]
    #[test]
    fn closed_worker_output_reports_signal_without_stderr() {
        let mut child = Command::new("sh")
            .args(["-c", "ulimit -c 0; kill -SEGV $$"])
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let message = describe_failure(&mut child, WorkerFailure::OutputClosed);
        assert!(message.contains("signal: 11"), "{message}");
    }

    #[test]
    fn only_successes_advance_progress_when_enforcing_output_count() {
        for ensure in [false, true] {
            let mut args = vec!["pf-sim", "--samples", "2"];
            if ensure {
                args.push("--ensure-samples");
            }
            let config = Config::try_parse_from(args).unwrap();
            let bar = ProgressBar::hidden();
            let mut counts = Counts::default();
            counts.record(Outcome::Skipped(Rejection::Gapped));
            counts.update_progress(&config, &bar);
            assert_eq!(bar.position(), if ensure { 0 } else { 1 });
            assert!(!counts.finished(&config));
            counts.record(Outcome::Written);
            counts.update_progress(&config, &bar);
            assert_eq!(bar.position(), if ensure { 1 } else { 2 });
            assert_eq!(counts.finished(&config), !ensure);
            if ensure {
                counts.record(Outcome::Written);
                counts.update_progress(&config, &bar);
                assert_eq!(bar.position(), 2);
                assert!(counts.finished(&config));
                assert_eq!(counts.completed(), 3);
            }
        }
    }
}
