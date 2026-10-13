//! Best-effort stderr reporting, separate from the server's binary responses.

use std::io::{self, IsTerminal, Write};
use std::sync::{Mutex, MutexGuard};
use std::time::{Duration, Instant};

pub(super) struct ServerProgress {
    inner: Mutex<Progress<io::Stderr>>,
}

impl ServerProgress {
    pub(super) fn new() -> Self {
        let stderr = io::stderr();
        let terminal = stderr.is_terminal();
        Self {
            inner: Mutex::new(Progress::new(stderr, terminal)),
        }
    }

    fn lock(&self) -> MutexGuard<'_, Progress<io::Stderr>> {
        self.inner.lock().unwrap_or_else(|error| error.into_inner())
    }

    pub(super) fn start(&self) {
        self.lock().start(Instant::now());
    }

    pub(super) fn completed(&self, sequences: usize, bases: u64) {
        // Logging failure must not disrupt a successfully written response.
        let _ = self.lock().completed(sequences, bases, Instant::now());
    }

    pub(super) fn message(&self, message: &str) {
        let _ = self.lock().message(message);
    }

    pub(super) fn finish(&self) {
        let _ = self.lock().finish(Instant::now());
    }

    pub(super) fn end_line(&self) {
        let _ = self.lock().end_line();
    }

    #[cfg(test)]
    pub(super) fn totals(&self) -> (u64, u64, u64) {
        let progress = self.lock();
        (progress.batches, progress.sequences, progress.bases)
    }
}

struct Progress<W> {
    writer: W,
    terminal: bool,
    started: Option<Instant>,
    last_report: Option<Instant>,
    batches: u64,
    sequences: u64,
    bases: u64,
    line_len: usize,
}

impl<W: Write> Progress<W> {
    fn new(writer: W, terminal: bool) -> Self {
        Self {
            writer,
            terminal,
            started: None,
            last_report: None,
            batches: 0,
            sequences: 0,
            bases: 0,
            line_len: 0,
        }
    }

    fn start(&mut self, now: Instant) {
        if self.started.is_none() {
            self.started = Some(now);
            self.last_report = Some(now);
        }
    }

    fn completed(&mut self, sequences: usize, bases: u64, now: Instant) -> io::Result<()> {
        self.batches += 1;
        self.sequences += sequences as u64;
        self.bases += bases;
        if self
            .last_report
            .is_some_and(|last| now.duration_since(last) >= Duration::from_secs(5))
        {
            self.last_report = Some(now);
            self.render("Processed", now, false)?;
        }
        Ok(())
    }

    fn render(&mut self, label: &str, now: Instant, finished: bool) -> io::Result<()> {
        let elapsed = self
            .started
            .map_or(0.0, |start| now.duration_since(start).as_secs_f64());
        let rate = if elapsed > 0.0 {
            self.bases as f64 / elapsed
        } else {
            0.0
        };
        let line = format!(
            "{label} {} {} | {} {} | {} | {}/s | elapsed {:.1}s",
            self.batches,
            if self.batches == 1 {
                "batch"
            } else {
                "batches"
            },
            self.sequences,
            if self.sequences == 1 {
                "sequence"
            } else {
                "sequences"
            },
            format_bases(self.bases as f64),
            format_bases(rate),
            elapsed,
        );
        if self.terminal {
            let padding = self.line_len.saturating_sub(line.len());
            write!(self.writer, "\r{line}{:padding$}", "")?;
            self.line_len = line.len();
            if finished {
                writeln!(self.writer)?;
                self.line_len = 0;
            }
        } else {
            writeln!(self.writer, "{line}")?;
        }
        self.writer.flush()
    }

    fn end_line(&mut self) -> io::Result<()> {
        if self.line_len > 0 {
            writeln!(self.writer)?;
            self.line_len = 0;
        }
        self.writer.flush()
    }

    fn message(&mut self, message: &str) -> io::Result<()> {
        self.end_line()?;
        writeln!(self.writer, "{message}")?;
        self.writer.flush()
    }

    fn finish(&mut self, now: Instant) -> io::Result<()> {
        self.render("Finished", now, true)
    }
}

fn format_bases(value: f64) -> String {
    if value >= 1_000_000_000.0 {
        format!("{:.1} Gb", value / 1_000_000_000.0)
    } else if value >= 1_000_000.0 {
        format!("{:.1} Mb", value / 1_000_000.0)
    } else if value >= 1_000.0 {
        format!("{:.1} kb", value / 1_000.0)
    } else {
        format!("{value:.0} bp")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn reports_cumulative_work_at_most_every_five_seconds() {
        let start = Instant::now();
        let mut progress = Progress::new(Vec::new(), false);
        progress.start(start);
        progress
            .completed(2, 1000, start + Duration::from_secs(4))
            .unwrap();
        assert!(progress.writer.is_empty());
        progress
            .completed(1, 1500, start + Duration::from_secs(5))
            .unwrap();
        progress
            .completed(0, 0, start + Duration::from_secs(9))
            .unwrap();
        progress
            .completed(3, 2500, start + Duration::from_secs(10))
            .unwrap();
        let output = String::from_utf8(progress.writer).unwrap();
        assert_eq!(output.lines().count(), 2);
        assert!(
            output.contains("Processed 2 batches | 3 sequences | 2.5 kb | 500 bp/s | elapsed 5.0s")
        );
        assert!(output
            .contains("Processed 4 batches | 6 sequences | 5.0 kb | 500 bp/s | elapsed 10.0s"));
        assert!(!output.contains('\r'));
        assert!(!output.contains('\x1b'));
    }

    #[test]
    fn terminal_updates_replace_one_line_and_finish_with_newline() {
        let start = Instant::now();
        let mut progress = Progress::new(Vec::new(), true);
        progress.start(start);
        progress
            .completed(1, 100, start + Duration::from_secs(5))
            .unwrap();
        progress
            .completed(1, 100, start + Duration::from_secs(10))
            .unwrap();
        assert!(!progress.writer.contains(&b'\n'));
        progress.finish(start + Duration::from_secs(11)).unwrap();
        let output = String::from_utf8(progress.writer).unwrap();
        assert_eq!(output.matches('\r').count(), 3);
        assert_eq!(output.matches('\n').count(), 1);
        assert!(output.ends_with('\n'));
        assert!(!output.contains('\x1b'));
    }

    #[test]
    fn shorter_replacement_clears_previous_characters() {
        let start = Instant::now();
        let mut progress = Progress::new(Vec::new(), true);
        progress.line_len = 200;
        progress.finish(start).unwrap();
        let output = String::from_utf8(progress.writer).unwrap();
        assert_eq!(output.len(), 202); // carriage return + padded line + newline
    }

    #[test]
    fn messages_end_an_active_progress_line_and_preserve_readiness() {
        let start = Instant::now();
        let mut progress = Progress::new(Vec::new(), true);
        progress.message("READY stdio").unwrap();
        progress.start(start);
        progress
            .completed(1, 100, start + Duration::from_secs(5))
            .unwrap();
        progress
            .message("Connection ended: truncated request")
            .unwrap();
        let output = String::from_utf8(progress.writer).unwrap();
        assert!(output.starts_with("READY stdio\n"));
        assert!(output.contains("\nConnection ended: truncated request\n"));
        assert_eq!(progress.line_len, 0);
    }

    #[test]
    fn short_and_empty_sessions_get_a_final_summary() {
        let start = Instant::now();
        let mut progress = Progress::new(Vec::new(), false);
        progress.start(start);
        progress
            .completed(2, 100, start + Duration::from_millis(20))
            .unwrap();
        progress.finish(start + Duration::from_secs(1)).unwrap();
        let output = String::from_utf8(progress.writer).unwrap();
        assert_eq!(output.lines().count(), 1);
        assert!(output.contains("Finished 1 batch | 2 sequences | 100 bp"));
        let mut empty = Progress::new(Vec::new(), false);
        empty.finish(start).unwrap();
        assert!(String::from_utf8(empty.writer)
            .unwrap()
            .contains("0 batches | 0 sequences | 0 bp"));
    }

    #[test]
    fn concurrent_connections_share_totals_and_serialized_updates() {
        let start = Instant::now();
        let mut state = Progress::new(Vec::new(), false);
        state.start(start);
        let progress = Arc::new(Mutex::new(state));
        let handles: Vec<_> = (0..8)
            .map(|_| {
                let progress = Arc::clone(&progress);
                std::thread::spawn(move || {
                    for _ in 0..100 {
                        progress
                            .lock()
                            .unwrap()
                            .completed(2, 1000, start + Duration::from_secs(5))
                            .unwrap();
                    }
                })
            })
            .collect();
        for handle in handles {
            handle.join().unwrap();
        }
        let mut progress = progress.lock().unwrap();
        progress.finish(start + Duration::from_secs(6)).unwrap();
        assert_eq!(progress.batches, 800);
        assert_eq!(progress.sequences, 1600);
        assert_eq!(progress.bases, 800_000);
        let output = String::from_utf8(progress.writer.clone()).unwrap();
        assert_eq!(output.lines().count(), 2);
        assert!(output
            .lines()
            .all(|line| line.starts_with("Processed") || line.starts_with("Finished")));
    }
}
