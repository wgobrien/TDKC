//! TDKC (Target Distilled K-mer Classifier) as a library.
//!
//! The crate was originally binary-only. Promoting the modules here lets the resident
//! [`engine::Engine`] be reused by the `tdkc serve` subcommand and by out-of-tree consumers
//! (e.g. PyO3 bindings) without duplicating the database load path.

pub mod build;
pub mod build_domain;
pub mod compression;
pub mod custom;
pub mod database;
pub mod engine;
pub mod hash;
pub mod inspect;
pub mod minimizer;
pub mod prep;
pub mod query;
pub mod serve;
pub mod taxonomy;
pub mod utils;

pub mod mphf;
pub mod querydb;
pub mod pgz;
