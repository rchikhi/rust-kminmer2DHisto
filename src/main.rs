// 2D histograms using kminmers v0.1.0
// Copyright 2020-2021 Baris Ekim, Rayan Chikhi.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.

#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
#![feature(iter_advance_by)]
use indicatif::ProgressBar;
use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use std::fs::{File};
use std::fs;
use structopt::StructOpt;
use std::sync::{Arc, Mutex};
use std::path::PathBuf;
use std::time::{Instant};
use std::mem::{MaybeUninit};
use seq_io::BaseRecord;
use lzzzz::lz4f::{WriteCompressor, BufReadDecompressor, Preferences};
use flate2::read::GzDecoder;
use dashmap::DashMap;
use std::cell::UnsafeCell;
use std::io::Result;
use core::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use crate::index::{Entry, Index};
use rust_seq2kminmers::Kminmer;
mod index;
mod closures;
mod mers;

type ThreadIdType = usize;
pub struct Params {
    k: usize,
    l: usize,
    density: f64,
}

/// Try to get memory usage (resident set size) in bytes using the `getrusage()` function from libc.
// from https://github.com/digama0/mm0/blob/bebd670c5a77a1400913ebddec2c6248e76f90fe/mm0-rs/src/util.rs
fn get_memory_rusage() -> usize {
  let usage = unsafe {
    let mut usage = MaybeUninit::uninit();
    assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
    usage.assume_init()
  };
  usage.ru_maxrss as usize * 1024
}

fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :Box<dyn BufRead + Send> = match filetype { 
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))), 
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ =>     Box::new(BufReader::new(file)), 
    }; 
    reader
}
#[derive(Debug, StructOpt)]
#[structopt(name = "kminmer2Dhisto")]
/// Original implementation of hifimap, a fast HiFi read mapper.
struct Opt {
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    ///
    /// Input file can be FASTA/FASTQ, as well as gzip-compressed (.gz) or
    /// lz4-compressed (.lz4). Lowercase bases are currently not supported;
    /// see documentation for formatting.
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    /// Output prefix 
    ///
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    ///
    /// The length of each node of the mdBG. If
    /// fewer l-mers than this value are obtained
    /// from a read, they will be ignored.
    #[structopt(short, long)]
    k: Option<usize>,
    /// l-mer (minimizer) length
    ///
    /// The length of each minimizer selected using
    /// the minimizer scheme from base-space sequences.
    #[structopt(short, long)]
    l: Option<usize>,
    /// Density threshold for density-based selection scheme
    /// 
    /// The density threshold is analogous to the
    /// fraction of l-mers that will be selected as
    /// minimizers from a read.
    #[structopt(short, long)]
    density: Option<f64>,
    /// Reference genome input
    ///
    /// Reference to be indexed and mapped to. 
    /// Allows multi-line FASTA and
    /// doesn't filter any kminmers.
    #[structopt(parse(from_os_str), long)]
    reference: Option<PathBuf>,
    /// Number of threads
    #[structopt(long)]
    threads: Option<usize>,
}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut filename = PathBuf::new();
    let mut ref_filename = PathBuf::new();
    let mut output_prefix;
    let mut k : usize = 5;
    let mut l : usize = 31;
    let mut density : f64 = 0.01;
    let mut threads : usize = 8;
    if opt.reads.is_some() {filename = opt.reads.unwrap();} 
    if opt.reference.is_some() {ref_filename = opt.reference.unwrap();} 
    if filename.as_os_str().is_empty() {panic!("Please specify an input file.");}
    if ref_filename.as_os_str().is_empty() {panic!("Please specify a reference file.");}
    let filename_str = filename.to_str().unwrap();
    let mut reads_are_fasta : bool = false;
    let mut ref_is_fasta    : bool = false;
    if filename_str.contains(".fasta.") || filename_str.contains(".fa.") || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") {
        reads_are_fasta = true;
        println!("Input file: {}", filename_str);
        println!("Format: FASTA");
    }
    let ref_filename_str = ref_filename.to_str().unwrap();
    if ref_filename_str.contains(".fasta.") || ref_filename_str.contains(".fa.") || ref_filename_str.ends_with(".fa") || ref_filename_str.ends_with(".fasta") {
        ref_is_fasta = true;
        println!("Reference file: {}", ref_filename_str);
        println!("Format: FASTA");
    }
    if opt.k.is_some() {k = opt.k.unwrap()} else {println!("Warning: Using default k value ({}).", k);} 
    if opt.l.is_some() {l = opt.l.unwrap()} else {println!("Warning: Using default l value ({}).", l);}
    if opt.density.is_some() {density = opt.density.unwrap()} else {println!("Warning: Using default density value ({}%).", density * 100.0);}
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    output_prefix = PathBuf::from(format!("2DHisto-k{}-d{}-l{}", k, density, l));
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
 
    let params = Params { 
        k,
        l,
        density,
    };

    let ref_threads = threads;
    let ref_queue_len = threads;
    let queue_len = 200; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing

    closures::run_mers(&filename, &ref_filename, &params, ref_threads, threads, ref_queue_len, queue_len, reads_are_fasta, ref_is_fasta, &output_prefix);
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}
