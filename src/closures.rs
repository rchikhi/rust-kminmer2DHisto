// closures.rs
// Functions for FASTA parsing and invoking all main functions 

use std::io::{self};
use std::error::Error;
use std::io::{BufRead, BufReader};
use std::path::Path;
use crate::BufReadDecompressor;
use std::fs::{File};
use std::sync::{Arc};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::DashMap;
use super::mers;
use std::path::PathBuf;
use super::Params;
use crate::get_reader;
use indicatif::ProgressBar;
use std::time::Instant;
use dashmap::DashSet;
use crate::index::{Entry, Index};
use std::borrow::Cow;
use std::io::Write;


// Main function for all FASTA parsing + mapping / alignment functions.
pub fn run_mers(filename: &PathBuf, ref_filename: &PathBuf, params: &Params, ref_threads: usize, threads: usize, ref_queue_len: usize, queue_len: usize, reads_are_fasta: bool, ref_is_fasta: bool, output_prefix: &PathBuf) {

    let ref_mers_index = Index::new(); // Index of reference k-min-mer entries
    let read_mers_index = Index::new(); // Index of read k-min-mer entries
    let lens : DashMap<String, usize> = DashMap::new(); // Sequence lengths per reference

    // Closure for indexing reference k-min-mers
    let index_mers = |seq_id: &str, seq: &[u8], params: &Params| -> usize {
        let nb_mers = mers::ref_extract(seq_id, seq, params, &ref_mers_index);
        lens.insert(seq_id.to_string(), seq.len());
        nb_mers
    };

    // Closures for obtaining k-min-mers from references

    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        let nb_mers = index_mers(ref_id, ref_str, params);
        println!("Indexed reference {}: {} k-min-mers.", ref_id, nb_mers);
        return Some(1)
    };

    let ref_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq().to_vec(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);

    };
    let ref_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap().to_string();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);
    };
    let ref_main_thread_mer = |found: &mut Option<u64>| { // runs in main thread
        None::<()>
    };

    // Closures for mapping queries to references

    let query_process_read_aux_mer = |seq_str: &[u8], seq_id: &str| -> bool {
        mers::process_read(&seq_id, seq_str.len(), &seq_str, &lens, &read_mers_index, params);
        return true;
    };
    let query_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut bool| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    
    };
    let query_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut bool| {
        let seq_str = record.seq(); 
        let seq_id = record.id().unwrap().to_string();
        *found = query_process_read_aux_mer(&seq_str, &seq_id);
    };
    let mut main_thread_mer = |found: &mut bool| { // runs in main thread
        None::<()>
    };

    // Start processing references

    let start = Instant::now();
    let buf = get_reader(&ref_filename);
    if ref_is_fasta {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, ref_threads as u32, ref_queue_len, ref_process_read_fasta_mer, |record, found| {ref_main_thread_mer(found)});
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, ref_threads as u32, ref_queue_len, ref_process_read_fastq_mer, |record, found| {ref_main_thread_mer(found)});
    }
    let duration = start.elapsed();
    println!("Indexed references in {:?}.", duration);

    // Done, start processing reads

    let query_start = Instant::now();
    let buf = get_reader(&filename);
    if reads_are_fasta {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, query_process_read_fasta_mer, |record, found| {main_thread_mer(found)});
        let query_duration = query_start.elapsed();
        println!("Processed reads in {:?}.", query_duration);
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, query_process_read_fastq_mer, |record, found| {main_thread_mer(found)});
        let query_duration = query_start.elapsed();
        println!("Processed reads in {:?}.", query_duration);
    }


    // Now produce the 2D histogram by iterating read kmers
    let mut hist = vec![vec![0u64; 10]; 10000];

    let hist_path = format!("{}{}", output_prefix.to_str().unwrap(), ".hist2D");
    let mut hist_file = match File::create(&hist_path) {
        Err(why) => panic!("Couldn't create {}: {}", hist_path, why.description()),
        Ok(hist_file) => hist_file,
    };

    println!("nb read kminmers {}",read_mers_index.index.len());
    println!("nb ref kminmers {}",ref_mers_index.index.len());

    for item in read_mers_index.index.iter() {
        let (node, entry) = item.pair();
        let kminmer_abundance = entry.counter;
        let ref_e = ref_mers_index.get(node);
        let ref_abundance = if let Some(m) = ref_e {
           m.counter
        } else {0};
        let i = if kminmer_abundance > 9999 { 9999 } else { kminmer_abundance } as usize;
        let j = if ref_abundance > 9 { 9 } else { ref_abundance } as usize;
        hist[i][j] += 1;
    } 

    // now do the edge case where reference kminmers aren't found in the reads
    for item in ref_mers_index.index.iter() {
        let (node, entry) = item.pair();
        let ref_abundance = entry.counter;
        let read_e = read_mers_index.get(node);
        let read_abundance = if let Some(m) = read_e {
           m.counter
        } else {0};
        if read_abundance == 0
        {
            let i = 0;
            let j = if ref_abundance > 9 { 9 } else { ref_abundance } as usize;
            hist[i][j] += 1;
        }
    } 
 
    for i in 0..10000 {
        for j in 0..10 {
            write!(hist_file, "{}\t", hist[i][j]).expect("Error writing hist file.");
        }
        write!(hist_file, "\n").expect("Error writing hist file.");
    }

}
