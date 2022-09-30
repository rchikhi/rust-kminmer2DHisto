// mers.rs

use crate::{File, Kminmer, Index, Params};
use std::borrow::Cow;
use std::cmp;
use std::collections::{hash_map::DefaultHasher, HashMap, HashSet, VecDeque};
use std::hash::{Hash, Hasher};
use std::io::Write;
use dashmap::{DashMap, DashSet};
use rust_seq2kminmers::KminmersIterator;

// Extract k-min-mers from reference. We don't store k-min-mer objects or hashes in a Vec, but rather immediately insert into the Index.
pub fn ref_extract(seq_id: &str, inp_seq_raw: &[u8], params: &Params, ref_mers_index: &Index) -> usize {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return 0;
    }
    let density = params.density;
    let iter = KminmersIterator::new(inp_seq_raw, l, k, density, false).unwrap();
    let mut count = 0;
    for kminmer in iter {
        // Add a reference k-min-mer to the Index.
        ref_mers_index.increment(kminmer.get_hash_u64());
        count += 1;
    }
    count
}

// Extract k-min-mers from the query. We need to store Kminmer objects for the query in order to compute Hits.
pub fn extract<'a>(seq_id: &str, inp_seq_raw: &'a [u8], params: &Params) -> Option<KminmersIterator<'a>> {
    let l = params.l;
    let k = params.k;
    if inp_seq_raw.len() < l+k-1 {
        return None;
    }
    let density = params.density;
    return Some(KminmersIterator::new(inp_seq_raw, l, k, density, false).unwrap());
}

// populate the hashtable with read kminmers
pub fn insert_kminmers(query_id: &str, query_it_raw: &mut Option<KminmersIterator>, index: &Index, params: &Params, q_len: usize)  {
    let l = params.l;
    let k = params.k;
    if query_it_raw.is_none() {return;}
    let mut query_it = query_it_raw.as_mut().unwrap();
    while let Some(q) = query_it.next() {
        index.increment(q.get_hash_u64());
    }
}


pub fn process_read(q_id: &str, q_len: usize, q_str: &[u8], ref_lens: &DashMap<String, usize>, read_mers_index: &Index, params: &Params) {
    let mut kminmers = extract(q_id, q_str, params);
    insert_kminmers(q_id, &mut kminmers, read_mers_index, params, q_len);
}
