// index.rs
// Contains the "Index" and "Entry" structs, which describe how reference k-min-mers are stored. 

use crate::Kminmer;
use dashmap::DashMap;
use std::sync::Arc;
use std::hash::BuildHasherDefault;
use fxhash::FxHasher64;


// An Entry object holds information for a reference k-min-mer without storing the minimizer hashes themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    //pub id: String, // Reference ID
    pub counter: u64,
}
impl Entry {

    // Create a new Entry.
    pub fn new(counter: u64) -> Self {
        Entry {counter: counter}
    }

    // An empty Entry.
    pub fn empty() -> Self {
        Entry {counter: 0}
    }

    // Check if this Entry is Empty.
    pub fn is_empty(&self) -> bool {
        self.counter == 0
    }
}

// An Index object is a mapping of k-min-mer hashes (see kminmer.rs) to a single Entry (multiple Entries are not allowed).
pub struct Index {
    pub index: Arc<DashMap<u64, Entry, BuildHasherDefault<FxHasher64>>>
}
impl Index {

    // Create a new Index.
    pub fn new() -> Self {
        let hasher = BuildHasherDefault::<FxHasher64>::default();
        Index {index: Arc::new(DashMap::with_hasher(hasher))}
    }


    // Return the Entry associated with the k-min-mer hash h, or None if none.
    pub fn get(&self, h: &u64) -> Option<Entry> {
        let e = self.index.get(h);
        if let Some(r) = e {
            if !r.is_empty() {
                return Some(r.clone());
            }
        }
        None
    }

    // Add an Entry to the Index. If an Entry for the hash h already exists, insert None to prevent duplicates.
    pub fn add(&self, h: u64, counter: u64) {
        let e = self.index.insert(h, Entry::new(counter));
        if e.is_some() {self.index.insert(h, Entry::empty());}
    }

    pub fn increment(&self, h: u64) {
        let e_mut = self.index.get_mut(&h);
        if let Some(mut r) = e_mut
        {
            r.counter += 1;
        }
        else
        {
            self.index.insert(h, Entry::new(1));
        }
    }

}
