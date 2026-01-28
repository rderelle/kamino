// Node-colored, node-based de Bruijn graph with variable alphabet.
//
// Nodes are (k-1)-mers packed with 3 bits per symbol (6-letter alphabet).
// Adjacency: node -> bitmask of outgoing symbols; mask is u32 (supports up to 32 symbols).
// Colors are stored as compact sorted species lists per node.

use bitvec::prelude::*;
use smallvec::SmallVec;

use crate::recode::{RecodeScheme, RECODE_BITS_PER_SYMBOL};

pub type AdjMask = u32;
pub type SpeciesSet = SmallVec<[u16; 4]>;

pub struct Graph {
    pub k: usize,
    pub recode_scheme: RecodeScheme,
    pub k_mask: u64,   // ((1 << (sym_bits*k)) - 1) or u64::MAX if overflow
    pub k1_mask: u64,  // ((1 << (sym_bits*(k-1))) - 1)
    pub adj: AdjTable, // node -> bitmask (compact storage)
    pub samples: NodeColorTable,
    pub n_species: usize, // number of species registered
    pub species_names: Vec<String>, // species names in exact sid order
                          //pub chain: HashMap<u64, Vec<u64>>,     // rep -> full chain of original nodes (including rep)
                          //pub end_kmer: HashMap<u64, u64>,       // rep -> last node in the chain
}

impl Graph {
    #[inline]
    fn compute_masks(k: usize) -> (u64, u64) {
        let kb = RECODE_BITS_PER_SYMBOL.saturating_mul(k as u32);
        let k1b = RECODE_BITS_PER_SYMBOL.saturating_mul((k - 1) as u32);
        let k_mask = if kb >= 64 { u64::MAX } else { (1u64 << kb) - 1 };
        let k1_mask = if k1b >= 64 {
            u64::MAX
        } else {
            (1u64 << k1b) - 1
        };
        (k_mask, k1_mask)
    }

    pub fn new(k: usize, recode_scheme: RecodeScheme) -> Self {
        let (k_mask, k1_mask) = Self::compute_masks(k);
        Self {
            k,
            recode_scheme,
            k_mask,
            k1_mask,
            adj: AdjTable::new(),
            samples: NodeColorTable::new(),
            n_species: 0,
            species_names: Vec::new(),
            //chain: HashMap::new(),
            //end_kmer: HashMap::new(),
        }
    }

    /// Optional: preset the expected number of species (kept for API compatibility).
    /// Names still must be registered via register_species to keep sid/name order aligned.
    pub fn init_species_len(&mut self, n: usize) {
        self.n_species = n;
    }

    /// Register a species name and return its SID (index in species_names).
    /// Call exactly once per species, in the same order files are processed.
    pub fn register_species(&mut self, name: String) -> usize {
        let sid = self.species_names.len();
        self.species_names.push(name);
        if sid + 1 > self.n_species {
            self.n_species = sid + 1;
        }
        sid
    }
}

#[derive(Clone)]
pub struct AdjTable {
    nodes: Vec<u64>,
    masks: MaskStorage,
}

#[derive(Clone)]
enum MaskStorage {
    U8(Vec<u8>),
}

impl AdjTable {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            masks: MaskStorage::U8(Vec::new()),
        }
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn nodes(&self) -> &[u64] {
        &self.nodes
    }

    pub fn get(&self, node: u64) -> Option<AdjMask> {
        self.nodes
            .binary_search(&node)
            .ok()
            .map(|idx| self.mask_at(idx))
    }

    pub fn iter(&self) -> AdjIter<'_> {
        AdjIter {
            table: self,
            index: 0,
        }
    }

    pub fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = (u64, u32)>,
    {
        let mut entries: Vec<(u64, u32)> = iter.into_iter().collect();
        entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        let mut nodes: Vec<u64> = Vec::with_capacity(entries.len());
        let mut masks_acc: Vec<u32> = Vec::with_capacity(entries.len());

        let mut iter = entries.into_iter();
        if let Some((mut cur_node, mut cur_mask)) = iter.next() {
            for (node, mask) in iter {
                if node == cur_node {
                    cur_mask |= mask;
                } else {
                    nodes.push(cur_node);
                    masks_acc.push(cur_mask);
                    cur_node = node;
                    cur_mask = mask;
                }
            }
            nodes.push(cur_node);
            masks_acc.push(cur_mask);
        }

        let masks = MaskStorage::from_masks(masks_acc);
        Self { nodes, masks }
    }

    fn mask_at(&self, idx: usize) -> AdjMask {
        let MaskStorage::U8(v) = &self.masks;
        v[idx] as AdjMask
    }
}

impl Default for AdjTable {
    fn default() -> Self {
        Self::new()
    }
}

impl MaskStorage {
    fn from_masks(masks: Vec<u32>) -> Self {
        let compact: Vec<u8> = masks.iter().map(|m| *m as u8).collect();
        MaskStorage::U8(compact)
    }
}

pub struct AdjIter<'a> {
    table: &'a AdjTable,
    index: usize,
}

impl<'a> Iterator for AdjIter<'a> {
    type Item = (u64, AdjMask);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.table.nodes.len() {
            return None;
        }
        let node = self.table.nodes[self.index];
        let mask = self.table.mask_at(self.index);
        self.index += 1;
        Some((node, mask))
    }
}

#[derive(Clone)]
pub struct NodeColorTable {
    nodes: Vec<u64>,
    offsets: Vec<u32>,
    lens: Vec<u32>,
    species: Vec<u32>,
    n_species: usize,
}

#[derive(Clone, Copy)]
pub struct NodeColorSlice<'a> {
    species: &'a [u32],
    n_species: usize,
}

pub struct NodeColorOnes<'a> {
    species: &'a [u32],
    index: usize,
}

impl NodeColorTable {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            offsets: Vec::new(),
            lens: Vec::new(),
            species: Vec::new(),
            n_species: 0,
        }
    }

    pub fn nodes(&self) -> &[u64] {
        &self.nodes
    }

    pub fn get(&self, node: u64) -> Option<NodeColorSlice<'_>> {
        self.nodes
            .binary_search(&node)
            .ok()
            .map(|idx| self.slice_at(idx))
    }

    pub fn from_iter<I>(iter: I, n_species: usize) -> Self
    where
        I: IntoIterator<Item = (u64, SpeciesSet)>,
    {
        let mut entries: Vec<(u64, SpeciesSet)> = iter.into_iter().collect();
        entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        let mut nodes: Vec<u64> = Vec::with_capacity(entries.len());
        let mut offsets: Vec<u32> = Vec::with_capacity(entries.len());
        let mut lens: Vec<u32> = Vec::with_capacity(entries.len());
        let mut species: Vec<u32> = Vec::new();

        let mut cur_node: Option<u64> = None;
        let mut cur_species: Vec<u32> = Vec::new();

        for (node, bits) in entries.into_iter() {
            if cur_node == Some(node) {
                Self::accumulate_species(&mut cur_species, &bits, n_species);
                continue;
            }

            if let Some(prev) = cur_node {
                Self::flush_node(
                    prev,
                    &mut cur_species,
                    &mut nodes,
                    &mut offsets,
                    &mut lens,
                    &mut species,
                );
            }

            cur_node = Some(node);
            Self::accumulate_species(&mut cur_species, &bits, n_species);
        }

        if let Some(last) = cur_node {
            Self::flush_node(
                last,
                &mut cur_species,
                &mut nodes,
                &mut offsets,
                &mut lens,
                &mut species,
            );
        }

        Self {
            nodes,
            offsets,
            lens,
            species,
            n_species,
        }
    }

    fn slice_at(&self, idx: usize) -> NodeColorSlice<'_> {
        let len = *self.lens.get(idx).unwrap_or(&0) as usize;
        if len == 0 {
            return NodeColorSlice {
                species: &[],
                n_species: self.n_species,
            };
        }
        let start = *self.offsets.get(idx).unwrap_or(&0) as usize;
        let end = start + len;
        NodeColorSlice {
            species: &self.species[start..end],
            n_species: self.n_species,
        }
    }

    fn accumulate_species(species: &mut Vec<u32>, bits: &SpeciesSet, n_species: usize) {
        if n_species == 0 {
            return;
        }
        for &sid in bits.iter() {
            let idx = sid as usize;
            if idx >= n_species {
                continue;
            }
            species.push(idx as u32);
        }
    }

    fn flush_node(
        node: u64,
        buf: &mut Vec<u32>,
        nodes: &mut Vec<u64>,
        offsets: &mut Vec<u32>,
        lens: &mut Vec<u32>,
        species: &mut Vec<u32>,
    ) {
        if !buf.is_empty() {
            buf.sort_unstable();
            buf.dedup();
        }
        nodes.push(node);
        offsets.push(species.len() as u32);
        lens.push(buf.len() as u32);
        species.extend_from_slice(buf);
        buf.clear();
    }
}

impl<'a> NodeColorSlice<'a> {
    pub fn iter_ones(self) -> NodeColorOnes<'a> {
        NodeColorOnes {
            species: self.species,
            index: 0,
        }
    }

    pub fn to_bitvec(self) -> BitVec<u32, Lsb0> {
        let mut bv = BitVec::<u32, Lsb0>::new();
        bv.resize(self.n_species, false);
        for &sid in self.species {
            if let Some(mut bit) = bv.get_mut(sid as usize) {
                *bit = true;
            }
        }
        bv
    }
}

impl<'a> Iterator for NodeColorOnes<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.species.len() {
            return None;
        }
        let sid = self.species[self.index] as usize;
        self.index += 1;
        Some(sid)
    }
}
