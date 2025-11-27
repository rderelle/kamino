use smallvec::SmallVec;

/// Mixed-radix base: packed_prov = pid * PROV_M + seq
pub const PROV_M: u32 = 130_000;

pub const ENSURE_SUFFIX_NODE: bool = true;

pub type ProvVec = SmallVec<[u32; 4]>;
pub type SpeciesId = u16;
pub type ColorId   = u32;
pub type SpeciesSmall = SmallVec<[SpeciesId; 4]>;

/// Pack a (k−1)-mer left-shifted by `sym_bits` and add `sym`; then mask width for (k−1).
#[inline]
pub fn roll_pack_add_k1(prefix_k1: u64, sym: u32, mask_k1: u64, sym_bits: u32) -> u64 {
    ((prefix_k1 << sym_bits) & mask_k1) | (sym as u64)
}

/// Unpack provenance u32 → (pid, seq).
#[inline]
pub fn unpack_prov(packed: u32) -> (u32, u32) {
    (packed / PROV_M, packed % PROV_M)
}

