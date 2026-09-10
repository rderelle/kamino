//! Validation and compact packing for canonical amino acids.

pub(crate) const BITS_PER_SYMBOL: usize = 5;
pub(crate) const SYMBOL_MASK: u64 = 0x1f;
pub(crate) const MAX_PACKED_K: usize = 12;
pub(crate) const INVALID_STATE: u8 = 255;

const CANONICAL: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

const fn build_codes() -> [u8; 256] {
    let mut codes = [INVALID_STATE; 256];
    let mut i = 0;
    while i < CANONICAL.len() {
        let upper = CANONICAL[i];
        let code = upper - b'A';
        codes[upper as usize] = code;
        codes[(upper + (b'a' - b'A')) as usize] = code;
        i += 1;
    }
    codes
}

static CODES: [u8; 256] = build_codes();

#[inline]
pub(crate) fn encode(byte: u8) -> u8 {
    CODES[byte as usize]
}

#[inline]
pub(crate) fn encode_valid_uppercase(byte: u8) -> u8 {
    debug_assert_ne!(encode(byte), INVALID_STATE);
    byte - b'A'
}

#[inline]
pub(crate) fn kmer_mask(k: usize) -> u64 {
    debug_assert!((1..=MAX_PACKED_K).contains(&k));
    (1u64 << (BITS_PER_SYMBOL * k)) - 1
}

#[inline]
pub(crate) fn roll(value: u64, code: u8, mask: u64) -> u64 {
    ((value << BITS_PER_SYMBOL) | code as u64) & mask
}

#[inline]
pub(crate) fn overlap_k(anchor_k: usize) -> usize {
    MAX_PACKED_K.min(2 * anchor_k + 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn canonical_amino_acids_have_unique_case_insensitive_five_bit_codes() {
        let mut states = HashSet::new();
        for &aa in CANONICAL {
            let code = encode(aa);
            assert!(code <= SYMBOL_MASK as u8);
            assert_eq!(code, encode(aa.to_ascii_lowercase()));
            assert!(states.insert(code), "duplicate code for {}", aa as char);
        }
        assert_eq!(states.len(), 20);
    }

    #[test]
    fn ambiguous_and_non_amino_acid_bytes_are_invalid() {
        for &byte in b"BJOUXZ*-?.\0" {
            assert_eq!(encode(byte), INVALID_STATE);
        }
    }

    #[test]
    fn overlap_length_is_limited_by_path_and_packing_capacity() {
        for (anchor, expected) in [
            (1, 3),
            (2, 5),
            (3, 7),
            (4, 9),
            (5, 11),
            (6, 12),
            (8, 12),
            (12, 12),
        ] {
            assert_eq!(overlap_k(anchor), expected);
        }
    }
}
