//! Amino-acid alphabets used to encode packed k-mers.

use clap::ValueEnum;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
/// Optional six-state amino-acid recoding schemes exposed by the CLI.
pub enum RecodeScheme {
    /// Classic Dayhoff six-group recoding.
    Dayhoff6,
    /// Susko-Roger six-group recoding.
    SR6,
    /// Kosiol-Goldman-Buttimore six-group recoding.
    KGB6,
}

impl std::fmt::Display for RecodeScheme {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Dayhoff6 => "Dayhoff6",
            Self::SR6 => "SR6",
            Self::KGB6 => "KGB6",
        })
    }
}

/// A compact, copyable description of the alphabet used throughout an analysis.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) struct Alphabet(Option<RecodeScheme>);

impl Alphabet {
    pub(crate) const fn new(recode: Option<RecodeScheme>) -> Self {
        Self(recode)
    }

    #[inline]
    pub(crate) fn encode(self, b: u8) -> u8 {
        match self.0 {
            None => AA20_CODES[b as usize],
            Some(RecodeScheme::Dayhoff6) => dayhoff6_code(b),
            Some(RecodeScheme::SR6) => sr6_code(b),
            Some(RecodeScheme::KGB6) => kgb6_code(b),
        }
    }

    /// Encode a byte already validated as an uppercase canonical amino acid.
    #[inline]
    pub(crate) fn encode_valid_uppercase(self, b: u8) -> u8 {
        match self.0 {
            None => b - b'A',
            Some(_) => self.encode(b),
        }
    }

    pub(crate) const fn n_states(self) -> usize {
        if self.0.is_some() {
            6
        } else {
            20
        }
    }

    pub(crate) const fn bits_per_state(self) -> u32 {
        if self.0.is_some() {
            3
        } else {
            5
        }
    }

    pub(crate) const fn max_packed_k(self) -> usize {
        64 / self.bits_per_state() as usize
    }

    pub(crate) const fn packed_mask(self, k: usize) -> u64 {
        let bits = k * self.bits_per_state() as usize;
        if bits >= 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        }
    }

    pub(crate) const fn symbol_mask(self) -> u64 {
        (1u64 << self.bits_per_state()) - 1
    }
}

const INVALID_STATE: u8 = 255;

const fn build_aa20_codes() -> [u8; 256] {
    let mut codes = [INVALID_STATE; 256];
    let canonical = *b"ACDEFGHIKLMNPQRSTVWY";
    let mut i = 0;
    while i < canonical.len() {
        let uppercase = canonical[i];
        let code = uppercase - b'A';
        codes[uppercase as usize] = code;
        codes[(uppercase + (b'a' - b'A')) as usize] = code;
        i += 1;
    }
    codes
}

static AA20_CODES: [u8; 256] = build_aa20_codes();

#[inline]
fn dayhoff6_code(b: u8) -> u8 {
    match b {
        b'C' | b'c' => 0,
        b'A' | b'a' | b'G' | b'g' | b'P' | b'p' | b'S' | b's' | b'T' | b't' => 1,
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'Q' | b'q' => 2,
        b'H' | b'h' | b'K' | b'k' | b'R' | b'r' => 3,
        b'I' | b'i' | b'L' | b'l' | b'M' | b'm' | b'V' | b'v' => 4,
        b'F' | b'f' | b'W' | b'w' | b'Y' | b'y' => 5,
        _ => 255,
    }
}

#[inline]
fn sr6_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' | b'P' | b'p' | b'S' | b's' | b'T' | b't' => 0,
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'G' | b'g' => 1,
        b'Q' | b'q' | b'K' | b'k' | b'R' | b'r' => 2,
        b'M' | b'm' | b'I' | b'i' | b'V' | b'v' | b'L' | b'l' => 3,
        b'W' | b'w' | b'C' | b'c' => 4,
        b'F' | b'f' | b'Y' | b'y' | b'H' | b'h' => 5,
        _ => 255,
    }
}

#[inline]
fn kgb6_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' | b'G' | b'g' | b'P' | b'p' | b'S' | b's' => 0,
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'Q' | b'q' | b'H' | b'h' | b'K' | b'k'
        | b'R' | b'r' | b'T' | b't' => 1,
        b'M' | b'm' | b'I' | b'i' | b'L' | b'l' => 2,
        b'W' | b'w' => 3,
        b'F' | b'f' | b'Y' | b'y' => 4,
        b'C' | b'c' | b'V' | b'v' => 5,
        _ => 255,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn aa20_has_twenty_distinct_case_insensitive_states() {
        let alphabet = Alphabet::new(None);
        let amino_acids = b"ACDEFGHIKLMNPQRSTVWY";
        assert_eq!(
            amino_acids
                .iter()
                .map(|&b| alphabet.encode(b))
                .collect::<Vec<_>>(),
            vec![0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 24]
        );
        assert_eq!(
            amino_acids
                .iter()
                .map(|&b| alphabet.encode(b))
                .collect::<HashSet<_>>()
                .len(),
            20
        );
        for &b in amino_acids {
            assert_eq!(alphabet.encode(b), alphabet.encode(b.to_ascii_lowercase()));
        }
    }

    #[test]
    fn aa20_packed_kmer_value_is_stable() {
        let alphabet = Alphabet::new(None);
        let packed = b"ACDE".iter().fold(0u64, |value, &b| {
            (value << alphabet.bits_per_state()) | alphabet.encode(b) as u64
        });
        assert_eq!(packed, 0x864);
    }

    #[test]
    fn alphabet_dimensions_and_invalid_residues() {
        let aa20 = Alphabet::new(None);
        assert_eq!(
            (aa20.n_states(), aa20.bits_per_state(), aa20.max_packed_k()),
            (20, 5, 12)
        );
        for scheme in [
            RecodeScheme::SR6,
            RecodeScheme::Dayhoff6,
            RecodeScheme::KGB6,
        ] {
            let alphabet = Alphabet::new(Some(scheme));
            assert_eq!(
                (
                    alphabet.n_states(),
                    alphabet.bits_per_state(),
                    alphabet.max_packed_k()
                ),
                (6, 3, 21)
            );
        }
        for b in b"BXZJUO-*?" {
            assert_eq!(aa20.encode(*b), 255);
        }
    }

    #[test]
    fn six_state_mappings_are_preserved() {
        let sr6 = Alphabet::new(Some(RecodeScheme::SR6));
        assert_eq!(
            b"ADEQMWFY"
                .iter()
                .map(|&b| sr6.encode(b))
                .collect::<Vec<_>>(),
            vec![0, 1, 1, 2, 3, 4, 5, 5]
        );
        let dayhoff = Alphabet::new(Some(RecodeScheme::Dayhoff6));
        assert_eq!(
            b"CADERIF"
                .iter()
                .map(|&b| dayhoff.encode(b))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 2, 3, 4, 5]
        );
        let kgb = Alphabet::new(Some(RecodeScheme::KGB6));
        assert_eq!(
            b"ADEIMWFVC"
                .iter()
                .map(|&b| kgb.encode(b))
                .collect::<Vec<_>>(),
            vec![0, 1, 1, 2, 2, 3, 4, 5, 5]
        );
    }
}
