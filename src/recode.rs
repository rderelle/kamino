//! Amino-acid recoding schemes (6-letter alphabet, 3 bits/symbol).

use clap::ValueEnum;

pub const RECODE_ALPHABET_SIZE: usize = 6;
pub const RECODE_BITS_PER_SYMBOL: u32 = 3;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum RecodeScheme {
    Dayhoff6,
    SR6,
    KGB6,
}

impl std::fmt::Display for RecodeScheme {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            RecodeScheme::Dayhoff6 => "Dayhoff6",
            RecodeScheme::SR6 => "SR6",
            RecodeScheme::KGB6 => "KGB6",
        };
        f.write_str(name)
    }
}

/// Return class code (0..=5) or 255 for unknown/ambiguous.
#[inline]
pub fn recode_byte(b: u8, scheme: RecodeScheme) -> u8 {
    match scheme {
        RecodeScheme::Dayhoff6 => dayhoff6_code(b),
        RecodeScheme::SR6 => sr6_code(b),
        RecodeScheme::KGB6 => kgb6_code(b),
    }
}

#[inline]
fn dayhoff6_code(b: u8) -> u8 {
    match b {
        b'C' | b'c' => 0,                                                         // C
        b'A' | b'a' | b'G' | b'g' | b'P' | b'p' | b'S' | b's' | b'T' | b't' => 1, // A/G/P/S/T
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'Q' | b'q' => 2,               // D/E/N/Q
        b'H' | b'h' | b'K' | b'k' | b'R' | b'r' => 3,                             // H/K/R
        b'I' | b'i' | b'L' | b'l' | b'M' | b'm' | b'V' | b'v' => 4,               // I/L/M/V
        b'F' | b'f' | b'W' | b'w' | b'Y' | b'y' => 5,                             // F/W/Y
        _ => 255,
    }
}

#[inline]
fn sr6_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' | b'P' | b'p' | b'S' | b's' | b'T' | b't' => 0, // A/P/S/T
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'G' | b'g' => 1, // D/E/N/G
        b'Q' | b'q' | b'K' | b'k' | b'R' | b'r' => 2,               // Q/K/R
        b'M' | b'm' | b'I' | b'i' | b'V' | b'v' | b'L' | b'l' => 3, // M/I/V/L
        b'W' | b'w' | b'C' | b'c' => 4,                             // W/C
        b'F' | b'f' | b'Y' | b'y' | b'H' | b'h' => 5,               // F/Y/H
        _ => 255,
    }
}

#[inline]
fn kgb6_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' | b'G' | b'g' | b'P' | b'p' | b'S' | b's' => 0, // A/G/P/S
        b'D' | b'd' | b'E' | b'e' | b'N' | b'n' | b'Q' | b'q' | b'H' | b'h' | b'K' | b'k'
        | b'R' | b'r' | b'T' | b't' => 1, // D/E/N/Q/H/K/R/T
        b'M' | b'm' | b'I' | b'i' | b'L' | b'l' => 2,               // M/I/L
        b'W' | b'w' => 3,                                           // W
        b'F' | b'f' | b'Y' | b'y' => 4,                             // F/Y
        b'C' | b'c' | b'V' | b'v' => 5,                             // C/V
        _ => 255,
    }
}
