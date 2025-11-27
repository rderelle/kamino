//! Dayhoff6 amino-acid recoding (6-letter alphabet).

pub const DAYHOFF6_ALPHABET_SIZE: usize = 6;
pub const DAYHOFF6_BITS_PER_SYMBOL: u32 = 3;

/// Return class code (0..DAYHOFF6_ALPHABET_SIZE-1) or 255 for unknown/ambiguous.
#[inline]
pub fn recode_byte(b: u8) -> u8 {
    dayhoff6_code(b)
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
