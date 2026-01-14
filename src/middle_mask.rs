pub(crate) fn mask_middle_if_diff_run(block: &mut [Vec<u8>], k1: usize, m: usize) {
    if block.is_empty() {
        return;
    }
    let bl = block[0].len();
    let middle_start = k1;
    let middle_end = bl - k1;
    let middle_len = middle_end - middle_start;
    if middle_len < m {
        return;
    }

    let mut consensus: Vec<u8> = vec![b'X'; middle_len];
    for (offset, slot) in consensus.iter_mut().enumerate() {
        let col = middle_start + offset;
        let mut counts = [0usize; 256];
        let mut valid_votes = 0usize;
        for row in block.iter() {
            let b = row[col];
            if b == b'-' || b == b'X' {
                continue;
            }
            counts[b as usize] += 1;
            valid_votes += 1;
        }
        if valid_votes == 0 {
            *slot = b'X';
            continue;
        }
        let mut best = b'X';
        let mut best_count = 0usize;
        for (aa, &count) in counts.iter().enumerate() {
            if count > best_count {
                best_count = count;
                best = aa as u8;
            }
        }
        // Majority-rule consensus: require ≥50% of valid votes (not strict).
        if best_count * 2 >= valid_votes {
            *slot = best;
        } else {
            *slot = b'X';
        }
    }

    for row in block.iter_mut() {
        // '-' or 'X' in either the sequence or the consensus is ignored: no difference
        // and reset the consecutive-mismatch counter to 0.
        let mut run = 0usize;
        let mut should_mask = false;
        for col in middle_start..middle_end {
            let a = row[col];
            let c = consensus[col - middle_start];
            
            if a == b'-' || a == b'X' || c == b'-' || c == b'X' {
                run = 0;
                continue;
            }
            
            if a != c {
                run += 1;
            } else {
                run = 0;
            }
            if run >= m {
                should_mask = true;
                break;
            }
        }
        
        if should_mask {
            for cell in row
                .iter_mut()
                .take(middle_end)
                .skip(middle_start)
            {
                *cell = b'X';
            }
        }
    }
}
