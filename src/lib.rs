use pyo3::prelude::*;
use pyo3::types::PyList;

/// Check if a cleavage site is enzymatic
fn is_enzymatic(aa1: char, aa2: char, pre: &[char], not_post: &[char], post: &[char]) -> bool {
    (pre.contains(&aa1) && !not_post.contains(&aa2)) || post.contains(&aa2)
}

/// Non-specific digestion
fn non_specific_digest(seq: &str, min_len: usize, max_len: usize) -> Vec<String> {
    let mut peptides = Vec::new();
    let seq_len = seq.len();

    for i in 0..=seq_len {
        for j in (i + min_len)..=(usize::min(seq_len, i + max_len)) {
            if j <= seq_len {
                peptides.push(seq[i..j].to_string());
            }
        }
    }
    peptides
}

/// Semi-specific digestion
fn semi_specific_digest(
    seq: &str,
    min_len: usize,
    max_len: usize,
    pre: &[char],
    not_post: &[char],
    post: &[char],
    miscleavages: usize,
    methionine_cleavage: bool,
) -> Vec<String> {
    let seq_chars: Vec<char> = seq.chars().collect();
    let seq_len = seq_chars.len();
    let mut peptides = Vec::new();
    let mut starts = vec![0];
    let methionine_cleavage = methionine_cleavage && seq_chars[0] == 'M';

    for i in 0..=seq_len-1 {
        let is_cleavage_site = i == seq_len-1
            || is_enzymatic(
                seq_chars[usize::min(seq_len - 1, i)],
                seq_chars[usize::min(seq_len - 1, i + 1)],
                pre,
                not_post,
                post,
            )
            || (i == 0 && methionine_cleavage);

        if is_cleavage_site {
            let start = starts[0];
            for j in start..std::cmp::min(i, seq_len - 1) {
                let pep_len = std::cmp::min(i, seq_len - 1) - j + 1;
                if (min_len..=max_len).contains(&pep_len) {
                    // assuming you're collecting peptides in a Vec
                    peptides.push(seq[j..=i].to_string()); // inclusive range like Python's [j:i+1]
                }
            }

            starts.push(i + 1);

            let methionine_cleaved = if starts[0] == 0 && methionine_cleavage {
                1
            } else {
                0
            };

            if starts.len() > miscleavages + 1 + methionine_cleaved || i == seq_len {
                starts.drain(0..1 + methionine_cleaved); // removes from front
            }
        } else {
            for &start in &starts {
                let pep_len = i - start + 1;
                if (min_len..=max_len).contains(&pep_len) && !starts.contains(&(i + 1)) {
                    peptides.push(seq[start..=i].to_string());
                }
            }
        }
    }

    peptides
}

/// Full digestion
fn full_digest(
    seq: &str,
    min_len: usize,
    max_len: usize,
    pre: &[char],
    not_post: &[char],
    post: &[char],
    miscleavages: usize,
    methionine_cleavage: bool,
) -> Vec<String> {
    let seq_chars: Vec<char> = seq.chars().collect();
    let seq_len = seq_chars.len();
    let mut peptides = Vec::new();
    let mut starts = vec![0];
    let methionine_cleavage = methionine_cleavage && seq_chars[0] == 'M';

    let check_pre = !pre.is_empty();
    let check_post = !post.is_empty();

    let mut cleavage_sites: Vec<usize> = if methionine_cleavage { vec![0] } else { vec![] };

    for i in 0..=seq_len-1 {
        if (check_pre
            && pre.contains(&seq_chars[i])
            && !not_post.contains(&seq_chars[usize::min(seq_len - 1, i + 1)]))
            || (check_post && post.contains(&seq_chars[usize::min(seq_len - 1, i + 1)]))
        {
            cleavage_sites.push(i);
        }
    }
    cleavage_sites.push(seq_len - 1);

    for &i in &cleavage_sites {
        for &start in &starts {
            let pep_len = 1 + i - start;
            if (min_len..=max_len).contains(&pep_len) {
                peptides.push(seq[start..=i].to_string());
            }
        }
        starts.push(i + 1);
        let methionine_cleaved = usize::from(starts[0] == 0 && methionine_cleavage);
        if starts.len() > miscleavages + 1 + methionine_cleaved {
            starts = starts.split_off(1 + methionine_cleaved);
        }
    }

    peptides
}

/// Python-exposed function
#[pyfunction]
fn get_digested_peptides(
    seq: &str,
    min_len: usize,
    max_len: usize,
    pre: &Bound<'_, PyList>,
    not_post: &Bound<'_, PyList>,
    post: &Bound<'_, PyList>,
    digestion: &str,
    miscleavages: usize,
    methionine_cleavage: bool,
) -> PyResult<Vec<String>> {
    let pre: Vec<char> = pre
        .extract::<Vec<String>>()?
        .into_iter()
        .flat_map(|s| s.chars().collect::<Vec<char>>())
        .collect();
    let not_post: Vec<char> = not_post
        .extract::<Vec<String>>()?
        .into_iter()
        .flat_map(|s| s.chars().collect::<Vec<char>>())
        .collect();
    let post: Vec<char> = post
        .extract::<Vec<String>>()?
        .into_iter()
        .flat_map(|s| s.chars().collect::<Vec<char>>())
        .collect();

    let peptides = match digestion {
        "none" => non_specific_digest(seq, min_len, max_len),
        "semi" => semi_specific_digest(
            seq,
            min_len,
            max_len,
            &pre,
            &not_post,
            &post,
            miscleavages,
            methionine_cleavage,
        ),
        _ => full_digest(
            seq,
            min_len,
            max_len,
            &pre,
            &not_post,
            &post,
            miscleavages,
            methionine_cleavage,
        ),
    };

    Ok(peptides)
}

#[pymodule]
fn protein_digest(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_digested_peptides, m)?)?;
    Ok(())
}
