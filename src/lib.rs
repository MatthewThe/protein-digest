use pyo3::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use pyo3::types::PyList;

fn parse_until_first_space(fasta_id: &str) -> &str {
    fasta_id.split(' ').next().unwrap_or(fasta_id)
}

fn swap_positions(seq: &mut Vec<char>, pos1: usize, pos2: usize) {
    seq.swap(pos1, pos2);
}

fn swap_special_aas(seq: &str, special_aas: &[char]) -> String {
    let mut chars: Vec<char> = seq.chars().collect();
    for i in 1..chars.len() {
        if special_aas.contains(&chars[i]) {
            swap_positions(&mut chars, i, i - 1);
        }
    }
    chars.into_iter().collect()
}

fn read_fasta_maxquant(
    file_path: &str,
    db: &str,
    special_aas: &[char],
    decoy_prefix: &str,
) -> Vec<(String, String)> {
    let file = File::open(file_path).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    let lines = reader.lines().map(|l| l.unwrap()).chain(std::iter::once(">".to_string()));

    for line in lines {
        let line = line.trim();
        if line.starts_with('>') {
            if let Some(ref name_val) = name {
                let sequence: String = seq.concat();

                if db == "target" || db == "concat" {
                    records.push((name_val.clone(), sequence.clone()));
                }

                if db == "decoy" || db == "concat" {
                    let mut rev_seq: String = sequence.chars().rev().collect();
                    if !special_aas.is_empty() {
                        rev_seq = swap_special_aas(&rev_seq, special_aas);
                    }
                    records.push((format!("{}{}", decoy_prefix, name_val), rev_seq));
                }
            }
            if line.len() > 1 {
                name = Some(parse_until_first_space(&line[1..]).to_string());
                seq = Vec::new();
            }
        } else {
            seq.push(line.to_string());
        }
    }
    records
}

#[pyfunction]
fn get_peptide_to_protein_map(
    fasta_file: &str,
    db: &str,
    min_len: usize,
    max_len: usize,
    pre: &Bound<'_, PyList>,
    not_post: &Bound<'_, PyList>,
    post: &Bound<'_, PyList>,
    digestion: &str,
    miscleavages: usize,
    methionine_cleavage: bool,
    special_aas: &Bound<'_, PyList>,
) -> PyResult<HashMap<String, Vec<String>>> {
    let mut peptide_to_protein_map: HashMap<String, Vec<String>> = HashMap::new();

    // Convert strings to chars for internal use
    let special_aas: Vec<char> = special_aas
        .extract::<Vec<String>>()?
        .into_iter()
        .flat_map(|s| s.chars().collect::<Vec<char>>())
        .collect();

    let records = read_fasta_maxquant(fasta_file, db, &special_aas, "REV__");

    for (idx, (protein, seq)) in records.iter().enumerate() {
        if idx % 10000 == 0 {
            println!("Digesting protein {}", idx);
        }

        let mut seen_peptides = HashSet::new();

        let peptides = get_digested_peptides(
            seq,
            min_len,
            max_len,
            &pre,
            &not_post,
            &post,
            digestion,
            miscleavages,
            methionine_cleavage,
        )?;

        for pep in peptides {
            if seen_peptides.insert(pep.clone()) {
                peptide_to_protein_map.entry(pep).or_default().push(protein.clone());
            }
        }
    }

    Ok(peptide_to_protein_map)
}


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
    m.add_function(wrap_pyfunction!(get_peptide_to_protein_map, m)?)?;
    Ok(())
}
