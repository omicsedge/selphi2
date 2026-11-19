//! Sample subsetting of a single-chromosome `.srp2` panel, SRP to SRP with no VCF/BCF in between.
//!
//! `--include-samples` keeps the listed samples, `--exclude-samples` drops them. The output keeps
//! the panel's own sample order, every site (also sites no kept haplotype carries), the original
//! variant IDs and the `##contig` line. Unknown IDs and an empty result are errors: a silently
//! ignored ID would leave a held-out sample in the panel. The file is written in a fresh hidden
//! directory next to `--out` and renamed into place when complete; the directory is removed on
//! any error or panic, so neither a truncated panel nor a stray `.srp2` is left behind.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use super::writer::{PanelVariant, SrpPanelWriter};
use super::SrpReader;

/// Sites read per batch (one tile stripe): bounds the dense rows to TILE_ROWS × kept haplotypes.
const BATCH: usize = super::TILE_ROWS;

/// Which samples a list selects.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SampleSelection {
    Include,
    Exclude,
}

/// Reads a sample list: one ID per line, surrounding whitespace and empty lines ignored.
pub fn read_sample_list(path: &Path) -> Result<Vec<String>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {}: {}", path.display(), e))?;
    Ok(text.lines().map(str::trim).filter(|l| !l.is_empty()).map(str::to_string).collect())
}

/// Indices (into `panel`) of the samples kept, in panel order.
pub fn kept_sample_indices(panel: &[String], list: &[String], mode: SampleSelection) -> Result<Vec<usize>, String> {
    let known: HashSet<&str> = panel.iter().map(String::as_str).collect();
    let unknown: Vec<&str> = list.iter().map(String::as_str).filter(|s| !known.contains(s)).collect();
    if !unknown.is_empty() {
        let shown: Vec<&str> = unknown.iter().take(5).copied().collect();
        return Err(format!("{} sample ID(s) in the list are not in the panel (first: {})", unknown.len(), shown.join(", ")));
    }
    let listed: HashSet<&str> = list.iter().map(String::as_str).collect();
    let kept: Vec<usize> = panel.iter().enumerate()
        .filter(|(_, s)| listed.contains(s.as_str()) == (mode == SampleSelection::Include))
        .map(|(i, _)| i)
        .collect();
    if kept.is_empty() {
        return Err("no sample would be left in the panel".into());
    }
    Ok(kept)
}

/// The `.srp2` path a subset is written to (the extension is added when missing).
pub fn subset_output_path(output: &str) -> PathBuf {
    super::helpers::with_panel_suffix(Path::new(output), super::PANEL_EXT)
}

/// True when `a` and `b` both exist and are the same file (same device and inode: hard and soft
/// links included), which a path comparison misses.
pub fn same_file(a: &Path, b: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    match (std::fs::metadata(a), std::fs::metadata(b)) {
        (Ok(x), Ok(y)) => x.dev() == y.dev() && x.ino() == y.ino(),
        _ => false,
    }
}

/// Writes `output` (.srp2) holding only the selected samples of the single-chromosome panel `source`.
/// Returns (samples kept, samples in the source).
pub fn subset_panel_samples(source: &Path, output: &Path, list: &[String], mode: SampleSelection)
    -> Result<(usize, usize), String> {
    let output = super::helpers::with_panel_suffix(output, super::PANEL_EXT);
    let name = output.file_name().ok_or("--out has no file name")?;
    if same_file(source, &output) {
        return Err(format!("{} is the source panel itself", output.display()));
    }
    let parent = output.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
    // Removed when dropped, also on a panic; panel directory scans do not enter it.
    let tmp_dir = tempfile::Builder::new().prefix(".selphi-subset-").tempdir_in(parent)
        .map_err(|e| format!("cannot create a temporary folder in {}: {}", parent.display(), e))?;
    let tmp = tmp_dir.path().join(name);
    let counts = write_subset(source, &tmp, list, mode)?;
    std::fs::rename(&tmp, &output).map_err(|e| format!("cannot move the new panel to {}: {}", output.display(), e))?;
    Ok(counts)
}

fn write_subset(source: &Path, output: &Path, list: &[String], mode: SampleSelection) -> Result<(usize, usize), String> {
    let mut reader = SrpReader::open(source, 0).map_err(|e| e.to_string())?;
    if !reader.load_tiled() {
        return Err(format!("{} is not a tiled panel; rebuild it with --prepare-reference-from", source.display()));
    }
    let kept = kept_sample_indices(&reader.sample_ids, list, mode)?;
    let n_var = reader.n_variants();
    let chrom = reader.chromosome().to_string();
    if let Some(i) = (0..n_var).find(|&i| reader.variants.get(i).chr != chrom) {
        return Err(format!("{} holds more than one chromosome ({} and {} at site {}); subset one chromosome at a time",
            source.display(), chrom, reader.variants.get(i).chr, i + 1));
    }
    let kept_names: Vec<String> = kept.iter().map(|&i| reader.sample_ids[i].clone()).collect();
    let kept_haps: Vec<usize> = kept.iter().flat_map(|&s| [2 * s, 2 * s + 1]).collect();

    let mut writer = {
        let pvs: Vec<PanelVariant> = (0..n_var).map(|i| {
            let v = reader.variants.get(i);
            PanelVariant { chrom: v.chr, pos: v.pos, ref_allele: v.ref_allele, alt_allele: v.alt_allele,
                           id: reader.original_ids.get(i).unwrap_or("") }
        }).collect();
        let contig = Some(reader.metadata.contig_field.as_str()).filter(|c| !c.is_empty());
        SrpPanelWriter::new_with_contig_field(&pvs, &kept_names, kept_haps.len(), output, contig)
            .map_err(|e| e.to_string())?
    };
    for base in (0..n_var).step_by(BATCH) {
        let idx: Vec<usize> = (base..(base + BATCH).min(n_var)).collect();
        let bm = reader.extract_ref_alleles_bitmatrix(&idx);
        let rows: Vec<Vec<u8>> = (0..idx.len()).into_par_iter()
            .map(|local| kept_haps.iter().map(|&h| bm.get(local, h) as u8).collect())
            .collect();
        for row in &rows { writer.push_row(row).map_err(|e| e.to_string())?; }
    }
    writer.finish().map_err(|e| e.to_string())?;
    Ok((kept.len(), reader.sample_ids.len()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::srp::writer::build_srp_from_panel;

    /// A panel of `n_samples` on chr22 where hap h carries ALT at variant v when (31v + 7h) % 5 == 0;
    /// every third site has an rsID, one site has a long REF, and the contig line carries a length.
    fn panel(dir: &Path, n_var: usize, n_samples: usize) -> (PathBuf, Vec<u8>) {
        let n_haps = 2 * n_samples;
        let phased: Vec<u8> = (0..n_var * n_haps).map(|i| (i / n_haps * 31 + i % n_haps * 7).is_multiple_of(5) as u8).collect();
        let chroms = vec!["22".to_string(); n_var];
        let refs: Vec<String> = (0..n_var).map(|v| if v == 7 { "ACGTACGTAC".into() } else { "A".into() }).collect();
        let alts = vec!["C".to_string(); n_var];
        let ids: Vec<String> = (0..n_var).map(|v| if v % 3 == 0 { format!("rs{}", 100 + v) } else { ".".into() }).collect();
        let pvs: Vec<PanelVariant> = (0..n_var).map(|v| PanelVariant {
            chrom: &chroms[v], pos: 1000 + 10 * v as i64, ref_allele: &refs[v], alt_allele: &alts[v], id: &ids[v],
        }).collect();
        let samples: Vec<String> = (0..n_samples).map(|s| format!("S{s}")).collect();
        let path = dir.join("panel.srp2");
        let mut w = SrpPanelWriter::new_with_contig_field(&pvs, &samples, n_haps, &path,
            Some("##contig=<ID=22,length=50818468>")).unwrap();
        for v in 0..n_var { w.push_row(&phased[v * n_haps..(v + 1) * n_haps]).unwrap(); }
        w.finish().unwrap();
        (path, phased)
    }

    /// Samples, per-site (chrom, pos, REF, ALT, original ID), the contig line and the allele bits.
    fn read_back(path: &Path) -> (Vec<String>, Vec<(String, i64, String, String, String)>, String, Vec<Vec<bool>>) {
        let mut r = SrpReader::open(path, 0).unwrap();
        r.load_tiled();
        let n = r.n_variants();
        let sites = (0..n).map(|i| { let v = r.variants.get(i);
            (v.chr.to_string(), v.pos, v.ref_allele.to_string(), v.alt_allele.to_string(), r.original_ids.get(i).unwrap_or("").to_string()) }).collect();
        let all: Vec<usize> = (0..n).collect();
        let bm = r.extract_ref_alleles_bitmatrix(&all);
        let rows = (0..n).map(|v| (0..r.n_haps()).map(|h| bm.get(v, h)).collect()).collect();
        (r.sample_ids.clone(), sites, r.metadata.contig_field.clone(), rows)
    }

    /// Exclude keeps the other samples' haplotypes bit for bit, in panel order, with every site and its
    /// metadata unchanged; include of the complement writes the same file. Two shapes: many sites
    /// (more than the writer's 16-stripe flush) and many haplotypes (past 64-bit words and 4,096-hap tiles).
    #[test]
    fn exclude_and_include_of_the_complement_keep_the_source_as_it_was() {
        for (n_var, n_samples, drop) in [(20_000usize, 6usize, vec![4usize, 1]), (1_100, 2_100, vec![0, 31, 32, 2047, 2099])] {
            let dir = tempfile::tempdir().unwrap();
            let (src, phased) = panel(dir.path(), n_var, n_samples);
            let n_haps = 2 * n_samples;
            let drop_ids: Vec<String> = drop.iter().map(|s| format!("S{s}")).collect();
            let keep: Vec<usize> = (0..n_samples).filter(|s| !drop.contains(s)).collect();
            let keep_ids: Vec<String> = keep.iter().rev().map(|s| format!("S{s}")).collect();
            let (ex, inc) = (dir.path().join("ex.srp2"), dir.path().join("inc.srp2"));

            assert_eq!(subset_panel_samples(&src, &ex, &drop_ids, SampleSelection::Exclude).unwrap(), (keep.len(), n_samples));
            assert_eq!(subset_panel_samples(&src, &inc, &keep_ids, SampleSelection::Include).unwrap(), (keep.len(), n_samples));

            let (src_names, src_sites, src_contig, _) = read_back(&src);
            let (names, sites, contig, rows) = read_back(&ex);
            let want: Vec<String> = keep.iter().map(|&s| src_names[s].clone()).collect();
            assert_eq!(names, want);
            assert_eq!(sites, src_sites);
            assert_eq!(contig, src_contig);
            assert_eq!(contig, "##contig=<ID=22,length=50818468>");
            assert_eq!(sites[3].4, "rs103");
            for v in 0..n_var {
                for (j, &s) in keep.iter().enumerate() {
                    for k in 0..2 {
                        assert_eq!(rows[v][2 * j + k], phased[v * n_haps + 2 * s + k] != 0, "variant {v}, sample {s}, hap {k}");
                    }
                }
            }
            assert_eq!(std::fs::read(&ex).unwrap(), std::fs::read(&inc).unwrap());
        }
    }

    #[test]
    fn an_unknown_id_is_an_error_and_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let (src, _) = panel(dir.path(), 10, 3);
        let out = dir.path().join("o.srp2");
        let err = subset_panel_samples(&src, &out, &["S9".into()], SampleSelection::Exclude).unwrap_err();
        assert!(err.contains("not in the panel") && err.contains("S9"), "{err}");
        let left: Vec<_> = std::fs::read_dir(dir.path()).unwrap().map(|e| e.unwrap().file_name()).collect();
        assert_eq!(left, [std::ffi::OsString::from("panel.srp2")]);
    }

    /// A hard link to the source at --out is the source: refused, and the source keeps its bytes.
    #[test]
    fn a_hard_link_to_the_source_at_out_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let (src, _) = panel(dir.path(), 50, 3);
        let before = std::fs::read(&src).unwrap();
        let link = dir.path().join("link.srp2");
        std::fs::hard_link(&src, &link).unwrap();
        assert!(subset_panel_samples(&src, &link, &["S1".into()], SampleSelection::Exclude).is_err());
        assert_eq!(std::fs::read(&src).unwrap(), before);
    }

    /// A failure after the panel is written (here the rename, --out being a folder) is an error and
    /// leaves nothing next to --out.
    #[test]
    fn a_failed_rename_is_an_error_and_leaves_no_temporary_file() {
        let dir = tempfile::tempdir().unwrap();
        let (src, _) = panel(dir.path(), 50, 3);
        let out = dir.path().join("taken.srp2");
        std::fs::create_dir(&out).unwrap();
        std::fs::write(out.join("keep"), "x").unwrap();
        let err = subset_panel_samples(&src, &out, &["S1".into()], SampleSelection::Exclude).unwrap_err();
        assert!(err.contains("cannot move"), "{err}");
        let mut left: Vec<_> = std::fs::read_dir(dir.path()).unwrap().map(|e| e.unwrap().file_name()).collect();
        left.sort();
        assert_eq!(left, [std::ffi::OsString::from("panel.srp2"), "taken.srp2".into()]);
    }

    /// A container whose sites name two chromosomes is refused rather than relabelled.
    #[test]
    fn a_panel_with_two_chromosomes_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let chroms = ["21".to_string(), "21".to_string(), "22".to_string()];
        let (a, c) = ("A".to_string(), "C".to_string());
        let pvs: Vec<PanelVariant> = (0..3).map(|v| PanelVariant {
            chrom: &chroms[v], pos: 100 + v as i64, ref_allele: &a, alt_allele: &c, id: ".",
        }).collect();
        let src = dir.path().join("mixed.srp2");
        build_srp_from_panel(&[0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1], &pvs, &["S0".into(), "S1".into()], 4, &src).unwrap();
        let err = subset_panel_samples(&src, &dir.path().join("o.srp2"), &["S1".into()], SampleSelection::Exclude).unwrap_err();
        assert!(err.contains("more than one chromosome"), "{err}");
    }

    #[test]
    fn the_source_itself_at_out_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let (src, _) = panel(dir.path(), 10, 3);
        let before = std::fs::read(&src).unwrap();
        assert!(subset_panel_samples(&src, &src, &["S1".into()], SampleSelection::Exclude).is_err());
        assert_eq!(std::fs::read(&src).unwrap(), before);
    }

    #[test]
    fn an_empty_result_is_an_error() {
        let panel_ids: Vec<String> = vec!["A".into(), "B".into()];
        assert!(kept_sample_indices(&panel_ids, &panel_ids, SampleSelection::Exclude).is_err());
        assert!(kept_sample_indices(&panel_ids, &[], SampleSelection::Include).is_err());
    }

    #[test]
    fn the_list_ignores_blank_lines_and_spaces() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("l.txt");
        std::fs::write(&p, " S1 \n\nS2\r\n").unwrap();
        assert_eq!(read_sample_list(&p).unwrap(), ["S1", "S2"]);
    }
}
