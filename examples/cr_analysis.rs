//! Run the CellRanger-style secondary analysis on a MatrixMarket directory.
//!
//! `cargo run --release --example cr_analysis -- <matrix_dir> <analysis_out_dir>`
use rustar_aligner::solo::cr_analysis;
use std::path::Path;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 3 {
        eprintln!("usage: cr_analysis <filtered_feature_bc_matrix dir> <analysis out dir>");
        std::process::exit(2);
    }
    let m = cr_analysis::read_mtx_dir(
        Path::new(&a[1]),
        "features.tsv",
        "barcodes.tsv",
        "matrix.mtx",
    )
    .expect("read matrix");
    let t = std::time::Instant::now();
    cr_analysis::run(&m, Path::new(&a[2])).expect("analysis");
    eprintln!("done in {:.1}s", t.elapsed().as_secs_f64());
}
