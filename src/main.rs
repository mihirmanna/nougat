use anyhow::{Context, Result};
use nougat::config::Config;
use nougat::run;
use std::env;

fn main() -> Result<()> {
    let path = env::args()
        .nth(1)
        .context("Usage: ./nougat <config.toml>")?;
    let mut config = Config::from_file(path)?;
    let results = run(&mut config);

    println!("{:<18} {}", "Photons:", config.photon.n_photons);
    println!(
        "{:<18} {} ({:.1}%)",
        "Escaped top:",
        results.n_escaped_top,
        100.0 * results.transmittance
    );
    println!(
        "{:<18} {} ({:.1}%)",
        "Escaped bottom:",
        results.n_escaped_bottom,
        100.0 * results.n_escaped_bottom as f64 / config.photon.n_photons as f64
    );
    println!(
        "{:<18} {} ({:.1}%)",
        "Escaped side:",
        results.n_escaped_side,
        100.0 * results.n_escaped_side as f64 / config.photon.n_photons as f64
    );
    println!(
        "{:<18} {:.2}",
        "Mean scatters:",
        results.total_scatters as f64 / config.photon.n_photons as f64
    );
    println!("{:<18} {:.4}", "Transmittance T:", results.transmittance);

    Ok(())
}
