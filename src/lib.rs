pub mod config;
mod photon;
mod random;
mod ray;
mod slab;

use config::Config;
use glam::DVec3 as Vec3;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use photon::{Photon, PhotonEscape, transport_photon};
use rand::prelude::StdRng;
use rand::{RngExt, SeedableRng};
use random::random_upward;
use ray::Ray;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use slab::{Clump, Slab, VoxelGrid};
use std::sync::atomic::{AtomicU32, Ordering};

/// Stores the result of a simulation run for display and calculations.
pub struct SimulationResult {
    pub n_escaped_top: u32,
    pub n_escaped_bottom: u32,
    pub n_escaped_side: u32,
    pub total_scatters: u32,
    pub transmittance: f64,
}

/// Runs the simulation with parameters from the config file.
pub fn run(config: &mut Config) -> SimulationResult {
    let mut rng = StdRng::seed_from_u64(config.simulation.seed);

    let mut slab = Slab::new(config.slab.width, config.slab.height);
    for _ in 0..config.clump.n_clumps {
        slab.add_clump(Clump {
            pos: Vec3::new(
                rng.random_range(0.0..slab.width),
                rng.random_range(0.0..slab.width),
                rng.random_range(0.0..slab.height),
            ),
            radius: config.clump.radius.sample(&mut rng),
            density: config.clump.density.sample(&mut rng),
            opacity: config.clump.opacity.sample(&mut rng),
        })
    }

    let grid = VoxelGrid::new(&slab, config.grid.nx, config.grid.ny, config.grid.nz);

    let n_photons = config.photon.n_photons;
    let phase_function = &config.photon.phase_function;
    let n_escaped_top = AtomicU32::new(0);
    let n_escaped_bottom = AtomicU32::new(0);
    let n_escaped_side = AtomicU32::new(0);
    let total_scatters = AtomicU32::new(0);

    let pb = ProgressBar::new(n_photons as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{elapsed_precise} [{bar:40}] {pos}/{len} ({eta})")
            .unwrap(),
    );

    (0..n_photons)
        .into_par_iter()
        .progress_with(pb.clone())
        .for_each(|i| {
            let mut rng = StdRng::seed_from_u64(config.simulation.seed + i as u64);

            // Inject photon at random position on the bottom of the slab, inside the margins
            let origin = Vec3::new(
                rng.random_range(config.slab.edge_margin..(slab.width - config.slab.edge_margin)),
                rng.random_range(config.slab.edge_margin..(slab.width - config.slab.edge_margin)),
                0.0,
            );
            let dir = random_upward(&mut rng);
            let mut photon = Photon::new(Ray::new(origin, dir));

            match transport_photon(&mut photon, phase_function, &slab, &grid, &mut rng) {
                PhotonEscape::Top(n) => {
                    n_escaped_top.fetch_add(1, Ordering::Relaxed);
                    total_scatters.fetch_add(n, Ordering::Relaxed);
                }
                PhotonEscape::Bottom(n) => {
                    n_escaped_bottom.fetch_add(1, Ordering::Relaxed);
                    total_scatters.fetch_add(n, Ordering::Relaxed);
                }
                PhotonEscape::Side(n) => {
                    n_escaped_side.fetch_add(1, Ordering::Relaxed);
                    total_scatters.fetch_add(n, Ordering::Relaxed);
                }
            }
        });

    let n_escaped_top = n_escaped_top.load(Ordering::Relaxed);
    let n_escaped_bottom = n_escaped_bottom.load(Ordering::Relaxed);
    let n_truncated = n_escaped_side.load(Ordering::Relaxed);
    let total_scatters = total_scatters.load(Ordering::Relaxed);

    SimulationResult {
        n_escaped_top,
        n_escaped_bottom,
        n_escaped_side: n_truncated,
        total_scatters,
        transmittance: n_escaped_top as f64 / n_photons as f64,
    }
}
