use crate::config::PhaseFunction;
use crate::random::random_tau;
use crate::ray::{Ray, ray_clump_intersection, ray_voxel_iter, t_exit_slab, t_exit_xy, t_exit_z};
use crate::slab::{Slab, VoxelGrid};
use glam::DVec3 as Vec3;
use rand::Rng;
use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

/// Particle that can scatter through clumps.
pub(crate) struct Photon {
    ray: Ray,
    n_scatters: u32,
}

impl Photon {
    /// Creates a new photon.
    pub fn new(ray: Ray) -> Self {
        Self { ray, n_scatters: 0 }
    }

    /// Scatters the photon from the point at `origin` along the direction provided by `dir`.
    pub fn scatter(&mut self, origin: Vec3, dir: Vec3) {
        self.ray = Ray::new(origin, dir);
        self.n_scatters += 1;
    }
}

/// Describes the point where a photon enters or exits a clump and provides an index into the slab
/// for that clump.
enum ClumpEvent {
    Enter { si: usize, t: f64 },
    Exit { si: usize, t: f64 },
}

impl ClumpEvent {
    /// Returns the `t` value at which a clump intersection event occurs.
    pub fn t(&self) -> f64 {
        match self {
            ClumpEvent::Enter { si: _, t } => *t,
            ClumpEvent::Exit { si: _, t } => *t,
        }
    }
}

impl PartialEq for ClumpEvent {
    fn eq(&self, other: &Self) -> bool {
        self.t() == other.t()
    }
}

impl Eq for ClumpEvent {}

impl PartialOrd for ClumpEvent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ClumpEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        self.t().total_cmp(&other.t())
    }
}

/// Describes how a photon exits the slab, including the number of scatters.
pub(crate) enum PhotonEscape {
    Top(u32),
    Bottom(u32),
    Side(u32),
}

/// Determines the points at which a given photon intersects clumps in the slab. Intersection points
/// are clamped to the boundaries of the slab on all sides.
///
/// # Returns
///
/// Returns an iterator over unique clump events ordered by their respective `t` values.
fn clump_events_iter(
    photon: &Photon,
    slab: &Slab,
    grid: &VoxelGrid,
) -> impl Iterator<Item = ClumpEvent> {
    let mut voxel_iter = ray_voxel_iter(&photon.ray, grid);
    let t_exit_slab = t_exit_slab(&photon.ray, slab);
    let mut t_boundary = 0.0; // Distance to exit current voxel
    let mut voxels_remaining = true;

    // Events are discovered per voxel, but must be emitted in order of t
    let mut events = BinaryHeap::<Reverse<ClumpEvent>>::new();

    // A clump may occupy multiple voxels, but should only be tested once
    let mut seen = vec![false; slab.clumps.len()];

    std::iter::from_fn(move || {
        loop {
            let has_eligible_event = events
                .peek()
                .is_some_and(|Reverse(event)| event.t() <= t_boundary);

            // Emit any events that occur before the next voxel boundary
            if has_eligible_event || !voxels_remaining {
                return events.pop().map(|Reverse(event)| event);
            }

            let Some((vi, t_boundary_next)) = voxel_iter.next() else {
                voxels_remaining = false;
                continue;
            };

            t_boundary = t_boundary_next;

            for &si in &grid.cells[vi] {
                if seen[si] {
                    continue;
                }
                seen[si] = true;
                let clump = &slab.clumps[si];

                if let Some((t_enter, t_exit)) = ray_clump_intersection(&photon.ray, clump) {
                    let t_enter = t_enter.max(0.0); // Photon may originate inside a clump
                    let t_exit = t_exit.min(t_exit_slab);

                    // Ignore intersections that lie entirely outside the slab
                    if t_enter < t_exit {
                        events.push(Reverse(ClumpEvent::Enter { si, t: t_enter }));
                        events.push(Reverse(ClumpEvent::Exit { si, t: t_exit }));
                    }
                }
            }
        }
    })
}

/// Inverts the optical depth equation along a photon's trajectory to calculate the physical
/// distance `t` at which it reaches an optical depth of `target_tau`.
///
/// # Arguments
///
/// - `slab`: The slab containing clumps
/// - `events`: An iterator of clump intersection events along the photon's trajectory
/// - `target_tau`: The target optical depth for the photon to reach
///
/// # Returns
///
/// Returns `Some(t)` if the photon reaches `target_tau` before exiting the slab. Otherwise,
/// returns `None`.
fn invert_optical_depth(
    slab: &Slab,
    events: impl Iterator<Item = ClumpEvent>,
    target_tau: f64,
) -> Option<f64> {
    let mut tau_acc = 0.0;
    let mut t_prev = 0.0;
    let mut active_clumps: Vec<usize> = Vec::new();

    for event in events {
        let t_next = event.t();

        if !active_clumps.is_empty() {
            let mean_density = active_clumps
                .iter()
                .map(|&si| slab.clumps[si].density)
                .sum::<f64>()
                / active_clumps.len() as f64;
            let mean_opacity = active_clumps
                .iter()
                .map(|&si| slab.clumps[si].opacity)
                .sum::<f64>()
                / active_clumps.len() as f64;

            let dt = t_next - t_prev;
            let dtau = dt * mean_density * mean_opacity;

            if tau_acc + dtau > target_tau {
                // target_tau falls within this segment, so solve linearly
                let t_hit = t_prev + (target_tau - tau_acc) / (mean_density * mean_opacity);
                return Some(t_hit);
            }

            tau_acc += dtau;
        }

        t_prev = t_next;

        // Update active clumps
        match event {
            ClumpEvent::Enter { si, t: _ } => active_clumps.push(si),
            ClumpEvent::Exit { si, t: _ } => active_clumps.retain(|&s| s != si),
        }
    }

    // target_tau is not reached, so the photon escapes
    None
}

/// Propagates a photon until it reaches an assigned optical depth or exits the slab.
pub(crate) fn transport_photon(
    photon: &mut Photon,
    phase_function: &PhaseFunction,
    slab: &Slab,
    grid: &VoxelGrid,
    rng: &mut impl Rng,
) -> PhotonEscape {
    loop {
        let target_tau = random_tau(rng);
        let events = clump_events_iter(photon, slab, grid);

        match invert_optical_depth(slab, events, target_tau) {
            Some(t_hit) => {
                // Scatter the photon in a new direction
                let new_origin = photon.ray.at(t_hit);
                let new_dir = phase_function.sample(photon.ray.dir, rng);
                photon.scatter(new_origin, new_dir);
            }
            None => {
                let t_exit_xy = t_exit_xy(&photon.ray, slab);
                let t_exit_z = t_exit_z(&photon.ray, slab);

                return if t_exit_xy < t_exit_z {
                    PhotonEscape::Side(photon.n_scatters)
                } else if photon.ray.dir.z > 0.0 {
                    PhotonEscape::Top(photon.n_scatters)
                } else {
                    PhotonEscape::Bottom(photon.n_scatters)
                };
            }
        }
    }
}
