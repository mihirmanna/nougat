use glam::DQuat as Quat;
use glam::DVec3 as Vec3;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, UnitCircle, UnitSphere};

/// Samples a random optical depth from an exponential distribution.
pub(crate) fn random_tau(rng: &mut impl Rng) -> f64 {
    let x: f64 = rng.random_range(f64::MIN_POSITIVE..=1.0);

    // Optical depth tau = -ln(x)
    -x.ln()
}

/// Uniformly samples a random 3-vector from the unit sphere.
pub(crate) fn random_direction(rng: &mut impl Rng) -> Vec3 {
    let dir = UnitSphere.sample(rng);

    Vec3::new(dir[0], dir[1], dir[2])
}

/// Uniformly samples a random 3-vector from the upper half (z > 0) of the unit sphere.
pub(crate) fn random_upward(rng: &mut impl Rng) -> Vec3 {
    loop {
        let dir = random_direction(rng);
        if dir.z > 0.0 {
            return dir;
        }
    }
}

/// Samples the Henyey-Greenstein (HG) phase function, using the cumulative distribution function
/// derived by Witt (1977)[^1].
///
/// # Arguments
///
/// - `g`: HG asymmetry factor, between -1 and 1 (inclusive)
/// - `incident`: The direction vector of the incident photon
///
/// # Returns
///
/// Returns the direction vector of the outgoing (scattered) photon.
///
/// [^1]: <https://doi.org/10.1086/190463>
pub(crate) fn henyey_greenstein(g: f64, incident: Vec3, rng: &mut impl Rng) -> Vec3 {
    let x: f64 = rng.random();

    let cos_theta = if g.abs() < 1e-15 {
        // g = 0 is isotropic, so sample uniformly
        rng.random_range(-1.0..1.0)
    } else {
        1.0 / (2.0 * g) * (1.0 + g * g - ((1.0 - g * g) / (1.0 - g + 2.0 * g * x)).powi(2))
    };
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let [cos_phi, sin_phi]: [f64; 2] = UnitCircle.sample(rng);

    // Scattered direction in local frame (incident along +z)
    let local = Vec3::new(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

    Quat::from_rotation_arc(Vec3::Z, incident.normalize()) * local
}

/// Samples a random value between `min` and `max` from a power-law distribution with exponent
/// `alpha`.
pub(crate) fn power_law(min: f64, max: f64, alpha: f64, rng: &mut impl Rng) -> f64 {
    let x: f64 = rng.random();
    let beta = alpha + 1.0;

    if beta.abs() < 1e-15 {
        // alpha = -1 => Reciprocal distribution
        min * (max / min).powf(x)
    } else {
        let x_min_b = min.powf(beta);
        let x_max_b = max.powf(beta);

        ((x_max_b - x_min_b) * x + x_min_b).powf(1.0 / beta)
    }
}
