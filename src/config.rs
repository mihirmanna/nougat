use crate::random::{henyey_greenstein, power_law, random_direction};
use anyhow::{Context, Result, anyhow};
use glam::DVec3 as Vec3;
use rand::{Rng, RngExt};
use rand_distr::Distribution;
use rand_distr::Normal;
use serde::Deserialize;
use std::fs;
use std::path::Path;
use toml;

#[derive(Deserialize, Debug)]
pub struct Config {
    pub(crate) simulation: SimulationConfig,
    pub photon: PhotonConfig,
    pub(crate) slab: SlabConfig,
    pub(crate) grid: GridConfig,
    pub(crate) clump: ClumpConfig,
}

impl Config {
    /// Loads and validates a simulation configuration from a TOML file.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be read
    /// - The TOML contents cannot be deserialized
    /// - Any configuration parameter fails validation
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let contents = fs::read_to_string(path)
            .with_context(|| format!("Could not read config file: {}", path.display()))?;
        let config: Config = toml::from_str(&contents)
            .with_context(|| format!("Invalid config file syntax: {}", path.display()))?;

        config.validate()?;

        Ok(config)
    }

    /// Validates the simulation values specified in the config.
    fn validate(&self) -> Result<()> {
        self.clump
            .validate()
            .context("Invalid clump configuration")?;
        self.photon
            .validate()
            .context("Invalid photon configuration")?;
        self.slab.validate().context("Invalid slab configuration")?;
        self.grid.validate().context("Invalid grid configuration")?;

        Ok(())
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct SimulationConfig {
    pub(crate) seed: u64,
}

#[derive(Deserialize, Debug)]
pub struct PhotonConfig {
    pub n_photons: u32,
    pub(crate) phase_function: PhaseFunction,
}

impl PhotonConfig {
    /// Validates the photon values specified in the config.
    fn validate(&self) -> Result<()> {
        self.phase_function.validate()
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct SlabConfig {
    pub(crate) width: f64,
    pub(crate) height: f64,
    pub(crate) edge_margin: f64,
    pub(crate) boundary_condition: BoundaryCondition,
}

impl SlabConfig {
    /// Validates the slab values specified in the config.
    fn validate(&self) -> Result<()> {
        if self.width <= 0.0 || self.height <= 0.0 || self.edge_margin <= 0.0 {
            return Err(anyhow!("Slab dimensions must be positive"));
        }

        Ok(())
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct GridConfig {
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    pub(crate) nz: usize,
}

impl GridConfig {
    /// Validates the grid values specified in the config.
    fn validate(&self) -> Result<()> {
        if self.nx == 0 || self.ny == 0 || self.nz == 0 {
            return Err(anyhow!("Grid dimensions must be nonzero"));
        }

        Ok(())
    }
}

#[derive(Deserialize, Debug)]
pub(crate) struct ClumpConfig {
    pub(crate) n_clumps: usize,
    pub(crate) radius: ValueDistribution,
    pub(crate) density: ValueDistribution,
    pub(crate) opacity: ValueDistribution,
}

impl ClumpConfig {
    /// Validates the clump values specified in the config.
    fn validate(&self) -> Result<()> {
        self.radius.validate().context("Invalid clump radius")?;
        self.density.validate().context("Invalid clump density")?;
        self.opacity.validate().context("Invalid clump opacity")?;

        Ok(())
    }
}

/// Defines the distribution used when sampling scalar values.
#[derive(Deserialize, Debug)]
#[serde(tag = "distribution", rename_all = "snake_case")]
pub(crate) enum ValueDistribution {
    Constant { value: f64 },
    Uniform { min: f64, max: f64 },
    Normal { mean: f64, stdev: f64 },
    PowerLaw { min: f64, max: f64, alpha: f64 },
}

impl ValueDistribution {
    /// Samples a random value from the distribution specified in the config.
    pub(crate) fn sample(&self, rng: &mut impl Rng) -> f64 {
        match self {
            Self::Constant { value } => *value,
            Self::Uniform { min, max } => rng.random_range(*min..*max),
            Self::Normal { mean, stdev } => {
                Normal::new(*mean, *stdev).unwrap().sample(rng).max(0.0)
            }
            Self::PowerLaw { min, max, alpha } => power_law(*min, *max, *alpha, rng),
        }
    }

    /// Validates the distribution values specified in the config.
    fn validate(&self) -> Result<()> {
        match self {
            Self::Constant { .. } => Ok(()),
            Self::Uniform { min, max } => {
                if min >= max {
                    Err(anyhow!("Uniform distribution requires min < max"))
                } else {
                    Ok(())
                }
            }
            Self::Normal { stdev, .. } => {
                if *stdev <= 0.0 {
                    Err(anyhow!("Normal distribution requires stdev > 0"))
                } else {
                    Ok(())
                }
            }
            Self::PowerLaw { min, max, .. } => {
                if min >= max {
                    Err(anyhow!("Power law requires min < max"))
                } else {
                    Ok(())
                }
            }
        }
    }
}

/// Defines the angular distribution used when sampling photon scattering directions.
#[derive(Deserialize, Debug)]
#[serde(tag = "phase_function", rename_all = "snake_case")]
pub(crate) enum PhaseFunction {
    Isotropic,
    HenyeyGreenstein { g: f64 },
}

impl PhaseFunction {
    /// Samples a random vector relative to the provided `incident` vector.
    pub(crate) fn sample(&self, incident: Vec3, rng: &mut impl Rng) -> Vec3 {
        match self {
            Self::Isotropic => random_direction(rng),
            Self::HenyeyGreenstein { g } => henyey_greenstein(*g, incident, rng),
        }
    }

    /// Validates the phase function values specified in the config.
    fn validate(&self) -> Result<()> {
        match self {
            Self::Isotropic => Ok(()),
            Self::HenyeyGreenstein { g } => {
                if !(-1.0..=1.0).contains(g) {
                    Err(anyhow!(
                        "HG asymmetry parameter g should be in range [-1.0, 1.0]"
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }
}

/// Defines the photon traversal condition enforced at the lateral boundaries of the slab.
#[derive(Deserialize, Debug)]
#[serde(tag = "boundary_condition")]
pub(crate) enum BoundaryCondition {
    Open,
    Reflective,
    Periodic,
}
