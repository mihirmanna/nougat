# Nougat

Nougat is a Monte Carlo radiative transfer code for simulating photon propagation through a plane-parallel slab populated by discrete, spherical clumps. 
It implements the algorithm described in [Townsend (2007)](https://arxiv.org/abs/0709.0860).

## Features

- Parallelized photon transport with voxel-based spatial partitioning for fast lookup
- Isotropic and Henyey-Greenstein scattering phase functions
- Configurable distributions for clump radius, density, and opacity:
  - Constant
  - Uniform
  - Normal (Gaussian)
  - Power law
- Human-readable TOML configuration files with parameter validation

## Building

Compile with Cargo:

```bash
cargo build --release
```

The executable will be located at `target/release/nougat`.

## Usage

Simulation parameters are specified in a TOML configuration file. 
An annotated example containing all available options is provided in [`examples/example.toml`](examples/example_config.toml).

Run with Cargo:

```bash
cargo run --release -- <CONFIG.toml>
```

Or run the compiled binary directly:

```bash
./nougat <CONFIG.toml>
```

## Example output

```
Photons:           1000000
Escaped top:       470340 (47.0%)
Escaped bottom:    448550 (44.9%)
Escaped side:      81110 (8.1%)
Mean scatters:     7.30
Transmittance T:   0.4703
```
