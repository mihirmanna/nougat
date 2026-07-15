use glam::{DVec3 as Vec3, IVec3};

/// Spherical, homogeneous clump
pub(crate) struct Clump {
    pub(crate) pos: Vec3,
    pub(crate) radius: f64,
    pub(crate) density: f64,
    pub(crate) opacity: f64,
}

/// Plane-parallel slab containing clumps
pub(crate) struct Slab {
    pub(crate) width: f64,
    pub(crate) height: f64,
    pub(crate) clumps: Vec<Clump>,
}

impl Slab {
    /// Creates a slab of size `width` in the x- and y-directions, and `height` in the z-direction.
    pub(crate) fn new(width: f64, height: f64) -> Self {
        Self {
            width,
            height,
            clumps: Vec::new(),
        }
    }

    /// Adds a single clump to the slab.
    pub fn add_clump(&mut self, clump: Clump) {
        self.clumps.push(clump);
    }
}

/// Voxel structure implemented over a slab to optimize ray-clump intersection tests
pub(crate) struct VoxelGrid {
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    pub(crate) nz: usize,
    pub(crate) dx: f64,
    pub(crate) dy: f64,
    pub(crate) dz: f64,
    pub(crate) cells: Vec<Vec<usize>>,
}

impl VoxelGrid {
    /// Overlays a voxel grid on the slab with `nx` cells in the x-direction, `ny` cells in the
    /// y-direction, and `nz` cells in the z-direction.
    pub(crate) fn new(slab: &Slab, nx: usize, ny: usize, nz: usize) -> Self {
        let dx = slab.width / nx as f64;
        let dy = slab.width / ny as f64;
        let dz = slab.height / nz as f64;

        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); nx * ny * nz];

        for (si, clump) in slab.clumps.iter().enumerate() {
            // AABB bounding box of clump in voxel coordinates
            let lower_corner = ((clump.pos - clump.radius) / Vec3::new(dx, dy, dz))
                .floor()
                .as_ivec3();
            let upper_corner = ((clump.pos + clump.radius) / Vec3::new(dx, dy, dz))
                .floor()
                .as_ivec3();

            // Clamp bounding box to grid dimensions
            let x_range = lower_corner.x.max(0)..upper_corner.x.min(nx as i32);
            let y_range = lower_corner.y.max(0)..upper_corner.y.min(ny as i32);
            let z_range = lower_corner.z.max(0)..upper_corner.z.min(nz as i32);

            for iz in z_range {
                for iy in y_range.clone() {
                    for ix in x_range.clone() {
                        cells[ix as usize + nx * (iy as usize + ny * iz as usize)].push(si);
                    }
                }
            }
        }

        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            cells,
        }
    }

    /// Flattens a 3-dimensional grid index (`ix`, `iy`, `iz`) into a 1-dimensional index.
    pub(crate) fn flat_idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }

    /// Returns `true` if the given voxel index lies inside the grid.
    pub(crate) fn contains(&self, voxel: IVec3) -> bool {
        let dims = IVec3::new(self.nx as i32, self.ny as i32, self.nz as i32);

        (voxel.cmpge(IVec3::ZERO) & voxel.cmplt(dims)).all()
    }
}
