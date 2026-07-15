use crate::slab::{Clump, Slab, VoxelGrid};
use glam::{DVec3 as Vec3, IVec3};

/// 3-dimensional ray
pub(crate) struct Ray {
    origin: Vec3,
    pub(crate) dir: Vec3,
}

impl Ray {
    /// Creates a ray and normalizes its direction vector.
    pub(crate) fn new(origin: Vec3, dir: Vec3) -> Self {
        Self {
            origin,
            dir: dir.normalize(),
        }
    }

    /// Returns the position of a point at distance `t` along the ray.
    pub(crate) fn at(&self, t: f64) -> Vec3 {
        self.origin + self.dir * t
    }
}

/// Calculates the time for a ray to exit the slab in the x or y direction.
pub(crate) fn t_exit_xy(ray: &Ray, slab: &Slab) -> f64 {
    let t_x = if ray.dir.x > 0.0 {
        (slab.width - ray.origin.x) / ray.dir.x
    } else if ray.dir.x < 0.0 {
        -ray.origin.x / ray.dir.x
    } else {
        f64::INFINITY
    };

    let t_y = if ray.dir.y > 0.0 {
        (slab.width - ray.origin.y) / ray.dir.y
    } else if ray.dir.y < 0.0 {
        -ray.origin.y / ray.dir.y
    } else {
        f64::INFINITY
    };

    t_x.min(t_y)
}

/// Calculates the time for a ray to exit the slab in the z direction.
pub(crate) fn t_exit_z(ray: &Ray, slab: &Slab) -> f64 {
    if ray.dir.z >= 0.0 {
        (slab.height - ray.origin.z) / ray.dir.z
    } else {
        -ray.origin.z / ray.dir.z
    }
}

/// Calculates the time for a ray to exit the slab in any direction.
pub(crate) fn t_exit_slab(ray: &Ray, slab: &Slab) -> f64 {
    t_exit_xy(ray, slab).min(t_exit_z(ray, slab))
}

/// Calculates the intersection points of a ray with a spherical clump.
///
/// # Returns
///
/// Returns `None` if:
/// - The clump is located behind the ray origin
/// - The ray does not pass through the clump (forwards or backwards)
/// - The ray only glances against the clump at one point
///
/// Otherwise, returns the two `t` values along the ray where it enters and exits the clump.
pub(crate) fn ray_clump_intersection(ray: &Ray, clump: &Clump) -> Option<(f64, f64)> {
    let v = ray.origin - clump.pos;
    let b = 2.0 * ray.dir.dot(v);
    let c = v.length_squared() - clump.radius * clump.radius;
    let discriminant = b * b - 4.0 * c;

    // Ray misses or glances against sphere
    if discriminant <= 0.0 {
        return None;
    }

    let sqrt_d = discriminant.sqrt();
    let t1 = (-b - sqrt_d) / 2.0;
    let t2 = (-b + sqrt_d) / 2.0;

    // Check if sphere is behind ray origin before returning
    if t2 < 0.0 { None } else { Some((t1, t2)) }
}

/// Calculates the voxels that a ray passes through before exiting the top or bottom of the slab,
/// using the algorithm developed by Amanatides and Woo (1987)[^1].
///
/// # Returns
///
/// Returns an iterator of tuples that each contain the index of the intersected voxel and the value
/// of `t` along the ray at the voxel boundary.
///
/// [^1]: <https://doi.org/10.2312/egtp.19871000>
pub(crate) fn ray_voxel_iter(ray: &Ray, grid: &VoxelGrid) -> impl Iterator<Item = (usize, f64)> {
    // Initialization phase
    let voxel_size = Vec3::new(grid.dx, grid.dy, grid.dz);
    let dir_step = ray.dir.signum().as_ivec3();
    let t_delta = (voxel_size / ray.dir).abs();
    let mut curr_voxel = (ray.origin / voxel_size).floor().as_ivec3();
    let next_boundary = (curr_voxel + dir_step.max(IVec3::ZERO)).as_dvec3() * voxel_size;
    let mut t_max = (next_boundary - ray.origin) / ray.dir;

    // Iterative traversal phase
    std::iter::from_fn(move || {
        if !grid.contains(curr_voxel) {
            return None; // Ray has exited the slab
        }

        let t_boundary = t_max.min_element(); // t value of current voxel boundary
        let result = (
            grid.flat_idx(
                curr_voxel.x as usize,
                curr_voxel.y as usize,
                curr_voxel.z as usize,
            ),
            t_boundary,
        );

        // Advance to next voxel
        let step_axis = t_max.min_position();
        curr_voxel[step_axis] += dir_step[step_axis];
        t_max[step_axis] += t_delta[step_axis];

        Some(result)
    })
}
