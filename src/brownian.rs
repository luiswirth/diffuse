use crate::Vector;
use crate::path::{PathIncrements, TimeGrid};

use rand::RngExt;
use rand_distr::StandardNormal;

pub fn standard_normal_vector(dim: usize, rng: &mut impl RngExt) -> Vector {
  Vector::from_fn(dim, |_, _| rng.sample::<f64, _>(StandardNormal))
}

pub fn sample_brownian(dim: usize, grid: &TimeGrid, rng: &mut impl RngExt) -> PathIncrements {
  let increments = grid
    .dif()
    .values()
    .iter()
    .map(|dt| standard_normal_vector(dim, rng) * dt.sqrt())
    .collect();
  PathIncrements::new(increments)
}
