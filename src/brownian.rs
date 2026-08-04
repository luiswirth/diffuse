//! The driving noise of an SDE, as increments on a uniform grid.

use crate::Vector;

use rand::RngExt;
use rand_distr::StandardNormal;

pub fn standard_normal_vector(n: usize, rng: &mut impl RngExt) -> Vector {
  Vector::from_fn(n, |_, _| rng.sample::<f64, _>(StandardNormal))
}

/// One realization of a Brownian motion in R^m,
/// held as its increments over a uniform grid of width `dt`.
/// Each increment is N(0, dt I).
pub struct BrownianPath {
  dt: f64,
  increments: Vec<Vector>,
}

impl BrownianPath {
  pub fn sample(dim: usize, nsteps: usize, dt: f64, rng: &mut impl RngExt) -> Self {
    let increments = (0..nsteps)
      .map(|_| standard_normal_vector(dim, rng) * dt.sqrt())
      .collect();
    Self { dt, increments }
  }

  pub fn dt(&self) -> f64 {
    self.dt
  }
  pub fn nsteps(&self) -> usize {
    self.increments.len()
  }
  pub fn dim(&self) -> usize {
    self.increments.first().map_or(0, Vector::len)
  }
  pub fn increments(&self) -> &[Vector] {
    &self.increments
  }

  pub fn coarsen(&self, factor: usize) -> Self {
    assert!(factor > 0);
    assert!(self.nsteps().is_multiple_of(factor));
    let dim = self.dim();
    let increments = self
      .increments
      .chunks(factor)
      .map(|block| block.iter().fold(Vector::zeros(dim), |sum, dw| sum + dw))
      .collect();
    Self {
      dt: self.dt * factor as f64,
      increments,
    }
  }
}
