use crate::Vector;
use crate::path::{Cochain, PathIncrements, TimeGrid};

use rand::RngExt;
use rand_distr::StandardNormal;

/// The law of a Wiener integral over one cell, conditioned on the recorded
/// increment over that same cell: the projection of the weight onto the
/// recorded direction, and the norm of what is left orthogonal to it.
pub struct CellWeight {
  pub average: f64,
  pub variation: f64,
}

pub type CellWeights = Cochain<1, CellWeight>;

impl CellWeight {
  pub fn new(average: f64, variation: f64) -> Self {
    Self { average, variation }
  }

  /// The weight given by its integral and square integral over a cell of length `dt`.
  pub fn from_moments(integral: f64, square_integral: f64, dt: f64) -> Self {
    let average = integral / dt;
    let variation = (square_integral - integral * average).sqrt();
    Self::new(average, variation)
  }
}

pub fn standard_normal_vector(dim: usize, rng: &mut impl RngExt) -> Vector {
  Vector::from_fn(dim, |_, _| rng.sample::<f64, _>(StandardNormal))
}

pub fn wiener_integral(weight: &CellWeight, record: &Vector, rng: &mut impl RngExt) -> Vector {
  weight.average * record + weight.variation * standard_normal_vector(record.len(), rng)
}

/// Samples one Wiener integral per cell against `weights`, conditioned on `record`.
///
/// The sampled parts are independent of `record` alone, so a second call against
/// the same record yields integrals uncorrelated with these, which is wrong
/// unless the two weights are orthogonal.
pub fn wiener_integrals(
  weights: &CellWeights,
  record: &PathIncrements,
  rng: &mut impl RngExt,
) -> PathIncrements {
  assert_eq!(weights.nsteps(), record.nsteps());
  let integrals = weights
    .values()
    .iter()
    .zip(record.values())
    .map(|(weight, increment)| wiener_integral(weight, increment, rng))
    .collect();
  PathIncrements::new(integrals)
}

pub fn sample_brownian(dim: usize, grid: &TimeGrid, rng: &mut impl RngExt) -> PathIncrements {
  let weights = CellWeights::new(
    grid
      .dif()
      .values()
      .iter()
      .map(|dt| CellWeight::new(0.0, dt.sqrt()))
      .collect(),
  );
  let record = PathIncrements::new(vec![Vector::zeros(dim); grid.nsteps()]);
  wiener_integrals(&weights, &record, rng)
}
