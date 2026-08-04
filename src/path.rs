use crate::{State, Vector};
use std::ops::{Add, Sub};

pub struct Cochain<const DEG: usize, V> {
  values: Vec<V>,
}

pub type TimeGrid = Cochain<0, f64>;
pub type TimeSteps = Cochain<1, f64>;
pub type PathStates = Cochain<0, State>;
pub type PathIncrements = Cochain<1, Vector>;

impl<const DEG: usize, V> Cochain<DEG, V> {
  pub fn new(values: Vec<V>) -> Self {
    Self { values }
  }
  pub fn values(&self) -> &[V] {
    &self.values
  }
}

impl<V> Cochain<0, V> {
  pub fn nsteps(&self) -> usize {
    self.values.len().saturating_sub(1)
  }
}

impl<V> Cochain<1, V> {
  pub fn nsteps(&self) -> usize {
    self.values.len()
  }
}

impl<V: Clone + Sub<Output = V>> Cochain<0, V> {
  pub fn dif(&self) -> Cochain<1, V> {
    let values = self
      .values
      .windows(2)
      .map(|pair| pair[1].clone() - pair[0].clone())
      .collect();
    Cochain::new(values)
  }
}

impl<V: Clone + Add<Output = V>> Cochain<1, V> {
  pub fn integrate(&self, initial: V) -> Cochain<0, V> {
    let mut values = Vec::with_capacity(self.values.len() + 1);
    let mut acc = initial;
    values.push(acc.clone());
    for increment in &self.values {
      acc = acc + increment.clone();
      values.push(acc.clone());
    }
    Cochain::new(values)
  }
}

impl<V: Clone> Cochain<0, V> {
  /// The de Rham map.
  pub fn coarsen(&self, factor: usize) -> Self {
    assert!(factor > 0);
    assert!(self.nsteps().is_multiple_of(factor));
    Self::new(self.values.iter().step_by(factor).cloned().collect())
  }
}

impl<V: Clone + Add<Output = V>> Cochain<1, V> {
  /// The de Rham map.
  pub fn coarsen(&self, factor: usize) -> Self {
    assert!(factor > 0);
    assert!(self.nsteps().is_multiple_of(factor));
    let values = self
      .values
      .chunks(factor)
      .map(|block| block.iter().cloned().reduce(|sum, dv| sum + dv).unwrap())
      .collect();
    Self::new(values)
  }
}

impl TimeGrid {
  pub fn uniform(nsteps: usize, dt: f64) -> Self {
    Self::new((0..=nsteps).map(|istep| istep as f64 * dt).collect())
  }
}
