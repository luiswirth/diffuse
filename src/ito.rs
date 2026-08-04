use crate::{State, Vector};
use std::ops::Mul;

pub trait ItoProcess {
  type Diffusion: for<'a> Mul<&'a Vector, Output = Vector>;

  fn drift(&self, t: f64, x: &State) -> State;
  fn diffusion(&self, t: f64, x: &State) -> Self::Diffusion;
}
