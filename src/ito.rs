use crate::{Matrix, State};

// TODO: dense diffusion wastes isotropy
pub trait ItoProcess {
  fn drift(&self, t: f64, x: &State) -> State;
  fn diffusion(&self, t: f64, x: &State) -> Matrix;
}
