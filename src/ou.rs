use crate::State;
use crate::ito::ItoProcess;
use crate::path::{PathIncrements, PathStates, TimeGrid};

pub struct OrnsteinUhlenbeck {
  pub theta: f64,
  pub sigma: f64,
}

impl Default for OrnsteinUhlenbeck {
  fn default() -> Self {
    Self::new(1.0, 2f64.sqrt())
  }
}

impl OrnsteinUhlenbeck {
  pub fn new(theta: f64, sigma: f64) -> Self {
    Self { theta, sigma }
  }
  pub fn exact_path(&self, x0: State, grid: &TimeGrid, brownian: &PathIncrements) -> PathStates {
    todo!()
  }
}

impl ItoProcess for OrnsteinUhlenbeck {
  type Diffusion = f64;

  fn drift(&self, _t: f64, x: &State) -> State {
    -self.theta * x
  }
  fn diffusion(&self, _t: f64, _x: &State) -> f64 {
    self.sigma
  }
}
