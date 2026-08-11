use crate::ito::ItoProcess;
use crate::path::{PathIncrements, PathStates, TimeGrid};
use crate::wiener::{CellWeight, CellWeights};
use crate::{Matrix, State};

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
  /// The weights of the stochastic integrals `exact_path` consumes,
  /// one per cell, from `f(s) = exp(-theta (t_(i+1) - s))`.
  pub fn weights(&self, grid: &TimeGrid) -> CellWeights {
    let weights = grid
      .dif()
      .values()
      .iter()
      .map(|dt| {
        let integral = -(-self.theta * dt).exp_m1() / self.theta;
        let square_integral = -(-2.0 * self.theta * dt).exp_m1() / (2.0 * self.theta);
        CellWeight::from_moments(integral, square_integral, *dt)
      })
      .collect();
    CellWeights::new(weights)
  }

  pub fn exact_path(&self, x0: State, grid: &TimeGrid, noise: &PathIncrements) -> PathStates {
    assert_eq!(grid.nsteps(), noise.nsteps());
    let dts = grid.dif();

    let mut states = Vec::with_capacity(grid.nsteps() + 1);
    states.push(x0);

    for istep in 0..grid.nsteps() {
      let x = &states[istep];
      let dt = dts.values()[istep];
      let dw = &noise.values()[istep];
      states.push((-self.theta * dt).exp() * x + self.sigma * dw);
    }

    PathStates::new(states)
  }
}

impl ItoProcess for OrnsteinUhlenbeck {
  fn drift(&self, _t: f64, x: &State) -> State {
    -self.theta * x
  }
  fn diffusion(&self, _t: f64, x: &State) -> Matrix {
    self.sigma * Matrix::identity(x.len(), x.len())
  }
}
