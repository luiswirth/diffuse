use crate::ito::ItoProcess;
use crate::path::{PathIncrements, PathStates, TimeGrid};
use crate::{State, Vector};
use std::ops::Mul;

pub fn euler_maruyama_step<Diffusion>(
  drift: Vector,
  diffusion: Diffusion,
  mut x: State,
  dt: f64,
  dw: &Vector,
) -> State
where
  Diffusion: for<'a> Mul<&'a Vector, Output = Vector>,
{
  x += drift * dt + diffusion * dw;
  x
}

pub fn euler_maruyama(
  process: &impl ItoProcess,
  x0: State,
  grid: &TimeGrid,
  brownian: &PathIncrements,
) -> PathStates {
  assert_eq!(grid.nsteps(), brownian.nsteps());
  let dts = grid.dif();

  let mut states = Vec::with_capacity(grid.nsteps() + 1);
  states.push(x0);

  for istep in 0..grid.nsteps() {
    let x = &states[istep];
    let t = grid.values()[istep];
    let dt = dts.values()[istep];
    let dw = &brownian.values()[istep];
    let drift = process.drift(t, x);
    let diffusion = process.diffusion(t, x);
    states.push(euler_maruyama_step(drift, diffusion, x.clone(), dt, dw));
  }

  PathStates::new(states)
}
