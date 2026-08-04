//! The Euler-Maruyama scheme.

use crate::brownian::BrownianPath;
use crate::{Matrix, Vector};

pub fn euler_maruyama_step(
  drift: Vector,
  diffusion: Matrix,
  mut x: Vector,
  dt: f64,
  dw: &Vector,
) -> Vector {
  x += drift * dt + diffusion * dw;
  x
}

pub fn euler_maruyama<DriftFn, DiffusionFn>(
  drift: DriftFn,
  diffusion: DiffusionFn,
  x: Vector,
  path: &BrownianPath,
) -> Vec<Vector>
where
  DriftFn: Fn(f64, &Vector) -> Vector,
  DiffusionFn: Fn(f64, &Vector) -> Matrix,
{
  let dt = path.dt();
  let mut xs = Vec::with_capacity(path.nsteps() + 1);
  xs.push(x);

  for (istep, dw) in path.increments().iter().enumerate() {
    let x = &xs[istep];
    let t = istep as f64 * dt;
    let drift = drift(t, x);
    let diffusion = diffusion(t, x);
    let x = euler_maruyama_step(drift, diffusion, x.clone(), dt, dw);
    xs.push(x);
  }

  xs
}
