//! Particle samplers for Ito SDEs and for the probability flow ODE.
//!
//! One state, one stepper: an ensemble of particles carried by a time-dependent
//! velocity field, with an optional diffusion coefficient. Setting the diffusion
//! to zero turns the SDE stepper into the deterministic flow, so both samplers
//! are the same code path.

pub mod brownian;
pub mod euler_maruyama;
pub mod ito;
pub mod ou;
pub mod path;

use nalgebra as na;

pub type Vector = na::DVector<f64>;
pub type Matrix = na::DMatrix<f64>;

pub type State = Vector;
