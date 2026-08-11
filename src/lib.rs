//! Samplers for Ito SDEs, with grids and paths as cochains on the time grid.

pub mod euler_maruyama;
pub mod ito;
pub mod ou;
pub mod path;
pub mod wiener;

use nalgebra as na;

pub type Vector = na::DVector<f64>;
pub type Matrix = na::DMatrix<f64>;

pub type State = Vector;
