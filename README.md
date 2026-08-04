# diffuse

Particle samplers for Ito SDEs and for the probability flow ODE, in 2D, with PNG output.

The point is to make the theory behind diffusion models computable at the smallest scale
where it is still the real thing:
a velocity field, an ensemble of particles, and a score.
There is no neural network here.
Where a score is needed it is available in closed form,
which leaves the time discretization as the only approximation
and makes the sampler error measurable rather than merely bounded.

## Plan

1. Euler-Maruyama on the Ornstein-Uhlenbeck process.
   OU has both an exact solution and an exact stationary distribution,
   so the empirical histogram and the measured convergence order are two independent checks.
   Milstein alongside it, to see the strong order rise from $1/2$ to $1$.
2. The reverse-time SDE on a Gaussian mixture.
   The forward noising of a mixture is again a mixture,
   so $grad log p_t$ is closed form.
   Integrating the reverse SDE from noise recovers the mixture.
3. The probability flow ODE.
   Same score, deterministic velocity
   $v = b - 1/2 sigma sigma^top grad log p_t$,
   run side by side with the reverse SDE on the same seed.
   The ODE trajectories are streamlines and do not cross; the SDE ones do.

## Building

Rust comes from the flake:

    nix develop
    cargo test
