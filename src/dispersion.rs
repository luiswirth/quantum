use crate::{
  exact::plane_wave::{HarmonicWave, PlaneWave},
  units::{ELECTRON_MASS, HBAR, LIGHT_SPEED},
};

/// Only the positive frequency branch, so the sign of the wavenumber alone
/// says which way the wave travels.
#[derive(Clone, Copy, Debug)]
pub enum Dispersion {
  /// `omega(k) = (planck k^2) / (2 m)`
  FreeParticle { mass: f64 },
  /// `omega(k) = c abs(k)`
  Light { speed: f64 },
}

impl Dispersion {
  pub fn electron() -> Self {
    Self::FreeParticle {
      mass: ELECTRON_MASS,
    }
  }

  pub fn vacuum_light() -> Self {
    Self::Light { speed: LIGHT_SPEED }
  }

  pub fn frequency(&self, wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR * wavenumber.powi(2) / (2.0 * mass),
      Self::Light { speed } => speed * wavenumber.abs(),
    }
  }

  /// `v_g = (dif omega) / (dif k)`, which light leaves undefined at `k = 0`.
  pub fn group_velocity(&self, wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR * wavenumber / mass,
      Self::Light { speed } => speed * wavenumber.signum(),
    }
  }

  /// `(dif^2 omega) / (dif k^2)`
  pub fn curvature(&self, _wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR / mass,
      Self::Light { .. } => 0.0,
    }
  }

  /// `m^* = planck / ((dif^2 omega) / (dif k^2))`
  pub fn effective_mass(&self, wavenumber: f64) -> f64 {
    HBAR / self.curvature(wavenumber)
  }

  /// `E(p) = planck omega(p / planck)`
  pub fn energy(&self, momentum: f64) -> f64 {
    HBAR * self.frequency(momentum / HBAR)
  }

  /// `v_p = omega / k`
  pub fn phase_velocity(&self, wavenumber: f64) -> f64 {
    self.frequency(wavenumber) / wavenumber
  }

  pub fn plane_wave(&self, wavenumber: f64) -> PlaneWave {
    PlaneWave::new(self.frequency(wavenumber), HarmonicWave::new(wavenumber))
  }
}
