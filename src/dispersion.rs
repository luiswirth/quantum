use crate::{
  plane_wave::PlaneWave,
  units::{ELECTRON_MASS, HBAR, LIGHT_SPEED},
};

/// Only the positive frequency branch, so the sign of the wavenumber alone
/// says which way the wave travels.
#[derive(Clone, Copy, Debug)]
pub enum Dispersion {
  FreeParticle { mass: f64 },
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

  /// Light has a kink at rest, where the derivative does not exist.
  pub fn group_velocity(&self, wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR * wavenumber / mass,
      Self::Light { speed } => speed * wavenumber.signum(),
    }
  }

  /// The second derivative of the dispersion relation.
  pub fn curvature(&self, _wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR / mass,
      Self::Light { .. } => 0.0,
    }
  }

  pub fn effective_mass(&self, wavenumber: f64) -> f64 {
    HBAR / self.curvature(wavenumber)
  }

  pub fn energy(&self, momentum: f64) -> f64 {
    HBAR * self.frequency(momentum / HBAR)
  }

  pub fn phase_velocity(&self, wavenumber: f64) -> f64 {
    self.frequency(wavenumber) / wavenumber
  }

  pub fn plane_wave(&self, wavenumber: f64) -> PlaneWave {
    PlaneWave::new(wavenumber, self.frequency(wavenumber))
  }
}
