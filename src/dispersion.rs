use crate::{
  plane_wave::PlaneWave,
  units::{ELECTRON_MASS, HBAR, LIGHT_SPEED},
};

/// The relation `omega(k)` tying a wave's frequency to its wavenumber, which
/// is where the physics of the medium enters and the only place it does as
/// long as the medium is homogeneous.
///
/// Only the positive frequency branch is carried, so the sign of the
/// wavenumber alone says which way the wave travels.
#[derive(Clone, Copy, Debug)]
pub enum Dispersion {
  /// `omega = hbar k^2 / 2 m`
  FreeParticle { mass: f64 },
  /// `omega = c |k|`, with `c` the speed in the medium
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

  /// The speed a packet centered on this wavenumber travels at, `d omega / d k`.
  ///
  /// The light branch has a kink at rest, where the two directions disagree
  /// and the derivative does not exist.
  pub fn group_velocity(&self, wavenumber: f64) -> f64 {
    match *self {
      Self::FreeParticle { mass } => HBAR * wavenumber / mass,
      Self::Light { speed } => speed * wavenumber.signum(),
    }
  }

  /// The speed the phase travels at, `omega / k`.
  pub fn phase_velocity(&self, wavenumber: f64) -> f64 {
    self.frequency(wavenumber) / wavenumber
  }

  pub fn plane_wave(&self, wavenumber: f64) -> PlaneWave {
    PlaneWave::new(wavenumber, self.frequency(wavenumber))
  }
}
