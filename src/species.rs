use crate::{plane_wave::PlaneWave, units::HBAR};

/// A kind of free non-relativistic particle, which is what ties a wavenumber
/// to a frequency: it supplies the dispersion relation `omega(k)`.
#[derive(Clone, Copy, Debug)]
pub struct Species {
  pub mass: f64,
}

impl Species {
  pub fn new(mass: f64) -> Self {
    Self { mass }
  }

  pub fn electron() -> Self {
    Self::new(crate::units::ELECTRON_MASS)
  }

  /// The dispersion relation, the kinetic energy `p^2 / 2 m` written in the
  /// wavenumber and divided by `hbar`.
  pub fn frequency(&self, wavenumber: f64) -> f64 {
    HBAR * wavenumber.powi(2) / (2.0 * self.mass)
  }

  /// The speed the envelope of a packet centered on this wavenumber travels
  /// at, `d omega / d k`, and the classical velocity `p / m`.
  pub fn group_velocity(&self, wavenumber: f64) -> f64 {
    HBAR * wavenumber / self.mass
  }

  pub fn plane_wave(&self, wavenumber: f64) -> PlaneWave {
    PlaneWave::new(wavenumber, self.frequency(wavenumber))
  }
}

#[cfg(test)]
mod tests {
  use super::Species;
  use crate::units::{ELECTRON_MASS, HBAR};

  #[test]
  fn group_velocity_is_the_slope_of_the_dispersion_relation() {
    let species = Species::electron();
    let step = 1e-6;
    for index in 1..20 {
      let wavenumber = index as f64 * 0.5;
      let slope = (species.frequency(wavenumber + step) - species.frequency(wavenumber - step))
        / (2.0 * step);
      assert!((slope - species.group_velocity(wavenumber)).abs() < 1e-6);
    }
  }

  #[test]
  fn the_phase_moves_at_half_the_particle() {
    let species = Species::electron();
    for index in 1..20 {
      let wavenumber = index as f64 * 0.5;
      let phase_velocity = species.plane_wave(wavenumber).phase_velocity();
      assert!((2.0 * phase_velocity - species.group_velocity(wavenumber)).abs() < 1e-12);
    }
  }

  /// The dispersion relation says that the wave's energy is the particle's
  /// kinetic energy.
  #[test]
  fn energy_is_kinetic() {
    let species = Species::electron();
    for index in 1..20 {
      let wavenumber = index as f64 * 0.5;
      let wave = species.plane_wave(wavenumber);
      let kinetic = wave.momentum().powi(2) / (2.0 * species.mass);
      assert!((wave.energy() - kinetic).abs() < 1e-12);
    }
  }

  /// The combination every effective mass model is quoted in.
  #[test]
  fn electron_energy_scale() {
    let scale = HBAR.powi(2) / (2.0 * ELECTRON_MASS);
    assert!((scale - 0.0380998).abs() < 1e-7);
  }
}
