use crate::{Complex, units::HBAR};
use std::f64::consts::TAU;

/// The wave `exp(i (k x - omega t))`, of unit modulus everywhere.
///
/// The sign of the exponent is the convention the whole crate is written in:
/// a positive wavenumber moves the wave towards positive `x`.
#[derive(Clone, Copy, Debug)]
pub struct PlaneWave {
  pub wavenumber: f64,
  pub frequency: f64,
}

impl PlaneWave {
  pub fn new(wavenumber: f64, frequency: f64) -> Self {
    Self {
      wavenumber,
      frequency,
    }
  }

  pub fn phase(&self, position: f64, time: f64) -> f64 {
    self.wavenumber * position - self.frequency * time
  }
  pub fn at(&self, position: f64, time: f64) -> Complex {
    Complex::from_polar(1.0, self.phase(position, time))
  }

  /// de Broglie.
  pub fn momentum(&self) -> f64 {
    HBAR * self.wavenumber
  }
  /// Planck-Einstein.
  pub fn energy(&self) -> f64 {
    HBAR * self.frequency
  }

  pub fn wavelength(&self) -> f64 {
    TAU / self.wavenumber
  }
  pub fn period(&self) -> f64 {
    TAU / self.frequency
  }

  /// The speed a surface of constant phase moves at.
  pub fn phase_velocity(&self) -> f64 {
    self.frequency / self.wavenumber
  }
}

#[cfg(test)]
mod tests {
  use super::PlaneWave;

  #[test]
  fn unit_modulus() {
    let wave = PlaneWave::new(3.0, 7.0);
    for step in 0..100 {
      let position = step as f64 * 0.1;
      let time = step as f64 * 0.07;
      assert!((wave.at(position, time).norm() - 1.0).abs() < 1e-12);
    }
  }

  #[test]
  fn phase_is_constant_along_the_phase_velocity() {
    let wave = PlaneWave::new(3.0, 7.0);
    let velocity = wave.phase_velocity();
    let reference = wave.at(0.0, 0.0);
    for step in 0..100 {
      let time = step as f64 * 0.13;
      let position = velocity * time;
      assert!((wave.at(position, time) - reference).norm() < 1e-12);
    }
  }
}
