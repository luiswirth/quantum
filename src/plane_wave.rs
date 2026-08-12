use crate::{Complex, units::HBAR};
use std::f64::consts::TAU;

/// `psi(x, t) = e^(i (k x - omega t))`, whose sign convention holds crate wide.
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

  /// `k x - omega t`
  pub fn phase(&self, position: f64, time: f64) -> f64 {
    self.wavenumber * position - self.frequency * time
  }
  pub fn at(&self, position: f64, time: f64) -> Complex {
    Complex::from_polar(1.0, self.phase(position, time))
  }

  /// `e^(-i omega t)`, the factor a wait of `time` turns the wave by.
  pub fn phase_factor(&self, time: f64) -> Complex {
    Complex::from_polar(1.0, -self.frequency * time)
  }

  /// `p = planck k`
  pub fn momentum(&self) -> f64 {
    HBAR * self.wavenumber
  }
  /// `E = planck omega`
  pub fn energy(&self) -> f64 {
    HBAR * self.frequency
  }

  /// `lambda = (2 pi) / k`
  pub fn wavelength(&self) -> f64 {
    TAU / self.wavenumber
  }
  /// `T = (2 pi) / omega`
  pub fn period(&self) -> f64 {
    TAU / self.frequency
  }

  /// `v_p = omega / k`
  pub fn phase_velocity(&self) -> f64 {
    self.frequency / self.wavenumber
  }
}
