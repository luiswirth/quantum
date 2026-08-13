use crate::{
  Complex,
  exact::combination::{Combination, Weighted},
  units::HBAR,
};
use std::f64::consts::TAU;
use std::ops::Mul;

/// `e^(i k x)`, a wave across space at one instant.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicWave {
  pub wavenumber: f64,
}

impl HarmonicWave {
  pub fn new(wavenumber: f64) -> Self {
    Self { wavenumber }
  }

  /// `k x`
  pub fn phase(&self, position: f64) -> f64 {
    self.wavenumber * position
  }
  pub fn at(&self, position: f64) -> Complex {
    Complex::from_polar(1.0, self.phase(position))
  }

  /// `p = planck k`
  pub fn momentum(&self) -> f64 {
    HBAR * self.wavenumber
  }

  /// `lambda = (2 pi) / k`
  pub fn wavelength(&self) -> f64 {
    TAU / self.wavenumber
  }
}

/// `psi(t, x) = e^(i (k x - omega t))`, whose sign convention holds crate wide.
#[derive(Clone, Copy, Debug)]
pub struct PlaneWave {
  pub frequency: f64,
  pub harmonic: HarmonicWave,
}

impl PlaneWave {
  pub fn new(frequency: f64, harmonic: HarmonicWave) -> Self {
    Self {
      frequency,
      harmonic,
    }
  }

  /// `k x - omega t`
  pub fn phase(&self, time: f64, position: f64) -> f64 {
    self.harmonic.phase(position) - self.frequency * time
  }
  pub fn at(&self, time: f64, position: f64) -> Complex {
    Complex::from_polar(1.0, self.phase(time, position))
  }

  /// `e^(-i omega t)`, the factor a wait of `time` turns the wave by.
  pub fn phase_factor(&self, time: f64) -> Complex {
    Complex::from_polar(1.0, -self.frequency * time)
  }

  /// `E = planck omega`
  pub fn energy(&self) -> f64 {
    HBAR * self.frequency
  }

  /// `T = (2 pi) / omega`
  pub fn period(&self) -> f64 {
    TAU / self.frequency
  }

  /// `v_p = omega / k`
  pub fn phase_velocity(&self) -> f64 {
    self.frequency / self.harmonic.wavenumber
  }
}

/// `psi(x) = sum_n c_n e^(i k_n x)`, a state across space at one instant.
pub type Waveform = Combination<HarmonicWave>;

impl Waveform {
  pub fn at(&self, position: f64) -> Complex {
    self
      .terms
      .iter()
      .map(|term| term.amplitude * term.element.at(position))
      .sum()
  }
}

/// `c e^(i (k x - omega t))`
pub type Mode = Weighted<PlaneWave>;

impl Mode {
  pub fn at(&self, time: f64, position: f64) -> Complex {
    self.amplitude * self.element.at(time, position)
  }

  pub fn evolve(&mut self, time: f64) {
    self.amplitude *= self.element.phase_factor(time);
  }
  pub fn evolved(mut self, time: f64) -> Self {
    self.evolve(time);
    self
  }
}

impl Mul<Complex> for HarmonicWave {
  type Output = Weighted<HarmonicWave>;
  fn mul(self, amplitude: Complex) -> Weighted<HarmonicWave> {
    Weighted::new(amplitude, self)
  }
}

impl Mul<Complex> for PlaneWave {
  type Output = Mode;
  fn mul(self, amplitude: Complex) -> Mode {
    Weighted::new(amplitude, self)
  }
}
