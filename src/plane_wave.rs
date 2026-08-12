use crate::{Complex, units::HBAR};
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

/// `c e^(i (k x - omega t))`
#[derive(Clone, Copy, Debug)]
pub struct Mode {
  pub amplitude: Complex,
  pub wave: PlaneWave,
}

impl Mode {
  pub fn new(amplitude: Complex, wave: PlaneWave) -> Self {
    Self { amplitude, wave }
  }

  pub fn at(&self, time: f64, position: f64) -> Complex {
    self.amplitude * self.wave.at(time, position)
  }

  /// `abs(c)^2`
  pub fn intensity(&self) -> f64 {
    self.amplitude.norm_sqr()
  }

  pub fn evolve(&mut self, time: f64) {
    self.amplitude *= self.wave.phase_factor(time);
  }
  pub fn evolved(mut self, time: f64) -> Self {
    self.evolve(time);
    self
  }
}

impl From<PlaneWave> for Mode {
  fn from(wave: PlaneWave) -> Self {
    Mode::new(Complex::ONE, wave)
  }
}

impl Mul<Complex> for PlaneWave {
  type Output = Mode;
  fn mul(self, amplitude: Complex) -> Mode {
    Mode::new(amplitude, self)
  }
}

impl Mul<Complex> for Mode {
  type Output = Mode;
  fn mul(self, factor: Complex) -> Mode {
    Mode::new(self.amplitude * factor, self.wave)
  }
}
