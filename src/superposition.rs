use crate::{Complex, dispersion::Dispersion, plane_wave::PlaneWave};
use std::ops::{Add, Mul};

/// `c e^(i (k x - omega t))`
#[derive(Clone, Copy, Debug)]
pub struct Mode {
  pub amplitude: Complex,
  pub wave: PlaneWave,
}

/// `psi(x, t) = sum_n c_n e^(i (k_n x - omega_n t))`
#[derive(Clone, Debug)]
pub struct Superposition {
  pub modes: Vec<Mode>,
}

impl Mode {
  pub fn new(amplitude: Complex, wave: PlaneWave) -> Self {
    Self { amplitude, wave }
  }

  pub fn at(&self, position: f64, time: f64) -> Complex {
    self.amplitude * self.wave.at(position, time)
  }

  pub fn evolve(&mut self, time: f64) {
    self.amplitude *= self.wave.phase_factor(time);
  }
  pub fn evolved(&self, time: f64) -> Self {
    let mut mode = *self;
    mode.evolve(time);
    mode
  }
}

impl Superposition {
  pub fn new(modes: Vec<Mode>) -> Self {
    Self { modes }
  }

  /// The modes a medium admits, whose frequencies it fixes through `omega(k)`.
  pub fn in_medium(
    dispersion: Dispersion,
    modes: impl IntoIterator<Item = (Complex, f64)>,
  ) -> Self {
    Self::new(
      modes
        .into_iter()
        .map(|(amplitude, wavenumber)| Mode::new(amplitude, dispersion.plane_wave(wavenumber)))
        .collect(),
    )
  }

  pub fn at(&self, position: f64, time: f64) -> Complex {
    self.modes.iter().map(|mode| mode.at(position, time)).sum()
  }

  pub fn evolve(&mut self, time: f64) {
    for mode in &mut self.modes {
      mode.evolve(time);
    }
  }
  pub fn evolved(&self, time: f64) -> Self {
    let mut superposition = self.clone();
    superposition.evolve(time);
    superposition
  }
}

impl From<PlaneWave> for Mode {
  fn from(wave: PlaneWave) -> Self {
    Mode::new(Complex::ONE, wave)
  }
}

impl From<Mode> for Superposition {
  fn from(mode: Mode) -> Self {
    Superposition::new(vec![mode])
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

impl Mul<Complex> for Superposition {
  type Output = Superposition;
  fn mul(self, factor: Complex) -> Superposition {
    Superposition::new(self.modes.into_iter().map(|mode| mode * factor).collect())
  }
}

impl Add for Mode {
  type Output = Superposition;
  fn add(self, other: Mode) -> Superposition {
    Superposition::new(vec![self, other])
  }
}

impl Add for Superposition {
  type Output = Superposition;
  fn add(mut self, other: Superposition) -> Superposition {
    self.modes.extend(other.modes);
    self
  }
}
