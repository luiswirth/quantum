use crate::{
  Complex,
  dispersion::Dispersion,
  plane_wave::{Mode, PlaneWave},
};
use std::ops::{Add, Mul};

/// `psi(x, t) = sum_n c_n e^(i (k_n x - omega_n t))`
#[derive(Clone, Debug)]
pub struct Superposition {
  pub modes: Vec<Mode>,
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

  pub fn total_intensity(&self) -> f64 {
    self.modes.iter().map(Mode::intensity).sum()
  }

  /// `angle.l f angle.r = (sum_n abs(c_n)^2 f(k_n)) / (sum_n abs(c_n)^2)`, for
  /// an observable the modes diagonalize.
  pub fn expectation(&self, observable: impl Fn(&PlaneWave) -> f64) -> f64 {
    let weighted: f64 = self
      .modes
      .iter()
      .map(|mode| mode.intensity() * observable(&mode.wave))
      .sum();
    weighted / self.total_intensity()
  }

  /// `angle.l f^2 angle.r - angle.l f angle.r^2`
  pub fn variance(&self, observable: impl Fn(&PlaneWave) -> f64) -> f64 {
    self.expectation(|wave| observable(wave).powi(2)) - self.expectation(&observable).powi(2)
  }

  /// `Delta f = sqrt(angle.l f^2 angle.r - angle.l f angle.r^2)`
  pub fn uncertainty(&self, observable: impl Fn(&PlaneWave) -> f64) -> f64 {
    self.variance(observable).sqrt()
  }

  pub fn evolve(&mut self, time: f64) {
    for mode in &mut self.modes {
      mode.evolve(time);
    }
  }
  pub fn evolved(mut self, time: f64) -> Self {
    self.evolve(time);
    self
  }
}

impl From<Mode> for Superposition {
  fn from(mode: Mode) -> Self {
    Superposition::new(vec![mode])
  }
}

impl Add for Mode {
  type Output = Superposition;
  fn add(self, other: Mode) -> Superposition {
    Superposition::new(vec![self, other])
  }
}

impl Mul<Complex> for Superposition {
  type Output = Superposition;
  fn mul(self, factor: Complex) -> Superposition {
    Superposition::new(self.modes.into_iter().map(|mode| mode * factor).collect())
  }
}

impl Add for Superposition {
  type Output = Superposition;
  fn add(mut self, other: Superposition) -> Superposition {
    self.modes.extend(other.modes);
    self
  }
}
