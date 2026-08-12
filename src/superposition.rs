use crate::{Complex, dispersion::Dispersion, plane_wave::PlaneWave};

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
}
