use crate::{
  Complex,
  dispersion::Dispersion,
  exact::{
    combination::Combination,
    plane_wave::{Mode, PlaneWave, Waveform},
  },
};

/// `psi(t, x) = sum_n c_n e^(i (k_n x - omega_n t))`
pub type Superposition = Combination<PlaneWave>;

impl Superposition {
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

  pub fn at(&self, time: f64, position: f64) -> Complex {
    self.terms.iter().map(|mode| mode.at(time, position)).sum()
  }

  /// The state at one instant, the frequencies having done their turning.
  pub fn snapshot(&self, time: f64) -> Waveform {
    Waveform::new(
      self
        .terms
        .iter()
        .map(|mode| {
          let turned = mode.evolved(time);
          turned.element.harmonic * turned.amplitude
        })
        .collect(),
    )
  }

  pub fn evolve(&mut self, time: f64) {
    for mode in &mut self.terms {
      mode.evolve(time);
    }
  }
  pub fn evolved(mut self, time: f64) -> Self {
    self.evolve(time);
    self
  }
}

impl Waveform {
  /// The history the medium gives this instant, one frequency per harmonic.
  pub fn lifted(&self, dispersion: Dispersion) -> Superposition {
    Superposition::in_medium(
      dispersion,
      self
        .terms
        .iter()
        .map(|term| (term.amplitude, term.element.wavenumber)),
    )
  }
}
