use crate::{
  Complex,
  combination::Combination,
  dispersion::Dispersion,
  grid::{Grid, Momentum, Position},
  grid_state::GridState,
  plane_wave::{HarmonicWave, Mode, PlaneWave, Waveform},
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

  /// `psi(t, x_j)`, which aliases whatever the grid does not resolve.
  pub fn on_grid(&self, time: f64, grid: Grid<Position>) -> GridState<Position> {
    self.snapshot(time).on_grid(grid)
  }

  /// The instant the grid holds, given a history by the medium.
  pub fn from_grid(
    spectrum: &GridState<Momentum>,
    grid: Grid<Momentum>,
    dispersion: Dispersion,
  ) -> Self {
    Waveform::from_grid(spectrum, grid).lifted(dispersion)
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
  /// `psi(x_j)`, which aliases whatever the grid does not resolve.
  pub fn on_grid(&self, grid: Grid<Position>) -> GridState<Position> {
    GridState::new(
      grid
        .coordinates()
        .map(|position| self.at(position))
        .collect(),
    )
  }

  /// One harmonic per grid wavenumber, needing no medium at all.
  ///
  /// The transform is unitary while the sum is not, so the amplitudes lose the
  /// `sqrt(N)` that sampling gave them.
  pub fn from_grid(spectrum: &GridState<Momentum>, grid: Grid<Momentum>) -> Self {
    assert_eq!(grid.npoints, spectrum.npoints());
    let scale = (grid.npoints as f64).sqrt().recip();
    Self::new(
      grid
        .coordinates()
        .zip(&spectrum.values)
        .map(|(wavenumber, amplitude)| HarmonicWave::new(wavenumber) * (amplitude * scale))
        .collect(),
    )
  }

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
