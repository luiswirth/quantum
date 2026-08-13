//! The bridge between the exact states and the grid, in both directions.
//!
//! Each crossing loses something, and the loss is the subject of this module.
//! Sampling forgets everything between the samples, and gathering forgets
//! where inside its cell a weight sat.

use crate::{
  Complex,
  discrete::{
    grid::{Grid, Momentum, Position},
    grid_state::GridState,
  },
  dispersion::Dispersion,
  exact::{
    dirac::{Comb, Dirac, Impulse},
    plane_wave::{HarmonicWave, Waveform},
    superposition::Superposition,
  },
};

impl Comb {
  /// A unit tooth on every sample, the operator that turns a function into its
  /// samples.
  pub fn sampling(grid: Grid<Position>) -> Self {
    Self::new(
      grid
        .coordinates()
        .map(|position| Impulse::new(Complex::ONE, Dirac::new(position)))
        .collect(),
    )
  }

  /// One tooth per sample, carrying its amplitude.
  pub fn from_grid(state: &GridState<Position>, grid: Grid<Position>) -> Self {
    assert_eq!(grid.npoints, state.npoints());
    Self::new(
      grid
        .coordinates()
        .zip(&state.values)
        .map(|(position, amplitude)| Dirac::new(position) * *amplitude)
        .collect(),
    )
  }

  /// The weights gathered onto the samples they fall nearest, which is what
  /// aliasing is on this side of the transform.
  pub fn on_grid(&self, grid: Grid<Position>) -> GridState<Position> {
    let mut values = vec![Complex::ZERO; grid.npoints];
    for term in &self.terms {
      values[grid.nearest_index(term.element.position)] += term.amplitude;
    }
    GridState::new(values)
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
}

impl Superposition {
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
}
