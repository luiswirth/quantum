use quantum::Complex;
use quantum::dispersion::Dispersion;
use quantum::grid::{Grid, Position};
use quantum::grid_state::GridState;
use quantum::superposition::Superposition;
use std::f64::consts::TAU;

fn grid() -> Grid<Position> {
  Grid::new(256, 0.1)
}

fn gaussian(width: f64) -> GridState<Position> {
  GridState::sample(grid(), |position| {
    Complex::from((-position.powi(2) / (4.0 * width.powi(2))).exp())
  })
}

/// The transform is unitary, so no probability is created or lost by it.
#[test]
fn the_transform_preserves_the_total_intensity() {
  let state = gaussian(1.0);
  let total = state.total_intensity();
  assert!((state.transformed().total_intensity() - total).abs() < 1e-9);
}

/// Transforming twice is the identity, matching the grid's own involution.
#[test]
fn transforming_twice_returns_the_state() {
  let state = gaussian(1.0);
  let returned = state.clone().transformed().transformed();
  for (before, after) in state.values.iter().zip(&returned.values) {
    assert!((before - after).norm() < 1e-12);
  }
}

/// A point source spreads over every wavenumber equally, which is the delta
/// and the plane wave standing on opposite sides of the transform.
#[test]
fn a_point_source_has_a_flat_spectrum() {
  let mut values = vec![Complex::ZERO; grid().npoints];
  values[0] = Complex::ONE;
  let spectrum = GridState::<Position>::new(values).transformed();
  let expected = (grid().npoints as f64).sqrt().recip();
  for value in &spectrum.values {
    assert!((value.norm() - expected).abs() < 1e-12);
  }
}

/// A Gaussian saturates the uncertainty relation, `Delta x Delta k = 1/2`.
#[test]
fn the_gaussian_sits_at_the_uncertainty_bound() {
  for width in [0.5, 1.0, 2.0] {
    let state = gaussian(width);
    let position_spread = state.uncertainty(grid(), |position| position);
    let wavenumber_spread = state
      .transformed()
      .uncertainty(grid().dual(), |wavenumber| wavenumber);
    assert!((position_spread - width).abs() < 1e-6);
    assert!((position_spread * wavenumber_spread - 0.5).abs() < 1e-6);
  }
}

/// Evolving on the grid and evolving the exact modes agree, since both are
/// the same multiplication by a phase.
#[test]
fn grid_evolution_matches_the_exact_modes() {
  let grid: Grid<Position> = Grid::new(64, 0.25);
  let dispersion = Dispersion::electron();
  let step = TAU / grid.length();
  let superposition = Superposition::in_medium(
    dispersion,
    [
      (Complex::ONE, 3.0 * step),
      (Complex::I, -5.0 * step),
      (Complex::new(0.5, 0.5), 12.0 * step),
    ],
  );
  let mut state = superposition.on_grid(0.0, grid);
  for time in [0.0, 0.4, 3.0] {
    let evolved = state.clone().evolved(time, dispersion, grid);
    let exact = superposition.on_grid(time, grid);
    for (here, there) in evolved.values.iter().zip(&exact.values) {
      assert!((here - there).norm() < 1e-10);
    }
  }
  state.evolve(1.0, dispersion, grid);
  assert!(
    (state.total_intensity() - superposition.on_grid(1.0, grid).total_intensity()).abs() < 1e-9
  );
}

/// A packet spreads in position and holds its shape in momentum, since
/// evolution only turns the phases there.
#[test]
fn evolution_spreads_the_packet_and_leaves_the_spectrum() {
  let dispersion = Dispersion::electron();
  let state = gaussian(1.0);
  let spread = state.uncertainty(grid(), |position| position);
  let wavenumber_spread = state
    .clone()
    .transformed()
    .uncertainty(grid().dual(), |k| k);
  let later = state.evolved(2.0, dispersion, grid());
  assert!(later.uncertainty(grid(), |position| position) > spread);
  assert!((later.transformed().uncertainty(grid().dual(), |k| k) - wavenumber_spread).abs() < 1e-9);
}
