use quantum::Complex;
use quantum::grid::{Grid, Position};
use quantum::grid_state::GridState;

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
