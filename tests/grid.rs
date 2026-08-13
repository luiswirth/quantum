use quantum::discrete::grid::{Grid, Position};
use std::f64::consts::TAU;

fn grid() -> Grid<Position> {
  Grid::new(8, 0.25)
}

#[test]
fn dualizing_twice_returns_the_grid() {
  assert_eq!(grid().dual().dual(), grid());
}

/// Discreteness on one side is periodicity on the other.
#[test]
fn the_dual_trades_spacing_for_length() {
  let grid = grid();
  let dual = grid.dual();
  assert!((dual.spacing - TAU / grid.length()).abs() < 1e-12);
  assert!((dual.length() - TAU / grid.spacing).abs() < 1e-12);
}

#[test]
fn the_samples_straddle_the_origin() {
  let grid = grid();
  let coordinates: Vec<f64> = grid.coordinates().collect();
  assert_eq!(coordinates[0], 0.0);
  assert!(coordinates.iter().cloned().fold(f64::NAN, f64::min) >= -grid.length() / 2.0);
  assert!(coordinates.iter().cloned().fold(f64::NAN, f64::max) < grid.length() / 2.0);
}

/// The largest wavenumber the grid resolves, where the samples alternate sign.
#[test]
fn the_dual_reaches_the_nyquist_wavenumber() {
  let grid = grid();
  let reach = grid
    .dual()
    .coordinates()
    .fold(f64::NAN, |extreme: f64, wavenumber| {
      extreme.max(wavenumber.abs())
    });
  assert!((reach - std::f64::consts::PI / grid.spacing).abs() < 1e-12);
}
