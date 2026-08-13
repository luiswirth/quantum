use quantum::Complex;
use quantum::discrete::grid::{Grid, Position};
use quantum::discrete::grid_state::GridState;
use quantum::discrete::hamiltonian::Hamiltonian;
use quantum::units::{ELECTRON_MASS, HBAR};
use std::f64::consts::PI;

fn grid() -> Grid<Position> {
  Grid::new(400, 0.05)
}

/// A wave the grid resolves is an eigenstate of the free operator, and its
/// energy is the discrete band, which the lattice bends away from the parabola
/// the continuum has.
#[test]
fn a_resolved_wave_is_a_free_eigenstate() {
  let hamiltonian = Hamiltonian::free(grid(), ELECTRON_MASS);
  let hopping = hamiltonian.hopping(grid());
  for index in [1, 5, 37] {
    let wavenumber = grid().dual().coordinate(index);
    let state = GridState::sample(grid(), |position| {
      Complex::from_polar(1.0, wavenumber * position)
    });
    let acted = hamiltonian.applied(state.clone(), grid());
    let energy = 2.0 * hopping * (1.0 - (wavenumber * grid().spacing).cos());
    for (after, before) in acted.values.iter().zip(&state.values) {
      assert!((after - before * energy).norm() < 1e-9);
    }
  }
}

/// The band is the parabola of the continuum where the wave is long against
/// the spacing, and falls below it where it is not.
#[test]
fn the_band_bends_away_from_the_parabola() {
  let hamiltonian = Hamiltonian::free(grid(), ELECTRON_MASS);
  let hopping = hamiltonian.hopping(grid());
  let band = |wavenumber: f64| 2.0 * hopping * (1.0 - (wavenumber * grid().spacing).cos());
  let parabola = |wavenumber: f64| HBAR.powi(2) * wavenumber.powi(2) / (2.0 * ELECTRON_MASS);

  let long = 0.1 / grid().spacing;
  assert!((band(long) - parabola(long)).abs() < 1e-3 * parabola(long));

  let edge = PI / grid().spacing;
  assert!(band(edge) < 0.5 * parabola(edge));
}

/// A potential the same everywhere is an energy offset and nothing else.
#[test]
fn a_constant_potential_shifts_the_spectrum() {
  let offset = 3.0;
  let free = Hamiltonian::free(grid(), ELECTRON_MASS);
  let shifted = Hamiltonian::sampled(grid(), ELECTRON_MASS, |_| offset);
  for (raised, plain) in shifted.spectrum(grid()).iter().zip(free.spectrum(grid())) {
    assert!((raised - plain - offset).abs() < 1e-9);
  }
}

/// A well with walls far above the levels it holds has the energies of the
/// infinitely deep one, `E_n = (planck^2 pi^2 n^2) / (2 m L^2)`.
///
/// They sit just below those, and both corrections push that way: the wall is
/// finite, so the state reaches past it and the well is effectively the wider
/// for it, and the band of the grid bends below the parabola of the continuum.
#[test]
fn a_deep_well_has_the_analytic_levels() {
  let width = 5.0;
  let hamiltonian = Hamiltonian::sampled(grid(), ELECTRON_MASS, |position| {
    if position.abs() < 0.5 * width {
      0.0
    } else {
      200.0
    }
  });
  let spectrum = hamiltonian.spectrum(grid());
  for level in 1..=3 {
    let analytic =
      HBAR.powi(2) * PI.powi(2) * (level * level) as f64 / (2.0 * ELECTRON_MASS * width.powi(2));
    let computed = spectrum[level - 1];
    assert!(computed < analytic);
    assert!(analytic - computed < 0.01 * analytic);
  }
}

/// The operator is symmetric, so it is its own transpose and its energies are
/// real.
#[test]
fn the_operator_is_symmetric() {
  let hamiltonian = Hamiltonian::sampled(Grid::new(16, 0.2), ELECTRON_MASS, |position| position);
  let matrix = hamiltonian.matrix(Grid::new(16, 0.2));
  assert!((matrix.clone() - matrix.transpose()).norm() < 1e-12);
}
