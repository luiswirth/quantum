use quantum::Complex;
use quantum::discrete::grid::{Grid, Position};
use quantum::dispersion::Dispersion;
use quantum::exact::operator::DiagonalOperator;
use quantum::exact::plane_wave::Mode;
use quantum::exact::superposition::Superposition;
use std::f64::consts::TAU;

fn positions() -> impl Iterator<Item = f64> {
  (-30..=30).map(|index| index as f64 * 0.1)
}

#[test]
fn an_empty_superposition_vanishes() {
  let superposition = Superposition::new(Vec::new());
  for position in positions() {
    assert_eq!(superposition.at(0.7, position), Complex::ZERO);
  }
}

#[test]
fn a_single_mode_is_its_own_wave() {
  let wave = Dispersion::electron().plane_wave(2.0);
  let superposition = Superposition::new(vec![Mode::new(Complex::ONE, wave)]);
  for position in positions() {
    assert!((superposition.at(0.3, position) - wave.at(0.3, position)).norm() < 1e-12);
  }
}

/// Two modes of opposite wavenumber share a frequency, so their sum only
/// breathes in place: the nodes never move.
#[test]
fn counter_propagating_modes_stand_still() {
  let dispersion = Dispersion::electron();
  let wavenumber = 2.0;
  let superposition = Superposition::in_medium(
    dispersion,
    [(Complex::ONE, wavenumber), (Complex::ONE, -wavenumber)],
  );
  for position in positions() {
    let envelope = 2.0 * (wavenumber * position).cos();
    for time in [0.0, 0.4, 1.3] {
      assert!((superposition.at(time, position).norm() - envelope.abs()).abs() < 1e-12);
    }
  }
}

/// Two modes of neighbouring wavenumber beat, and the beat travels at the
/// group velocity while the wave inside it travels at the phase velocity.
#[test]
fn a_beat_travels_at_the_group_velocity() {
  let dispersion = Dispersion::electron();
  let center = 2.0;
  let offset = 0.05;
  let superposition = Superposition::in_medium(
    dispersion,
    [
      (Complex::ONE, center - offset),
      (Complex::ONE, center + offset),
    ],
  );
  let velocity = dispersion.group_velocity(center);
  for time in [0.0, 0.5, 2.0, 5.0] {
    let crest = velocity * time;
    assert!((superposition.at(time, crest).norm() - 2.0).abs() < 1e-9);
  }
}

/// Evolution is a rotation of the amplitudes, and nothing else.
#[test]
fn evolving_the_amplitudes_is_waiting() {
  let superposition = Superposition::in_medium(
    Dispersion::electron(),
    [(Complex::ONE, 1.5), (Complex::I, -0.5), (Complex::ONE, 3.0)],
  );
  for time in [0.0, 0.2, 1.7] {
    let waited = superposition.clone().evolved(time);
    for position in positions() {
      assert!((waited.at(0.0, position) - superposition.at(time, position)).norm() < 1e-12);
    }
  }
}

/// Waiting twice is waiting once for the sum of the two waits.
#[test]
fn evolution_composes() {
  let superposition = Superposition::in_medium(
    Dispersion::electron(),
    [(Complex::ONE, 1.5), (Complex::I, -0.5)],
  );
  let stepwise = superposition.clone().evolved(0.3).evolved(0.9);
  let direct = superposition.evolved(1.2);
  for position in positions() {
    assert!((stepwise.at(0.0, position) - direct.at(0.0, position)).norm() < 1e-12);
  }
}

#[test]
fn the_superposition_is_linear() {
  let dispersion = Dispersion::electron();
  let left = Superposition::in_medium(dispersion, [(Complex::ONE, 1.0)]);
  let right = Superposition::in_medium(dispersion, [(Complex::I, 2.0)]);
  let factor = Complex::new(0.3, -1.4);
  let combined = left.clone() * factor + right.clone();
  for position in positions() {
    let expected = factor * left.at(0.4, position) + right.at(0.4, position);
    assert!((combined.at(0.4, position) - expected).norm() < 1e-12);
  }
}

/// A wave scaled is a mode, and two modes added are a superposition.
#[test]
fn the_levels_build_on_each_other() {
  let dispersion = Dispersion::electron();
  let built =
    dispersion.plane_wave(1.0) * Complex::ONE + dispersion.plane_wave(-1.0) * Complex::ONE;
  let direct = Superposition::in_medium(dispersion, [(Complex::ONE, 1.0), (Complex::ONE, -1.0)]);
  for position in positions() {
    assert!((built.at(0.6, position) - direct.at(0.6, position)).norm() < 1e-12);
  }
}

/// The moduli are all the expectations see, so evolution leaves them alone.
#[test]
fn momentum_and_energy_are_conserved() {
  let superposition = Superposition::in_medium(
    Dispersion::electron(),
    [
      (Complex::ONE, 1.0),
      (Complex::new(0.0, 2.0), 2.5),
      (Complex::ONE, -0.5),
    ],
  );
  let momentum = || DiagonalOperator::momentum().on_plane_waves();
  let expected_momentum = momentum().expectation(&superposition);
  let expected_energy = DiagonalOperator::energy().expectation(&superposition);
  let spread = momentum().uncertainty(&superposition);
  for time in [0.3, 2.0, 11.0] {
    let waited = superposition.clone().evolved(time);
    assert!((momentum().expectation(&waited) - expected_momentum).norm() < 1e-12);
    assert!((DiagonalOperator::energy().expectation(&waited) - expected_energy).norm() < 1e-12);
    assert!((momentum().uncertainty(&waited) - spread).abs() < 1e-12);
  }
}

/// A single mode is a state of definite momentum.
#[test]
fn one_mode_has_no_uncertainty() {
  let superposition =
    Superposition::in_medium(Dispersion::electron(), [(Complex::new(3.0, 4.0), 2.0)]);
  assert!((superposition.total_intensity() - 25.0).abs() < 1e-12);
  let wavenumber = || DiagonalOperator::wavenumber().on_plane_waves();
  assert!((wavenumber().expectation(&superposition) - 2.0).norm() < 1e-12);
  assert!(wavenumber().uncertainty(&superposition).abs() < 1e-12);
}

/// Two modes weighted equally sit halfway between them, and the uncertainty is
/// their half separation.
#[test]
fn two_modes_straddle_their_center() {
  let superposition = Superposition::in_medium(
    Dispersion::electron(),
    [(Complex::ONE, 1.0), (Complex::ONE, 3.0)],
  );
  let wavenumber = || DiagonalOperator::wavenumber().on_plane_waves();
  assert!((wavenumber().expectation(&superposition) - 2.0).norm() < 1e-12);
  assert!((wavenumber().uncertainty(&superposition) - 1.0).abs() < 1e-12);
}

/// Sampling a superposition and transforming it recovers the amplitudes it
/// was built from, as long as the grid resolves every mode.
#[test]
fn the_grid_recovers_the_modes_it_resolves() {
  let grid: Grid<Position> = Grid::new(64, 0.25);
  let dispersion = Dispersion::electron();
  let wavenumbers: Vec<f64> = [3, -5, 12]
    .iter()
    .map(|n| *n as f64 * TAU / grid.length())
    .collect();
  let amplitudes = [
    Complex::ONE,
    Complex::new(0.0, -2.0),
    Complex::new(0.5, 0.5),
  ];
  let superposition = Superposition::in_medium(
    dispersion,
    amplitudes.iter().copied().zip(wavenumbers.iter().copied()),
  );

  let spectrum = superposition.on_grid(0.0, grid).transformed();
  let recovered = Superposition::from_grid(&spectrum, grid.dual(), dispersion);

  for position in positions() {
    assert!((recovered.at(0.9, position) - superposition.at(0.9, position)).norm() < 1e-10);
  }
}

/// A mode beyond the Nyquist wavenumber is indistinguishable on the grid from
/// one inside it.
#[test]
fn the_grid_aliases_what_it_cannot_resolve() {
  let grid: Grid<Position> = Grid::new(64, 0.25);
  let dispersion = Dispersion::electron();
  let step = TAU / grid.length();
  let inside = Superposition::in_medium(dispersion, [(Complex::ONE, 3.0 * step)]);
  let outside = Superposition::in_medium(dispersion, [(Complex::ONE, 67.0 * step)]);
  for (here, there) in inside
    .on_grid(0.0, grid)
    .values
    .iter()
    .zip(&outside.on_grid(0.0, grid).values)
  {
    assert!((here - there).norm() < 1e-10);
  }
}

/// A history seen at one instant, and an instant given a history by the
/// medium, are inverse to each other.
#[test]
fn the_instant_and_the_history_determine_each_other() {
  let dispersion = Dispersion::electron();
  let superposition = Superposition::in_medium(
    dispersion,
    [(Complex::ONE, 1.5), (Complex::I, -0.5), (Complex::ONE, 3.0)],
  );
  for time in [0.0, 0.7, 4.0] {
    let waveform = superposition.snapshot(time);
    let lifted = waveform.lifted(dispersion);
    for position in positions() {
      assert!((waveform.at(position) - superposition.at(time, position)).norm() < 1e-12);
      assert!((lifted.at(0.0, position) - superposition.at(time, position)).norm() < 1e-12);
    }
  }
}
