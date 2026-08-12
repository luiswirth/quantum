use quantum::Complex;
use quantum::dirac::{Comb, Dirac, Event, Excitation, Impulse, Source};
use quantum::grid::{Grid, Position};
use quantum::plane_wave::{HarmonicWave, PlaneWave};
use std::f64::consts::TAU;

fn grid() -> Grid<Position> {
  Grid::new(32, 0.5)
}

/// A Dirac reads a function at one point, which is all it does.
#[test]
fn an_impulse_samples_the_test_function() {
  let impulse = Impulse::new(Complex::new(2.0, -1.0), Dirac::new(1.25));
  let harmonic = HarmonicWave::new(3.0);
  let taken = impulse.apply(|position| harmonic.at(position));
  assert!((taken - Complex::new(2.0, -1.0) * harmonic.at(1.25)).norm() < 1e-12);
}

/// The transform of a point source is flat, the mirror of a plane wave
/// transforming to a point.
#[test]
fn a_point_source_has_a_flat_transform() {
  let impulse = Impulse::new(Complex::ONE, Dirac::new(-0.75));
  for step in -20..=20 {
    let wavenumber = step as f64 * 0.3;
    assert!((impulse.transform_at(wavenumber).norm() - 1.0).abs() < 1e-12);
  }
}

/// The sampling comb transforms to the comb of the dual grid: it is flat on
/// the dual lattice and cancels everywhere else.
#[test]
fn the_comb_transforms_to_the_dual_comb() {
  let grid = grid();
  let comb = Comb::sampling(grid);
  let teeth = grid.npoints as f64;
  for step in -3..=3 {
    let on_lattice = step as f64 * TAU / grid.spacing;
    assert!((comb.transform_at(on_lattice).norm() - teeth).abs() < 1e-9);
    let between = on_lattice + TAU / grid.length();
    assert!(comb.transform_at(between).norm() < 1e-9);
  }
}

/// Weights land on the sample they fall nearest, and off-lattice ones snap.
#[test]
fn the_comb_gathers_onto_the_grid() {
  let grid = grid();
  let comb = Comb::new(vec![
    Impulse::new(Complex::ONE, Dirac::new(0.0)),
    Impulse::new(Complex::I, Dirac::new(1.02)),
  ]);
  let state = comb.on_grid(grid);
  assert_eq!(state.values[0], Complex::ONE);
  assert_eq!(state.values[2], Complex::I);
  assert_eq!(state.total_intensity(), 2.0);
}

/// An event and a plane wave are each other's transform, so pairing one with
/// the other gives the same number read either way.
#[test]
fn the_event_and_the_plane_wave_pair_to_the_same_number() {
  let event = Event::new(0.6, Dirac::new(-1.4));
  let wave = PlaneWave::new(2.2, HarmonicWave::new(1.7));
  let taken = event.apply(|time, position| wave.at(time, position));
  let transformed = event.transform_at(wave.frequency, wave.harmonic.wavenumber);
  assert!((taken - transformed.conj()).norm() < 1e-12);
}

/// An excitation firing twice at the same place is a source that beats in
/// time: its transform vanishes where the two firings cancel.
#[test]
fn two_firings_interfere_in_frequency() {
  let excitation = Excitation::new(vec![
    Source::new(Complex::ONE, Event::new(0.0, Dirac::new(0.0))),
    Source::new(Complex::ONE, Event::new(1.0, Dirac::new(0.0))),
  ]);
  assert!((excitation.transform_at(0.0, 0.0).norm() - 2.0).abs() < 1e-12);
  assert!(excitation.transform_at(std::f64::consts::PI, 0.0).norm() < 1e-12);
}
