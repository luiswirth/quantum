use quantum::dispersion::Dispersion;
use quantum::superposition::{Mode, Superposition};
use quantum::Complex;

fn positions() -> impl Iterator<Item = f64> {
  (-30..=30).map(|index| index as f64 * 0.1)
}

#[test]
fn an_empty_superposition_vanishes() {
  let superposition = Superposition::new(Vec::new());
  for position in positions() {
    assert_eq!(superposition.at(position, 0.7), Complex::ZERO);
  }
}

#[test]
fn a_single_mode_is_its_own_wave() {
  let wave = Dispersion::electron().plane_wave(2.0);
  let superposition = Superposition::new(vec![Mode::new(Complex::ONE, wave)]);
  for position in positions() {
    assert!((superposition.at(position, 0.3) - wave.at(position, 0.3)).norm() < 1e-12);
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
      assert!((superposition.at(position, time).norm() - envelope.abs()).abs() < 1e-12);
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
    assert!((superposition.at(crest, time).norm() - 2.0).abs() < 1e-9);
  }
}
