use quantum::dispersion::Dispersion;
use quantum::plane_wave::Mode;
use quantum::superposition::Superposition;
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

/// Evolution is a rotation of the amplitudes, and nothing else.
#[test]
fn evolving_the_amplitudes_is_waiting() {
  let superposition = Superposition::in_medium(
    Dispersion::electron(),
    [(Complex::ONE, 1.5), (Complex::I, -0.5), (Complex::ONE, 3.0)],
  );
  for time in [0.0, 0.2, 1.7] {
    let waited = superposition.evolved(time);
    for position in positions() {
      assert!((waited.at(position, 0.0) - superposition.at(position, time)).norm() < 1e-12);
    }
  }
}

/// Waiting twice is waiting once for the sum of the two waits.
#[test]
fn evolution_composes() {
  let superposition =
    Superposition::in_medium(Dispersion::electron(), [(Complex::ONE, 1.5), (Complex::I, -0.5)]);
  let stepwise = superposition.evolved(0.3).evolved(0.9);
  let direct = superposition.evolved(1.2);
  for position in positions() {
    assert!((stepwise.at(position, 0.0) - direct.at(position, 0.0)).norm() < 1e-12);
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
    let expected = factor * left.at(position, 0.4) + right.at(position, 0.4);
    assert!((combined.at(position, 0.4) - expected).norm() < 1e-12);
  }
}

/// A wave scaled is a mode, and two modes added are a superposition.
#[test]
fn the_levels_build_on_each_other() {
  let dispersion = Dispersion::electron();
  let built = dispersion.plane_wave(1.0) * Complex::ONE + dispersion.plane_wave(-1.0) * Complex::ONE;
  let direct = Superposition::in_medium(dispersion, [(Complex::ONE, 1.0), (Complex::ONE, -1.0)]);
  for position in positions() {
    assert!((built.at(position, 0.6) - direct.at(position, 0.6)).norm() < 1e-12);
  }
}
