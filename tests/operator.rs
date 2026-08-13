use quantum::Complex;
use quantum::dispersion::Dispersion;
use quantum::exact::combination::Combination;
use quantum::exact::dirac::Dirac;
use quantum::exact::operator::DiagonalOperator;
use quantum::exact::plane_wave::HarmonicWave;
use quantum::units::{ELECTRON_MASS, HBAR};
use std::f64::consts::TAU;

fn harmonic(wavenumber: f64) -> Combination<HarmonicWave> {
  (HarmonicWave::new(wavenumber) * Complex::ONE).into()
}

/// An eigenvector of an operator has that operator's eigenvalue and no spread.
#[test]
fn one_harmonic_has_a_sharp_momentum() {
  let state = harmonic(3.0);
  let momentum = DiagonalOperator::momentum();
  assert!((momentum.expectation(&state) - HBAR * 3.0).norm() < 1e-30);
  assert!(momentum.uncertainty(&state) < 1e-30);
}

/// The basis diagonalizes the operator, so acting only rescales the amplitudes
/// and never leaves the span.
#[test]
fn acting_rescales_the_amplitudes() {
  let state = Combination::new(vec![
    HarmonicWave::new(1.0) * Complex::ONE,
    HarmonicWave::new(-2.0) * Complex::I,
  ]);
  let acted = DiagonalOperator::wavenumber().applied(&state);
  assert_eq!(acted.terms.len(), state.terms.len());
  for (after, before) in acted.terms.iter().zip(&state.terms) {
    assert_eq!(after.element.wavenumber, before.element.wavenumber);
    assert!((after.amplitude - before.amplitude * before.element.wavenumber).norm() < 1e-30);
  }
}

/// A free particle carries no potential of its own, so the whole energy of a
/// mode is the kinetic one.
#[test]
fn a_free_medium_is_kinetic_energy_alone() {
  let dispersion = Dispersion::electron();
  let kinetic = DiagonalOperator::kinetic(ELECTRON_MASS);
  for wavenumber in [-2.0, 0.0, 0.5, 7.0] {
    let harmonic = HarmonicWave::new(wavenumber);
    let expected = dispersion.hamiltonian().eigenvalue(&harmonic);
    let difference = expected - kinetic.eigenvalue(&harmonic);
    assert!(difference.norm() <= 1e-12 * expected.norm());
  }
}

/// The velocity is the one the packets of the same medium travel at.
#[test]
fn the_velocity_operator_is_the_group_velocity() {
  let dispersion = Dispersion::electron();
  let velocity = dispersion.velocity();
  for wavenumber in [-3.0, 1.0, 4.0] {
    let expected = dispersion.group_velocity(wavenumber);
    assert!((velocity.eigenvalue(&HarmonicWave::new(wavenumber)) - expected).norm() < 1e-30);
  }
}

/// Translation is unitary, and its eigenvalue is the phase Bloch's theorem
/// gives a wave carried one period along a lattice. A wave whose wavelength is
/// that period comes back to itself.
#[test]
fn translation_carries_the_bloch_phase() {
  let period = 2.5;
  let translation = DiagonalOperator::translation(period);
  for wavenumber in [-1.0, 0.3, 6.0] {
    assert!(
      (translation
        .eigenvalue(&HarmonicWave::new(wavenumber))
        .norm()
        - 1.0)
        .abs()
        < 1e-12
    );
  }
  let commensurate = translation.eigenvalue(&HarmonicWave::new(TAU / period));
  assert!((commensurate - Complex::ONE).norm() < 1e-12);
}

/// The potential is a function of the position operator, so it reads the
/// profile at the point each tooth sits on.
#[test]
fn the_potential_reads_the_profile_at_the_tooth() {
  let well = DiagonalOperator::potential(|position| if position.abs() < 1.0 { 0.0 } else { 5.0 });
  assert!((well.eigenvalue(&Dirac::new(0.5)) - 0.0).norm() < 1e-30);
  assert!((well.eigenvalue(&Dirac::new(3.0)) - 5.0).norm() < 1e-30);
  assert!((DiagonalOperator::position().eigenvalue(&Dirac::new(3.0)) - 3.0).norm() < 1e-30);
}

/// The identity is diagonal in every basis and leaves the state alone.
#[test]
fn the_identity_leaves_the_state_alone() {
  let state = harmonic(2.0);
  let acted = DiagonalOperator::identity().applied(&state);
  assert!((acted.terms[0].amplitude - state.terms[0].amplitude).norm() < 1e-30);
}
