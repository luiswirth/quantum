use quantum::dispersion::Dispersion;
use quantum::units::{ELECTRON_MASS, HBAR};

fn wavenumbers() -> impl Iterator<Item = f64> {
  (-20..=20)
    .filter(|index| *index != 0)
    .map(|index| index as f64 * 0.5)
}

#[test]
fn group_velocity_is_the_slope_of_the_dispersion_relation() {
  let step = 1e-6;
  for dispersion in [Dispersion::electron(), Dispersion::vacuum_light()] {
    for wavenumber in wavenumbers() {
      let slope = (dispersion.frequency(wavenumber + step)
        - dispersion.frequency(wavenumber - step))
        / (2.0 * step);
      assert!((slope - dispersion.group_velocity(wavenumber)).abs() < 1e-6);
    }
  }
}

/// The two velocities part ways exactly where the relation bends.
#[test]
fn light_carries_phase_and_packet_together_and_the_electron_does_not() {
  let light = Dispersion::vacuum_light();
  let electron = Dispersion::electron();
  for wavenumber in wavenumbers() {
    assert!((light.phase_velocity(wavenumber) - light.group_velocity(wavenumber)).abs() < 1e-12);
    assert!(
      (2.0 * electron.phase_velocity(wavenumber) - electron.group_velocity(wavenumber)).abs()
        < 1e-12
    );
  }
}

/// The dispersion relation says that the wave's energy is the particle's
/// kinetic energy.
#[test]
fn energy_is_kinetic() {
  let mass = ELECTRON_MASS;
  let dispersion = Dispersion::FreeParticle { mass };
  for wavenumber in wavenumbers() {
    let wave = dispersion.plane_wave(wavenumber);
    let kinetic = wave.momentum().powi(2) / (2.0 * mass);
    assert!((wave.energy() - kinetic).abs() < 1e-12);
  }
}

/// The combination every effective mass model is quoted in.
#[test]
fn electron_energy_scale() {
  let scale = HBAR.powi(2) / (2.0 * ELECTRON_MASS);
  assert!((scale - 0.0380998).abs() < 1e-7);
}

#[test]
fn curvature_is_the_bend_of_the_dispersion_relation() {
  let step = 1e-3;
  for dispersion in [Dispersion::electron(), Dispersion::vacuum_light()] {
    for wavenumber in wavenumbers() {
      let bend = (dispersion.frequency(wavenumber + step) - 2.0 * dispersion.frequency(wavenumber)
        + dispersion.frequency(wavenumber - step))
        / step.powi(2);
      assert!((bend - dispersion.curvature(wavenumber)).abs() < 1e-6);
    }
  }
}

#[test]
fn the_dispersion_relation_is_the_energy_momentum_relation() {
  let dispersion = Dispersion::electron();
  for wavenumber in wavenumbers() {
    let wave = dispersion.plane_wave(wavenumber);
    assert!((dispersion.energy(wave.momentum()) - wave.energy()).abs() < 1e-12);
  }
}

#[test]
fn the_free_particle_responds_with_its_own_mass() {
  let dispersion = Dispersion::FreeParticle {
    mass: ELECTRON_MASS,
  };
  for wavenumber in wavenumbers() {
    assert!((dispersion.effective_mass(wavenumber) - ELECTRON_MASS).abs() < 1e-12);
    assert!(
      Dispersion::vacuum_light()
        .effective_mass(wavenumber)
        .is_infinite()
    );
  }
}
