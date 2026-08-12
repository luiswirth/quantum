use quantum::plane_wave::PlaneWave;

#[test]
fn unit_modulus() {
  let wave = PlaneWave::new(3.0, 7.0);
  for step in 0..100 {
    let position = step as f64 * 0.1;
    let time = step as f64 * 0.07;
    assert!((wave.at(position, time).norm() - 1.0).abs() < 1e-12);
  }
}

#[test]
fn phase_is_constant_along_the_phase_velocity() {
  let wave = PlaneWave::new(3.0, 7.0);
  let velocity = wave.phase_velocity();
  let reference = wave.at(0.0, 0.0);
  for step in 0..100 {
    let time = step as f64 * 0.13;
    let position = velocity * time;
    assert!((wave.at(position, time) - reference).norm() < 1e-12);
  }
}
