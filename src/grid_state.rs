use crate::{
  Complex,
  dispersion::Dispersion,
  grid::{Domain, Grid, Momentum, Position},
};
use rustfft::{FftDirection, FftPlanner};
use std::marker::PhantomData;

/// One amplitude per grid point, a vector in `CC^N` that a grid gives
/// coordinates to but does not otherwise constrain.
#[derive(Clone, Debug)]
pub struct GridState<D: Domain> {
  pub values: Vec<Complex>,
  domain: PhantomData<D>,
}

impl<D: Domain> GridState<D> {
  pub fn new(values: Vec<Complex>) -> Self {
    Self {
      values,
      domain: PhantomData,
    }
  }

  pub fn sample(grid: Grid<D>, amplitude: impl Fn(f64) -> Complex) -> Self {
    Self::new(grid.coordinates().map(amplitude).collect())
  }

  pub fn npoints(&self) -> usize {
    self.values.len()
  }

  /// `abs(psi_j)^2`
  pub fn intensities(&self) -> impl Iterator<Item = f64> + '_ {
    self.values.iter().map(|value| value.norm_sqr())
  }
  pub fn total_intensity(&self) -> f64 {
    self.intensities().sum()
  }

  /// `expect(f) = (sum_j abs(psi_j)^2 f(x_j)) / (sum_j abs(psi_j)^2)`
  pub fn expectation(&self, grid: Grid<D>, observable: impl Fn(f64) -> f64) -> f64 {
    assert_eq!(grid.npoints, self.npoints());
    let weighted: f64 = grid
      .coordinates()
      .zip(self.intensities())
      .map(|(coordinate, intensity)| intensity * observable(coordinate))
      .sum();
    weighted / self.total_intensity()
  }

  /// `variance(f) = expect(f^2) - expect(f)^2`
  pub fn variance(&self, grid: Grid<D>, observable: impl Fn(f64) -> f64) -> f64 {
    self.expectation(grid, |coordinate| observable(coordinate).powi(2))
      - self.expectation(grid, &observable).powi(2)
  }

  /// `Delta f = sqrt(variance(f))`
  pub fn uncertainty(&self, grid: Grid<D>, observable: impl Fn(f64) -> f64) -> f64 {
    self.variance(grid, observable).sqrt()
  }

  /// The change of basis into the dual domain, unitary in either direction.
  pub fn transformed(mut self) -> GridState<D::Dual> {
    let direction = if D::SIGN < 0.0 {
      FftDirection::Forward
    } else {
      FftDirection::Inverse
    };
    FftPlanner::new()
      .plan_fft(self.npoints(), direction)
      .process(&mut self.values);
    let scale = (self.npoints() as f64).sqrt().recip();
    for value in &mut self.values {
      *value *= scale;
    }
    GridState::new(self.values)
  }
}

impl GridState<Momentum> {
  /// Evolution is diagonal here, one turn of the phase per wavenumber.
  pub fn evolve(&mut self, time: f64, dispersion: Dispersion, grid: Grid<Momentum>) {
    assert_eq!(grid.npoints, self.npoints());
    for (value, wavenumber) in self.values.iter_mut().zip(grid.coordinates()) {
      *value *= Complex::from_polar(1.0, -dispersion.frequency(wavenumber) * time);
    }
  }
  pub fn evolved(mut self, time: f64, dispersion: Dispersion, grid: Grid<Momentum>) -> Self {
    self.evolve(time, dispersion, grid);
    self
  }
}

impl GridState<Position> {
  /// The exact free propagator: into the domain that diagonalizes evolution,
  /// and back.
  pub fn evolve(&mut self, time: f64, dispersion: Dispersion, grid: Grid<Position>) {
    let here = GridState::<Position>::new(std::mem::take(&mut self.values));
    let there = here.transformed().evolved(time, dispersion, grid.dual());
    self.values = there.transformed().values;
  }
  pub fn evolved(mut self, time: f64, dispersion: Dispersion, grid: Grid<Position>) -> Self {
    self.evolve(time, dispersion, grid);
    self
  }
}
