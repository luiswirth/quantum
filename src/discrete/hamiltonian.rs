use crate::{
  Complex,
  discrete::{
    grid::{Grid, Position},
    grid_state::GridState,
  },
  units::HBAR,
};
use nalgebra::{DMatrix, SymmetricEigen};

/// `hat(H) = - planck^2 / (2 m) (dif^2) / (dif x^2) + V(x)`, sampled on a grid.
///
/// The second derivative is the three point stencil, so the operator couples
/// each point to its two neighbors and to nothing else. The grid is periodic,
/// so the first point and the last are neighbors too, and a wall is made of
/// potential rather than of a boundary condition.
#[derive(Clone, Debug)]
pub struct Hamiltonian {
  pub mass: f64,
  pub potential: Vec<f64>,
}

impl Hamiltonian {
  pub fn new(mass: f64, potential: Vec<f64>) -> Self {
    Self { mass, potential }
  }

  /// The potential read off a profile at the sample points.
  pub fn sampled(grid: Grid<Position>, mass: f64, profile: impl Fn(f64) -> f64) -> Self {
    Self::new(mass, grid.coordinates().map(profile).collect())
  }

  pub fn free(grid: Grid<Position>, mass: f64) -> Self {
    Self::new(mass, vec![0.0; grid.npoints])
  }

  pub fn npoints(&self) -> usize {
    self.potential.len()
  }

  /// `t = planck^2 / (2 m (dif x)^2)`, the amplitude for a step between
  /// neighboring points, in which the whole discretization sits.
  pub fn hopping(&self, grid: Grid<Position>) -> f64 {
    HBAR.powi(2) / (2.0 * self.mass * grid.spacing.powi(2))
  }

  /// `(hat(H) psi)_j = (2 t + V_j) psi_j - t psi_(j+1) - t psi_(j-1)`
  pub fn apply(&self, state: &mut GridState<Position>, grid: Grid<Position>) {
    assert_eq!(grid.npoints, self.npoints());
    assert_eq!(grid.npoints, state.npoints());
    let hopping = self.hopping(grid);
    let npoints = self.npoints();
    let before = std::mem::take(&mut state.values);
    state.values = (0..npoints)
      .map(|index| {
        let left = before[(index + npoints - 1) % npoints];
        let right = before[(index + 1) % npoints];
        (2.0 * hopping + self.potential[index]) * before[index] - hopping * (left + right)
      })
      .collect();
  }
  pub fn applied(
    &self,
    mut state: GridState<Position>,
    grid: Grid<Position>,
  ) -> GridState<Position> {
    self.apply(&mut state, grid);
    state
  }

  /// The operator written out, one column per point, each column being the
  /// stencil applied to the state that is one at that point.
  pub fn matrix(&self, grid: Grid<Position>) -> DMatrix<Complex> {
    let npoints = self.npoints();
    DMatrix::from_columns(
      &(0..npoints)
        .map(|column| {
          let mut unit = vec![Complex::ZERO; npoints];
          unit[column] = Complex::ONE;
          let acted = self.applied(GridState::new(unit), grid);
          nalgebra::DVector::from_vec(acted.values)
        })
        .collect::<Vec<_>>(),
    )
  }

  /// The eigenstates, in increasing energy.
  ///
  /// The operator is real and symmetric, a real potential coupling each point
  /// to its neighbors by the same real amplitude in both directions, so the
  /// energies are real and the states can be chosen real.
  pub fn eigenstates(&self, grid: Grid<Position>) -> Vec<(f64, GridState<Position>)> {
    let complex = self.matrix(grid);
    assert!(complex.iter().all(|entry| entry.im.abs() < f64::EPSILON));
    let decomposition = SymmetricEigen::new(complex.map(|entry| entry.re));

    let mut order: Vec<usize> = (0..self.npoints()).collect();
    order.sort_by(|&left, &right| {
      decomposition.eigenvalues[left].total_cmp(&decomposition.eigenvalues[right])
    });
    order
      .into_iter()
      .map(|index| {
        let state = decomposition
          .eigenvectors
          .column(index)
          .iter()
          .map(|&amplitude| Complex::from(amplitude))
          .collect();
        (decomposition.eigenvalues[index], GridState::new(state))
      })
      .collect()
  }

  /// The energies alone, in increasing order.
  pub fn spectrum(&self, grid: Grid<Position>) -> Vec<f64> {
    self
      .eigenstates(grid)
      .into_iter()
      .map(|(energy, _)| energy)
      .collect()
  }
}
