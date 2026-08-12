use std::f64::consts::TAU;
use std::marker::PhantomData;

/// One side of the transform, whose dual is the other.
pub trait Domain {
  type Dual: Domain<Dual = Self>;
  /// The sign of `s` in `e^(i s k x)` when leaving this domain.
  const SIGN: f64;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Position;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Momentum;

impl Domain for Position {
  type Dual = Momentum;
  const SIGN: f64 = -1.0;
}
impl Domain for Momentum {
  type Dual = Position;
  const SIGN: f64 = 1.0;
}

/// `N` samples of spacing `dif x`, periodic over `L = N dif x`.
///
/// The samples are indexed by `ZZ_N`, whose identity `0` is the origin and
/// whose upper half are the negative residues, so that index arithmetic is
/// already the grid's own translation.
#[derive(Debug, PartialEq)]
pub struct Grid<D: Domain> {
  pub npoints: usize,
  pub spacing: f64,
  domain: PhantomData<D>,
}

impl<D: Domain> Clone for Grid<D> {
  fn clone(&self) -> Self {
    *self
  }
}
impl<D: Domain> Copy for Grid<D> {}

impl<D: Domain> Grid<D> {
  pub fn new(npoints: usize, spacing: f64) -> Self {
    assert!(npoints > 0);
    assert!(spacing > 0.0);
    Self {
      npoints,
      spacing,
      domain: PhantomData,
    }
  }

  /// `L = N dif x`
  pub fn length(&self) -> f64 {
    self.npoints as f64 * self.spacing
  }

  pub fn offset(&self, index: usize) -> i64 {
    assert!(index < self.npoints);
    if 2 * index < self.npoints {
      index as i64
    } else {
      index as i64 - self.npoints as i64
    }
  }

  /// `x_j = j dif x`
  pub fn coordinate(&self, index: usize) -> f64 {
    self.offset(index) as f64 * self.spacing
  }
  pub fn coordinates(&self) -> impl Iterator<Item = f64> + '_ {
    (0..self.npoints).map(|index| self.coordinate(index))
  }

  /// The grid of the other domain, on which discreteness and periodicity trade
  /// places: `dif k = (2 pi) / L` and `L_k = (2 pi) / (dif x)`.
  pub fn dual(&self) -> Grid<D::Dual> {
    Grid::new(self.npoints, TAU / self.length())
  }
}
