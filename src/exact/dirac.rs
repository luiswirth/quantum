use crate::{
  Complex,
  exact::{
    combination::{Combination, Weighted},
    operator::DiagonalOperator,
  },
};
use std::ops::Mul;

/// `delta(x - x_0)`, which has no value anywhere and is known only by what it
/// takes from a test function.
#[derive(Clone, Copy, Debug)]
pub struct Dirac {
  pub position: f64,
}

impl Dirac {
  pub fn new(position: f64) -> Self {
    Self { position }
  }

  /// `integral delta(x - x_0) f(x) dif x = f(x_0)`
  pub fn apply(&self, test: impl Fn(f64) -> Complex) -> Complex {
    test(self.position)
  }

  /// `e^(-i k x_0)`
  pub fn transform_at(&self, wavenumber: f64) -> Complex {
    Complex::from_polar(1.0, -wavenumber * self.position)
  }
}

/// `delta(t - t_0) delta(x - x_0)`, a point of spacetime.
#[derive(Clone, Copy, Debug)]
pub struct Event {
  pub time: f64,
  pub dirac: Dirac,
}

impl Event {
  pub fn new(time: f64, dirac: Dirac) -> Self {
    Self { time, dirac }
  }

  /// `f(t_0, x_0)`
  pub fn apply(&self, test: impl Fn(f64, f64) -> Complex) -> Complex {
    test(self.time, self.dirac.position)
  }

  /// `e^(-i (k x_0 - omega t_0))`
  pub fn transform_at(&self, frequency: f64, wavenumber: f64) -> Complex {
    self.dirac.transform_at(wavenumber) * Complex::from_polar(1.0, frequency * self.time)
  }
}

/// `c delta(x - x_0)`
pub type Impulse = Weighted<Dirac>;

impl Impulse {
  pub fn apply(&self, test: impl Fn(f64) -> Complex) -> Complex {
    self.amplitude * self.element.apply(test)
  }

  pub fn transform_at(&self, wavenumber: f64) -> Complex {
    self.amplitude * self.element.transform_at(wavenumber)
  }
}

/// `sum_j c_j delta(x - x_j)`
pub type Comb = Combination<Dirac>;

impl Comb {
  /// `sum_j c_j f(x_j)`
  pub fn apply(&self, test: impl Fn(f64) -> Complex + Copy) -> Complex {
    self.terms.iter().map(|term| term.apply(test)).sum()
  }

  /// `sum_j c_j e^(-i k x_j)`, an ordinary function of the wavenumber.
  pub fn transform_at(&self, wavenumber: f64) -> Complex {
    self
      .terms
      .iter()
      .map(|term| term.transform_at(wavenumber))
      .sum()
  }
}

/// `c delta(t - t_0) delta(x - x_0)`
pub type Source = Weighted<Event>;

impl Source {
  pub fn apply(&self, test: impl Fn(f64, f64) -> Complex) -> Complex {
    self.amplitude * self.element.apply(test)
  }

  pub fn transform_at(&self, frequency: f64, wavenumber: f64) -> Complex {
    self.amplitude * self.element.transform_at(frequency, wavenumber)
  }
}

/// `sum_j c_j delta(t - t_j) delta(x - x_j)`
pub type Excitation = Combination<Event>;

impl Excitation {
  pub fn apply(&self, test: impl Fn(f64, f64) -> Complex + Copy) -> Complex {
    self.terms.iter().map(|term| term.apply(test)).sum()
  }

  pub fn transform_at(&self, frequency: f64, wavenumber: f64) -> Complex {
    self
      .terms
      .iter()
      .map(|term| term.transform_at(frequency, wavenumber))
      .sum()
  }
}

impl Mul<Complex> for Dirac {
  type Output = Impulse;
  fn mul(self, amplitude: Complex) -> Impulse {
    Weighted::new(amplitude, self)
  }
}

impl Mul<Complex> for Event {
  type Output = Source;
  fn mul(self, amplitude: Complex) -> Source {
    Weighted::new(amplitude, self)
  }
}

impl DiagonalOperator<Dirac> {
  /// `hat(x) = x_0`
  pub fn position() -> Self {
    Self::new(|dirac: &Dirac| dirac.position.into())
  }

  /// `hat(V) = V(hat(x))`, a function of the position operator and so diagonal
  /// wherever that one is.
  pub fn potential(profile: impl Fn(f64) -> f64 + 'static) -> Self {
    Self::new(move |dirac: &Dirac| profile(dirac.position).into())
  }
}
