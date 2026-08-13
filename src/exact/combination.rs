use crate::Complex;
use std::ops::{Add, Mul};

/// A basis element with a complex weight.
#[derive(Clone, Copy, Debug)]
pub struct Weighted<B> {
  pub amplitude: Complex,
  pub element: B,
}

/// The free complex vector space over the basis elements.
#[derive(Clone, Debug)]
pub struct Combination<B> {
  pub terms: Vec<Weighted<B>>,
}

impl<B> Weighted<B> {
  pub fn new(amplitude: Complex, element: B) -> Self {
    Self { amplitude, element }
  }

  /// `abs(c)^2`
  pub fn intensity(&self) -> f64 {
    self.amplitude.norm_sqr()
  }
}

impl<B> Combination<B> {
  pub fn new(terms: Vec<Weighted<B>>) -> Self {
    Self { terms }
  }

  pub fn total_intensity(&self) -> f64 {
    self.terms.iter().map(Weighted::intensity).sum()
  }
}

impl<B> From<Weighted<B>> for Combination<B> {
  fn from(term: Weighted<B>) -> Self {
    Combination::new(vec![term])
  }
}

impl<B> Mul<Complex> for Weighted<B> {
  type Output = Weighted<B>;
  fn mul(self, factor: Complex) -> Weighted<B> {
    Weighted::new(self.amplitude * factor, self.element)
  }
}

impl<B> Mul<Complex> for Combination<B> {
  type Output = Combination<B>;
  fn mul(self, factor: Complex) -> Combination<B> {
    Combination::new(self.terms.into_iter().map(|term| term * factor).collect())
  }
}

impl<B> Add for Weighted<B> {
  type Output = Combination<B>;
  fn add(self, other: Weighted<B>) -> Combination<B> {
    Combination::new(vec![self, other])
  }
}

impl<B> Add for Combination<B> {
  type Output = Combination<B>;
  fn add(mut self, other: Combination<B>) -> Combination<B> {
    self.terms.extend(other.terms);
    self
  }
}
