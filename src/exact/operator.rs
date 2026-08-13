use crate::{Complex, exact::combination::Combination};

/// An operator the basis `B` diagonalizes, which is to say the basis of its own
/// eigenvectors.
///
/// Diagonality is relative to a basis and is what the parameter records. The
/// operator acts by scaling each amplitude by the eigenvalue of the element it
/// weights, and that is closed on a finite combination only here. In any other
/// basis the same operator differentiates instead, and leaves the span.
pub struct DiagonalOperator<B> {
  eigenvalue: Box<dyn Fn(&B) -> Complex>,
}

impl<B: 'static> DiagonalOperator<B> {
  pub fn new(eigenvalue: impl Fn(&B) -> Complex + 'static) -> Self {
    Self {
      eigenvalue: Box::new(eigenvalue),
    }
  }

  /// `hat(1)`, which every basis diagonalizes.
  pub fn identity() -> Self {
    Self::new(|_| Complex::ONE)
  }

  /// `a_n`
  pub fn eigenvalue(&self, element: &B) -> Complex {
    (self.eigenvalue)(element)
  }

  /// `f(hat(A))`, which the same basis diagonalizes, its eigenvalues carried
  /// through `f`.
  pub fn map(&mut self, function: impl Fn(Complex) -> Complex + 'static) {
    let inner = std::mem::replace(&mut self.eigenvalue, Box::new(|_| Complex::ZERO));
    self.eigenvalue = Box::new(move |element| function(inner(element)));
  }
  pub fn mapped(mut self, function: impl Fn(Complex) -> Complex + 'static) -> Self {
    self.map(function);
    self
  }

  /// `hat(A) psi = sum_n a_n c_n b_n`
  pub fn applied(&self, state: &Combination<B>) -> Combination<B>
  where
    B: Clone,
  {
    Combination::new(
      state
        .terms
        .iter()
        .map(|term| term.clone() * self.eigenvalue(&term.element))
        .collect(),
    )
  }

  /// `expect(hat(A)) = (sum_n abs(c_n)^2 a_n) / (sum_n abs(c_n)^2)`
  pub fn expectation(&self, state: &Combination<B>) -> Complex {
    let weighted: Complex = state
      .terms
      .iter()
      .map(|term| self.eigenvalue(&term.element) * term.intensity())
      .sum();
    weighted / state.total_intensity()
  }

  /// `variance(hat(A)) = expect(hat(A)^2) - expect(hat(A))^2`
  pub fn variance(&self, state: &Combination<B>) -> Complex {
    let squared: Complex = state
      .terms
      .iter()
      .map(|term| self.eigenvalue(&term.element).powi(2) * term.intensity())
      .sum();
    squared / state.total_intensity() - self.expectation(state).powi(2)
  }

  /// `Delta A = sqrt(variance(hat(A)))`, which a Hermitian operator leaves
  /// real.
  pub fn uncertainty(&self, state: &Combination<B>) -> f64 {
    self.variance(state).re.sqrt()
  }
}
