# quantum

Elementary quantum mechanics encoded in Rust, written to be read: each type is
an object of the theory, and each test is a statement about it rather than a
regression guard. Public at <https://github.com/luiswirth/quantum>.

The toolchain comes from the flake, so `nix develop --command cargo test` runs
from anywhere.

## Layout

Two dual bases, each in a space-only and a spacetime form, over one vector
space.

- `combination.rs`:
  `Weighted<B>` and `Combination<B>`, the free complex vector space over a
  basis type, carrying everything that holds whatever `B` is.
- `plane_wave.rs`, `dirac.rs`:
  the two bases and the aliases naming their weighted and combined forms.
- `dispersion.rs`:
  the relation the medium ties frequency to wavenumber by, and the only place a
  medium enters.
- `grid.rs`, `grid_state.rs`:
  the finite dimensional picture, and the transform between the two domains.
- `superposition.rs`:
  the plane wave combination, and what joins the exact objects to the grid.
- `units.rs`:
  the unit system every number is expressed in.

## Conventions

- The generic layer never learns which basis it serves.
  Anything naming a basis lives in that basis's module, which is also where an
  operator producing a weighted or combined form of it goes.
- Exact and discrete coexist, and neither replaces the other.
  The combinations are exact on the whole line, the grid types are the finite
  picture, and each bridge loses exactly what the transform says it must.
- Amplitudes are never normalized.
  A normalizing constructor would break the vector space, so every expectation
  divides by the total intensity where it is taken.
- Storage is in the grid's own index order.
  The monotone reordering exists for a reader and belongs to whatever draws the
  picture.
- A trait appears only where a closed duality needs one, as `Domain` does.
  Everything else is data and inherent impls.
