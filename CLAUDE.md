# quantum

Elementary quantum mechanics encoded in Rust, written to be read: each type is
an object of the theory, and each test is a statement about it rather than a
regression guard. Public at <https://github.com/luiswirth/quantum>.

The toolchain comes from the flake, which also installs the pre-commit hook.

## Layout

The exact and the discrete pictures are separate module trees, and what both
consult sits above them.

- `exact/`:
  the states on the whole line. `combination.rs` is `Weighted<B>` and
  `Combination<B>`, the free complex vector space over a basis type, carrying
  everything that holds whatever `B` is; `plane_wave.rs` and `dirac.rs` are the
  two dual bases, each in a space-only and a spacetime form, with the aliases
  naming their weighted and combined forms; `superposition.rs` is the plane
  wave combination.
- `discrete/`:
  the finite dimensional picture. `grid.rs` is the sample points and their
  dual, `grid_state.rs` the values on them and the transform between the two
  domains.
- `sampling.rs`:
  the bridge, in both directions, and the only module naming both sides.
- `dispersion.rs`:
  the relation the medium ties frequency to wavenumber by, and the only place a
  medium enters.
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
