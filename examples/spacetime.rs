//! A Gaussian packet let go in free space, written out as `psi(t, x)`.

mod npy;

use quantum::{
  Complex,
  discrete::{
    grid::{Grid, Position},
    grid_state::GridState,
  },
  dispersion::Dispersion,
};

fn main() -> std::io::Result<()> {
  let grid: Grid<Position> = Grid::new(1024, 0.1);
  let dispersion = Dispersion::electron();

  let start = -30.0;
  let width = 2.0f64;
  let wavenumber = 5.0;
  let mut state = GridState::sample(grid, |position| {
    let offset = position - start;
    Complex::from_polar(
      (-offset.powi(2) / (4.0 * width.powi(2))).exp(),
      wavenumber * position,
    )
  });

  let nframes = 240;
  let step = 0.3;
  let mut frames = Vec::with_capacity(nframes * grid.npoints);
  for _ in 0..nframes {
    frames.extend_from_slice(&state.values);
    state.evolve(step, dispersion, grid);
  }

  npy::write_complexes("out/psi.npy", &[nframes, grid.npoints], &frames)?;
  npy::write_reals(
    "out/position.npy",
    &[grid.npoints],
    &grid.coordinates().collect::<Vec<f64>>(),
  )?;
  npy::write_reals(
    "out/time.npy",
    &[nframes],
    &(0..nframes)
      .map(|frame| frame as f64 * step)
      .collect::<Vec<f64>>(),
  )?;

  println!(
    "group velocity {:.3} nm/fs over {:.0} fs",
    dispersion.group_velocity(wavenumber),
    nframes as f64 * step
  );
  Ok(())
}
