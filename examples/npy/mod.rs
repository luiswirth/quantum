//! Writing arrays in numpy's own format, which needs no dependency.

use quantum::Complex;
use std::io::{Result, Write};
use std::path::Path;

pub fn write_reals(path: impl AsRef<Path>, shape: &[usize], values: &[f64]) -> Result<()> {
  write(
    path,
    "<f8",
    shape,
    values.iter().flat_map(|value| value.to_le_bytes()),
  )
}

pub fn write_complexes(path: impl AsRef<Path>, shape: &[usize], values: &[Complex]) -> Result<()> {
  write(
    path,
    "<c16",
    shape,
    values
      .iter()
      .flat_map(|value| [value.re, value.im])
      .flat_map(|part| part.to_le_bytes()),
  )
}

fn write(
  path: impl AsRef<Path>,
  descr: &str,
  shape: &[usize],
  bytes: impl Iterator<Item = u8>,
) -> Result<()> {
  let dimensions: Vec<String> = shape.iter().map(|length| format!("{length},")).collect();
  let header = format!(
    "{{'descr': '{descr}', 'fortran_order': False, 'shape': ({}), }}",
    dimensions.join("")
  );
  // The header is padded so that the data begins on a 64 byte boundary.
  let padding = 64 - (10 + header.len() + 1) % 64;
  let header = format!("{header}{}\n", " ".repeat(padding));

  if let Some(parent) = path.as_ref().parent() {
    std::fs::create_dir_all(parent)?;
  }
  let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
  file.write_all(b"\x93NUMPY\x01\x00")?;
  file.write_all(&(header.len() as u16).to_le_bytes())?;
  file.write_all(header.as_bytes())?;
  for byte in bytes {
    file.write_all(&[byte])?;
  }
  file.flush()
}
