//! The unit system every number in the crate is expressed in: `eV`, `nm`, `fs`,
//! which leaves a mass in `eV fs^2 / nm^2`.

/// `eV fs`
pub const HBAR: f64 = 0.658_211_956_9;

/// `nm/fs`, exact by the definition of the meter.
pub const LIGHT_SPEED: f64 = 299.792_458;

/// `eV`
pub const ELECTRON_REST_ENERGY: f64 = 510_998.950;

/// `eV fs^2 / nm^2`
pub const ELECTRON_MASS: f64 = ELECTRON_REST_ENERGY / (LIGHT_SPEED * LIGHT_SPEED);
