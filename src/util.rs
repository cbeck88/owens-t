//! Constants and polyfills for functions like f64::exp and f64::sqrt which
//! dissappear on no-std mode.

pub use core::f64::consts::{FRAC_1_PI, FRAC_1_SQRT_2, FRAC_2_SQRT_PI};

pub const ONE_DIV_TWO_PI: f64 = FRAC_1_PI * 0.5;
pub const ONE_DIV_ROOT_TWO_PI: f64 = FRAC_2_SQRT_PI * 0.5 * FRAC_1_SQRT_2;

// When std is available, the built-in f64::exp uses intrinsics and is like 5 nanos
#[cfg(feature = "std")]
#[inline(always)]
pub fn exp(x: f64) -> f64 {
    x.exp()
}

#[cfg(not(feature = "std"))]
#[inline(always)]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

#[cfg(feature = "std")]
#[inline(always)]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

#[cfg(not(feature = "std"))]
#[inline(always)]
pub fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}
