//! This crate provides an implementation of:
//!
//! * [Owen's T function](https://en.wikipedia.org/wiki/Owen%27s_T_function)
//! * The bivariate normal CDF function (using Owen's T)

#![cfg_attr(not(feature = "std"), no_std)]

mod bivariate_normal;
mod owens_t;
mod util;

pub use bivariate_normal::{BivNormArg, BivNormRho, biv_norm, biv_norm_inner};
pub use owens_t::{owens_t, owens_t_inner};
