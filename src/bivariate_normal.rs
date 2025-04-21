use crate::{
    owens_t::{owens_t_dispatch, owens_t_znorm1, owens_t_znorm2},
    util::*,
};
use libm::{asin, erfc};

// Compute 1 - normal CDF using erf
fn one_minus_phi(x: f64) -> f64 {
    0.5 * erfc(x * FRAC_1_SQRT_2)
}

/// Compute bivariate normal CDF, using Owens' T function.
/// This is Pr[ X > x, Y > y] when X and Y are standard normals of correlation coefficient `rho`.
///
/// Accurate to ~15 decimals.
///
/// Preconditions:
///   -1 <= rho <= 1
pub fn biv_norm(x: f64, y: f64, rho: f64) -> f64 {
    biv_norm_inner(x.into(), y.into(), rho.into())
}

/// Argument to `biv_norm` which will be used repeatedly, and caches some re-used values.
#[derive(Copy, Clone, Default, Debug)]
pub struct BivNormArg {
    val: f64,
    val_recip: f64,
    one_minus_phi_val: f64,
}

impl From<f64> for BivNormArg {
    fn from(src: f64) -> Self {
        Self {
            one_minus_phi_val: one_minus_phi(src),
            val_recip: src.recip(),
            val: src,
        }
    }
}

impl AsRef<f64> for BivNormArg {
    fn as_ref(&self) -> &f64 {
        &self.val
    }
}

/// Rho argument to biv_norm which will be used repeatedly, and caches some re-used values.
#[derive(Copy, Clone, Default, Debug)]
pub struct BivNormRho {
    rho: f64,
    sqrt_1_minus_rho_sq_recip: f64,
}

impl From<f64> for BivNormRho {
    fn from(rho: f64) -> Self {
        Self {
            sqrt_1_minus_rho_sq_recip: sqrt(1.0 - rho * rho).recip(),
            rho,
        }
    }
}

/// Version of biv_norm which allows re-using some computation when the same value of x or y is used repeatedly.
pub fn biv_norm_inner(
    BivNormArg {
        val: x,
        val_recip: x_recip,
        one_minus_phi_val: one_minus_phi_x,
    }: BivNormArg,
    BivNormArg {
        val: y,
        val_recip: y_recip,
        one_minus_phi_val: one_minus_phi_y,
    }: BivNormArg,
    BivNormRho {
        rho,
        sqrt_1_minus_rho_sq_recip,
    }: BivNormRho,
) -> f64 {
    debug_assert!(-1.0 <= rho, "{rho}");
    debug_assert!(1.0 >= rho, "{rho}");

    if x == 0.0 && y == 0.0 {
        return 0.25 + ONE_DIV_TWO_PI * asin(rho);
    }

    if rho == 0.0 {
        return one_minus_phi_x * one_minus_phi_y;
    } else if rho == 1.0 {
        return if x < y {
            one_minus_phi_y
        } else {
            one_minus_phi_x
        };
    } else if rho == -1.0 {
        return f64::max(one_minus_phi_x + one_minus_phi_y - 1.0, 0.0);
    };

    // Nonzero
    // let sqrt_1_minus_rho_sq = (1.0 - rho * rho).sqrt();

    // Here, q_x = Q(x, a_x)
    // In notation of Owen '56, equation 2.1.
    // See also Pages 15-16 of Patefield Tandy.
    //
    // Or see section 2 of: https://www.scirp.org/journal/paperinformation?paperid=128377
    // although beware because they invert the CDF, so x must be negated etc.
    //
    // Read Patefield Tandy closely:
    // > Care should be taken when computing (11) as, although the function T is accurate to at least
    // > 14 significant figures, the subtraction can lead to a loss of relative accuracy in computing Q
    // > and hence in the resultant bivariate normal probability
    //
    // To do what they are saying, we don't call `owens_t` as a black-box, instead we build a function that
    // computes Q, which calls owens_t_dispatch directly.
    let (q_x, c_x) = if x == 0.0 {
        // We get an indeterminate value for r_x when x == 0.
        // However, Owen's T is defined even at infinity, and
        // Owen gives T(h, infinity) = 1/2 (1 - Phi(h)) if h >= 0, 1/2 Phi(h) if h <= 0.
        //
        // We also have T(h, a) = T(h, -a) always.
        //
        // So we can evaluate:
        // Q(h, infinity) := 1/2 [1 - Phi(h)] - T(h, infinity)
        // = 0 if h >= 0
        //   1/2 - Phi(h) if h < 0
        //
        // Since h = x = 0, either way the result is zero.
        (0.0, 0.0)
    } else {
        //let r_x = (y - rho * x) / (x * sqrt_1_minus_rho_sq);
        let r_x = (y * x_recip - rho) * sqrt_1_minus_rho_sq_recip;
        // znorm1 = 0.5 erf (x /sqrt(2)). Also erf is odd, so |znorm1(x)| = znorm1(|x|)
        // let znorm1_abs_x = (0.5 - one_minus_phi_x).abs();
        //0.5 * one_minus_phi_x - owens_t_inner(x, r_x, None /*Some(znorm1_abs_x)*/)

        q(x, r_x, one_minus_phi_x)
    };
    let (q_y, c_y) = if y == 0.0 {
        (0.0, 0.0)
    } else {
        //let r_y = (x - rho * y) / (y * sqrt_1_minus_rho_sq);
        let r_y = (x * y_recip - rho) * sqrt_1_minus_rho_sq_recip;
        // znorm1 = 0.5 erf (x /sqrt(2)). Also erf is odd, so |znorm1(y)| = znorm1(|y|)
        // let znorm1_abs_y = (0.5 - one_minus_phi_y).abs();
        //0.5 * one_minus_phi_y - owens_t_inner(y, r_y, None /*Some(znorm1_abs_y)*/)

        q(y, r_y, one_minus_phi_y)
    };

    let beta = if (x * y) > 0.0 {
        0.0
    } else if (x * y) < 0.0 {
        0.5
    } else {
        0.0
    };

    q_x + q_y + (c_x + c_y - beta)
}

// Compute Q(h, a), described in Owens '54 and in Patefield-Tandy.
// See also their Q routine, which is limited to h >= 0, a >= 0
//
// Ours is different in that, the possible 1/2 correction factor is held on the side
// to try to prevent loss of precision. It's not the only term that can have a nasty
// cancellation though.
fn q(h: f64, a: f64, one_minus_phi_h: f64) -> (f64, f64) {
    let s_a = a.signum();
    let h_abs = h.abs();
    let a_abs = a.abs();
    let ah_abs = h_abs * a_abs;

    if a_abs <= 1.0 {
        // Patefield-Tandy eq (11), citing Owens (2.1):
        //
        // Q(h, a) := 1/2 [1 - Phi(h)] - T(h, a)
        //
        // Their equation is valid for all h and a, but owens_t_dispatch has requirements.
        // We're using T(-h, a) = T(h, a)
        // and T(h, -a) = -T(h, a)
        (
            0.5 * one_minus_phi_h - (s_a * owens_t_dispatch(h_abs, a_abs, ah_abs, None)),
            0.0,
        )
    } else {
        // Patefield-Tandy page 16, citing Owens (2.3):
        //
        // Q(h, a) = T(ah, 1/a) - (Phi(h) - 1/2)(1- Phi(ah))
        //
        // This is valid for a > 0 and h,
        // but see Owens (3.5) for the correction factor, when a is negative:
        //
        // Q(h, a) = T(ah, 1/a) - (Phi(h) - 1/2)(1- Phi(ah)) + 1/2
        //
        // We then have:
        //
        // T(ah, 1/a)
        // = signum(a) * T(|ah|, |1/a|)
        //
        // For large h, we can use erf(h/sqrt(2)) = 1 - erfc(h/sqrt(2))
        // and won't lose precision, since erf will be large anyways, and save some cycles.
        // For small h, we should recompute erf to get the low order decimals right.
        // This is the same threshold for h that boost::math owens impl was using.
        let phi_h_minus_half = if h_abs <= 0.67 {
            owens_t_znorm1(h)
        } else {
            0.5 - one_minus_phi_h
        };
        let one_minus_phi_ah = owens_t_znorm2(a * h);
        let one_half_if_a_negative = (s_a - 1.0) * -0.25;
        (
            // Note: owens_t_znorm1 = erf(...) is an odd function, and we have |erf(x)| = erf(|x|).
            (s_a * owens_t_dispatch(ah_abs, a_abs.recip(), h_abs, Some(phi_h_minus_half.abs())))
                - phi_h_minus_half * one_minus_phi_ah,
            one_half_if_a_negative,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::{BvndTestPoint, get_axis_test_points, get_burkardt_nbs_test_points};
    use assert_within::assert_within;
    use rand::Rng;
    use rand_pcg::{Pcg64Mcg, rand_core::SeedableRng};

    #[test]
    fn spot_check_biv_norm_against_axis_points() {
        for (n, BvndTestPoint { x, y, r, expected }) in get_axis_test_points().enumerate() {
            let eps = 1e-15;

            let val = biv_norm(x, y, r);
            //eprintln!("n = {n}: biv_norm({x}, {y}, {r}) = {val}: expected: {fxy}");
            assert_within!(+eps, biv_norm(y, x, r), val, "n = {n}, x = {x}, y = {y}, rho = {r}");
            assert_within!(+eps, val, expected, "n = {n}, x = {x}, y = {y}, rho = {r}")
        }
    }

    #[test]
    fn spot_check_biv_norm_against_burkardt_nbs_points() {
        for (n, BvndTestPoint { x, y, r, expected }) in get_burkardt_nbs_test_points().enumerate() {
            let val = biv_norm(x, y, r);

            // FIXME: Precision should be a little better than this...
            // I strongly suspect these test vectors, since tvpack agrees with this one to 15 decimals.
            let eps = if x > 0.0 && y > 0.0 {
                1e-10
            } else if x == 0.0 || y == 0.0 {
                1e-6
            } else if x < 0.0 && y < 0.0 {
                1e-9
            } else {
                1e-7
            };
            //eprintln!("n = {n}: biv_norm({x}, {y}, {r}) = {val}: expected: {fxy}");
            assert_within!(+eps, biv_norm(y,x,r), val);
            assert_within!(+eps, val, expected, "n = {n}, x = {x}, y = {y}, rho = {r}")
        }
    }

    fn to_three_decimals(x: f64) -> f64 {
        (x * 1000.0).round() / 1000.0
    }

    #[test]
    fn check_symmetry_conditions() {
        let mut rng = Pcg64Mcg::seed_from_u64(9);

        for n in 0..10000 {
            let x = to_three_decimals(4.0 * rng.random::<f64>() - 2.0);
            let y = to_three_decimals(4.0 * rng.random::<f64>() - 2.0);
            let r = to_three_decimals(rng.random::<f64>());

            let val = biv_norm(x, y, r);
            // Phi_2(x,y,r) = Phi_2(y,x,r);
            let eps = 1e-15;
            assert_within!(+eps, val, biv_norm(y,x,r), "n = {n}, x = {x}, y = {y}, rho = {r}");

            // Pr[ X > x, Y > y ] = Pr[X > x] - Pr[ X > x, Y < y ]
            // { Y < y } iff { -Y > -y }, and correlation of X and -Y  is -rho.
            let eps = 1e-15;
            assert_within!(+eps, val, one_minus_phi(x) - biv_norm(x,-y,-r), "n = {n}, x = {x}, y = {y}, rho = {r}");
            let eps = 1e-15;
            assert_within!(+eps, val, one_minus_phi(x) - one_minus_phi(-y) + biv_norm(-x,-y,r), "n = {n}, x = {x}, y = {y}, rho = {r}");
        }
    }

    #[test]
    fn check_symmetry_conditions_wider_range() {
        let mut rng = Pcg64Mcg::seed_from_u64(9);

        for n in 0..10000 {
            let x = to_three_decimals(8.0 * rng.random::<f64>() - 4.0);
            let y = to_three_decimals(8.0 * rng.random::<f64>() - 4.0);
            let r = to_three_decimals(rng.random::<f64>());

            let val = biv_norm(x, y, r);
            // Phi_2(x,y,r) = Phi_2(y,x,r);
            let eps = 1e-15;
            assert_within!(+eps, val, biv_norm(y,x,r), "n = {n}, x = {x}, y = {y}, rho = {r}");

            // Pr[ X > x, Y > y ] = Pr[X > x] - Pr[ X > x, Y < y ]
            // { Y < y } iff { -Y > -y }, and correlation of X and -Y  is -rho.
            let eps = 1e-15;
            assert_within!(+eps, val, one_minus_phi(x) - biv_norm(x,-y,-r), "n = {n}, x = {x}, y = {y}, rho = {r}");
            let eps = 1e-15;
            assert_within!(+eps, val, one_minus_phi(x) - one_minus_phi(-y) + biv_norm(-x,-y,r), "n = {n}, x = {x}, y = {y}, rho = {r}");
        }
    }
}
