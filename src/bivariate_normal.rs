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
/// Accurate to ~15 decimals for x * y > 0.0
/// Only 6 or 7 decimals when x * y < 0.0
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
    } else if
    /*(x == 0.0 || y== 0.0) &&*/
    x + y >= 0.0 {
        -0.5
    } else {
        0.0
    };

    // The correction factors USUALLY but don't always, add up to 0.0 or 0.5
    //debug_assert!(c_x + c_y - beta == 0.5 || c_x + c_y - beta == 0.0);
    // If x == 0 or y == 0, then the correction factors add up to 0.0
    /*debug_assert!(
        x * y != 0.0 || c_x + c_y - beta == 0.0,
        "{x}, {y}, {c_x}, {c_y}, {beta}"
    );*/

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
        // This is also valid for a > 0 and h,
        // but see Owens (3.5) for the correction factor, when a is negative:
        //
        // Q(h, a) = T(ah, 1/a) - (Phi(h) - 1/2)(1- Phi(ah)) + 1/2
        //
        // We then have:
        //
        // T(ah, 1/a)
        // = signum(a) * T(|ah|, |1/a|)
        let phi_h_minus_half = owens_t_znorm1(h);
        let one_minus_phi_ah = owens_t_znorm2(a * h);
        let one_half_if_a_negative = (s_a - 1.0) * -0.25;
        (
            (s_a * owens_t_dispatch(ah_abs, a_abs.recip(), h_abs, Some(phi_h_minus_half.abs())))
                - phi_h_minus_half * one_minus_phi_ah,
            one_half_if_a_negative,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_within::assert_within;

    // These values from: https://people.math.sc.edu/Burkardt/cpp_src/test_values/test_values.cpp
    // who says they come from Mathematica
    //
    //void bivariate_normal_cdf_values ( int &n_data, double &x, double &y,
    //  double &r, double &fxy )
    //
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
    //
    //  Discussion:
    //
    //    FXY is the probability that two variables A and B, which are
    //    related by a bivariate normal distribution with correlation R,
    //    respectively satisfy A <= X and B <= Y.
    //
    //    Mathematica can evaluate the bivariate normal CDF via the commands:
    //
    //      <<MultivariateStatistics`
    //      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    23 November 2010
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    National Bureau of Standards,
    //    Tables of the Bivariate Normal Distribution and Related Functions,
    //    NBS, Applied Mathematics Series, Number 50, 1959.
    //
    //  Parameters:
    //
    //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
    //    first call.  On each call, the routine increments N_DATA by 1, and
    //    returns the corresponding data; when there is no more data, the
    //    output value of N_DATA will be 0 again.
    //
    //    Output, double &X, &Y, the parameters of the function.
    //
    //    Output, double &R, the correlation value.
    //
    //    Output, double &FXY, the value of the function.
    //

    const N_MAX: usize = 41;
    const FXY_VEC: [f64; N_MAX] = [
        0.02260327218569867E+00,
        0.1548729518584100E+00,
        0.4687428083352184E+00,
        0.7452035868929476E+00,
        0.8318608306874188E+00,
        0.8410314261134202E+00,
        0.1377019384919464E+00,
        0.1621749501739030E+00,
        0.1827411243233119E+00,
        0.2010067421506235E+00,
        0.2177751155265290E+00,
        0.2335088436446962E+00,
        0.2485057781834286E+00,
        0.2629747825154868E+00,
        0.2770729823404738E+00,
        0.2909261168683812E+00,
        0.3046406378726738E+00,
        0.3183113449213638E+00,
        0.3320262544108028E+00,
        0.3458686754647614E+00,
        0.3599150462310668E+00,
        0.3742210899871168E+00,
        0.3887706405282320E+00,
        0.4032765198361344E+00,
        0.4162100291953678E+00,
        0.6508271498838664E+00,
        0.8318608306874188E+00,
        0.0000000000000000,
        0.1666666666539970,
        0.2500000000000000,
        0.3333333333328906,
        0.5000000000000000,
        0.7452035868929476,
        0.1548729518584100,
        0.1548729518584100,
        0.06251409470431653,
        0.7452035868929476,
        0.1548729518584100,
        0.1548729518584100,
        0.06251409470431653,
        0.6337020457912916,
    ];
    const R_VEC: [f64; N_MAX] = [
        0.500, 0.500, 0.500, 0.500, 0.500, 0.500, -0.900, -0.800, -0.700, -0.600, -0.500, -0.400,
        -0.300, -0.200, -0.100, 0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800,
        0.900, 0.673, 0.500, -1.000, -0.500, 0.000, 0.500, 1.000, 0.500, 0.500, 0.500, 0.500,
        0.500, 0.500, 0.500, 0.500, 0.500,
    ];
    const X_VEC: [f64; N_MAX] = [
        -2.0,
        -1.0,
        0.0,
        1.0,
        2.0,
        3.0,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        -0.2,
        1.0,
        2.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        -1.0,
        -1.0,
        1.0,
        1.0,
        -1.0,
        -1.0,
        0.7071067811865475,
    ];
    const Y_VEC: [f64; N_MAX] = [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        -1.0,
        1.0,
        -1.0,
        1.0,
        -1.0,
        1.0,
        -1.0,
        0.7071067811865475,
    ];

    #[test]
    fn spot_check_phi2() {
        for n in 0..N_MAX {
            // These are negated, see documentation of these test points
            let x = -X_VEC[n];
            let y = -Y_VEC[n];
            let r = R_VEC[n];
            let expected = FXY_VEC[n];
            let val = biv_norm(x, y, r);

            // FIXME: Precision should be a little better than this...
            // I strongly suspect these test vectors, since tvpack agrees with this one.
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
            assert_within!(~eps, biv_norm(y,x,r), val);
            assert_within!(~eps, val, expected, "n = {n}, x = {x}, y = {y}, rho = {r}")
        }
    }
}
