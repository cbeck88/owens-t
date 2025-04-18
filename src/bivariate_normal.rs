use crate::{owens_t_inner, util::*};
use libm::{asin, erf};

// Compute normal CDF using erf
fn phi_1(x: f64) -> f64 {
    0.5 * (1.0 + erf(x * FRAC_1_SQRT_2))
}

/// Compute bivariate normal CDF, using owens' t.
///
/// This is from Owen, see section 2 of: https://www.scirp.org/journal/paperinformation?paperid=128377
pub fn biv_norm(x: f64, y: f64, rho: f64) -> f64 {
    biv_norm_inner(x.into(), y.into(), rho.into())
}

/// Argument to `biv_norm` which will be used repeatedly, and caches some re-used values.
#[derive(Copy, Clone, Default, Debug)]
pub struct BivNormArg {
    val: f64,
    val_recip: f64,
    phi_1_val: f64,
}

impl From<f64> for BivNormArg {
    fn from(src: f64) -> Self {
        Self {
            phi_1_val: phi_1(src),
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
        phi_1_val: phi_1_x,
    }: BivNormArg,
    BivNormArg {
        val: y,
        val_recip: y_recip,
        phi_1_val: phi_1_y,
    }: BivNormArg,
    BivNormRho {
        rho,
        sqrt_1_minus_rho_sq_recip,
    }: BivNormRho,
) -> f64 {
    debug_assert!(-1.0 <= rho);
    debug_assert!(1.0 >= rho);

    if x == 0.0 && y == 0.0 {
        return 0.25 + ONE_DIV_TWO_PI * asin(rho);
    }

    if rho == 0.0 {
        return phi_1_x * phi_1_y;
    } else if rho == 1.0 {
        // return phi_1(f64::min(x, y));
        return if x < y { phi_1_x } else { phi_1_y };
    } else if rho == -1.0 {
        return f64::max(phi_1_x + phi_1_y - 1.0, 0.0);
    };

    // Nonzero
    // let sqrt_1_minus_rho_sq = (1.0 - rho * rho).sqrt();

    let x_contrib = if x == 0.0 {
        0.0
    } else {
        //let r_x = (y - rho * x) / (x * sqrt_1_minus_rho_sq);
        //(0.5 * phi_1(x)) - owens_t(x, r_x)
        let r_x = (y * x_recip - rho) * sqrt_1_minus_rho_sq_recip;
        (0.5 * phi_1_x) - owens_t_inner(x, r_x, Some((phi_1_x-0.5).abs()))
    };
    let y_contrib = if y == 0.0 {
        0.0
    } else {
        //let r_y = (x - rho * y) / (y * sqrt_1_minus_rho_sq);
        //(0.5 * phi_1(y)) - owens_t(y, r_y)
        let r_y = (x * y_recip - rho) * sqrt_1_minus_rho_sq_recip;
        (0.5 * phi_1_y) - owens_t_inner(y, r_y, Some((phi_1_y-0.5).abs()))
    };

    // I know you dont like it clippy, but I just want to follow the paper.
    #[allow(clippy::if_same_then_else)]
    let beta = if (x * y) > 0.0 {
        0.0
    } else if (x * y) == 0.0 && x + y >= 0.0 {
        0.0
    } else {
        0.5
    };

    //0.5 * (phi(x) + phi(y)) - owens_t(x, r_x) - owens_t(y, r_y) - beta
    x_contrib + y_contrib - beta
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
        // FIXME: Precision should be a little better than this, probably we should be more careful with
        // erf vs. erfc
        let eps = 0.000001;
        for n in 0..N_MAX {
            let x = X_VEC[n];
            let y = Y_VEC[n];
            let r = R_VEC[n];
            let fxy = FXY_VEC[n];
            let val = biv_norm(x, y, r);
            //eprintln!("n = {n}: biv_norm({x}, {y}, {r}) = {val}: expected: {fxy}");
            assert_within!(~eps, biv_norm(y,x,r), val);
            assert_within!(~eps, val, fxy, "n = {n}, x = {x}, y = {y}, rho = {r}")
        }
    }
}
