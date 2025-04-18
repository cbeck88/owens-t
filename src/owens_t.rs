//! Implement owens_t
#![allow(non_snake_case)]
#![allow(clippy::excessive_precision)]

use crate::util::*;
use libm::{atan, atan2, erf, erfc, expm1};

/// Port from lib boost math (https://live.boost.org/doc/libs/1_81_0/boost/math/special_functions/owens_t.hpp)
///
/// Here a is the bounds of integration and h is the parameter to the integrand:
///
/// T(h,a) = 1/2π · ∫₀^a exp(½h²·(1+x²)) / (1+x²) dx
//
// compute Owen's T function, T(h,a), for arbitrary values of h and a
//         template<typename RealType, class Policy>
//         inline RealType owens_t(RealType h, RealType a, const Policy& pol)
pub fn owens_t(h: f64, a: f64) -> f64 {
    owens_t_inner(h, a, None)
}

/// Same as `owens_t`, but if you already have an evaluation of znorm1(|h|), then
/// this will save some cycles in the cases when a >= 1.
///
/// Here znorm1(x) := 0.5 · erf (x/sqrt(2))
/// and `erf` is odd, so |erf(x)| = erf(|x|)
/// If you know Phi(x) or Phi(-x), where Phi is normal CDF, then you have this number
/// at least to some accuracy, and can avoid an extra call to erf here.
pub fn owens_t_inner(mut h: f64, a: f64, znorm1_abs_h: Option<f64>) -> f64 {
    // exploit that T(-h,a) == T(h,a)
    h = h.abs();

    // Use equation (2) in the paper to remap the arguments
    // such that h>=0 and 0<=a<=1 for the call of the actual
    // computation routine.

    let fabs_a = a.abs();
    let fabs_ah = fabs_a * h;

    #[allow(clippy::collapsible_else_if)]
    let val: f64 = if fabs_a <= 1.0 {
        owens_t_dispatch(h, fabs_a, fabs_ah, None)
    }
    // if(fabs_a <= 1.0)
    else {
        if h <= 0.67 {
            let normh = znorm1_abs_h.unwrap_or_else(|| owens_t_znorm1(h));
            let normah = owens_t_znorm1(fabs_ah);
            0.25 - normh * normah - owens_t_dispatch(fabs_ah, fabs_a.recip(), h, Some(normh))
        }
        // if( h <= 0.67 )
        else {
            let normh = owens_t_znorm2(h);
            let normah = owens_t_znorm2(fabs_ah);

            // TODO: Is 0.5 - znorm2 = znorm1 reasonable here?
            // It seems like it should be fine but it's not what tatefield pandy did, maybe it loses some precision.
            //
            // In testing, haven't been able to find precision loss for large a and h when we do this, I think it's fine because
            // even if we lose some precision, this is only used to initialize some of the series computation, so the
            // lower order bits gets lost anyways.
            let znorm1_abs_h = znorm1_abs_h.unwrap_or_else(|| 0.5 - normh);
            //let znorm1_abs_h = znorm1_abs_h.unwrap_or_else(|| owens_t_znorm1(h));

            0.5 * (normh + normah) - normh * normah - owens_t_dispatch(fabs_ah, fabs_a.recip(), h, Some(znorm1_abs_h))
        } // else [if( h <= 0.67 )]
    }; // else [if(fabs_a <= 1)]

    // exploit that T(h,-a) == -T(h,a)
    if a.is_sign_negative() { -val } else { val }
}

// This routine dispatches the call to one of six subroutines, depending on the values
// of h and a.
// preconditions: h >= 0, 0<=a<=1, ah=a*h, znorm1_ah = owens_t_znorm1(ah) if present
fn owens_t_dispatch(h: f64, a: f64, ah: f64, znorm1_ah: Option<f64>) -> f64 {
    // Simple main case for 64-bit precision or less, this is as per the Patefield-Tandy paper:
    //
    // Handle some special cases first, these are from
    // page 1077 of Owen's original paper:
    //
    if h == 0.0 {
        return atan(a) * ONE_DIV_TWO_PI;
    }
    if a == 0.0 {
        return 0.0;
    }
    if a == 1.0 {
        return owens_t_znorm2(-h) * owens_t_znorm2(h) * 0.5;
    }
    if a >= f64::MAX {
        return owens_t_znorm2(h.abs());
    }

    // Decide which of the 18 regions the pair (h,a) falls into
    let icode = owens_t_compute_code(h, a);
    // Select appropriate algorithm and order to use for the approximation
    let (method, m) = owens_t_get_method_and_order(icode);

    match method {
        1 => owens_t_T1(h, a, m),
        2 => owens_t_T2(h, a, m, ah, znorm1_ah.unwrap_or_else(|| owens_t_znorm1(ah))),
        3 => owens_t_T3(h, a, ah, znorm1_ah.unwrap_or_else(|| owens_t_znorm1(ah))),
        4 => owens_t_T4(h, a, m),
        5 => owens_t_T5(h, a),
        6 => owens_t_T6(h, a),
        _ => {
            panic!("Owen's T algorithm selection failed: icode = {icode}, h = {h}, a = {a}")
        }
    }
}

#[inline(always)]
fn owens_t_znorm1(x: f64) -> f64 {
    0.5 * erf(x * FRAC_1_SQRT_2)
}

#[inline(always)]
fn owens_t_znorm2(x: f64) -> f64 {
    0.5 * erfc(x * FRAC_1_SQRT_2)
}

// compute the value of Owen's T function with method T1 from the reference paper
pub(crate) fn owens_t_T1(h: f64, a: f64, m: u16) -> f64 {
    let hs = -0.5 * (h * h);
    let dhs = exp(hs);
    let a2 = a * a;

    let mut j: u16 = 1;
    let mut jj = 1.0;
    let mut aj = a * ONE_DIV_TWO_PI;
    let mut dj = expm1(hs);
    let mut gj = hs * dhs;

    let mut val = atan(a) * ONE_DIV_TWO_PI;

    loop {
        val += dj * aj / jj;

        if m <= j {
            return val;
        }

        j += 1;
        jj += 2.0;
        aj *= a2;
        dj = gj - dj;
        gj *= hs / (j as f64);
    } // while( true )
}

// compute the value of Owen's T function with method T2 from the reference paper
pub(crate) fn owens_t_T2(h: f64, a: f64, m: u16, ah: f64, znorm1_ah: f64) -> f64 {
    let maxii = m + m + 1;
    let hs = h * h;
    let a2 = -a * a;
    let y = hs.recip();

    let mut ii = 1;
    let mut val = 0.0;
    let mut vi = a * exp(-0.5 * (ah * ah)) * ONE_DIV_ROOT_TWO_PI;
    let mut z = znorm1_ah / h;

    loop {
        val += z;
        if maxii <= ii {
            val *= exp(-0.5 * hs) * ONE_DIV_ROOT_TWO_PI;
            return val;
        } // if( maxii <= ii )
        z = y * (vi - (ii as f64) * z);
        vi *= a2;
        ii += 2;
    } // while( true )
}

// compute the value of Owen's T function with method T3 from the reference paper
// template<class RealType, class Policy>
// inline RealType owens_t_T3_imp(const RealType h, const RealType a, const RealType ah, const std::integral_constant<int, 64>&, const Policy& pol)
pub(crate) fn owens_t_T3(h: f64, a: f64, ah: f64, znorm1_ah: f64) -> f64 {
    const M: usize = 30;

    const C2: [f64; M + 1] = [
        0.99999999999999999999999729978162447266851932041876728736094298092917625009873,
        -0.99999999999999999999467056379678391810626533251885323416799874878563998732905968,
        0.99999999999999999824849349313270659391127814689133077036298754586814091034842536,
        -0.9999999999999997703859616213643405880166422891953033591551179153879839440241685,
        0.99999999999998394883415238173334565554173013941245103172035286759201504179038147,
        -0.9999999999993063616095509371081203145247992197457263066869044528823599399470977,
        0.9999999999797336340409464429599229870590160411238245275855903767652432017766116267,
        -0.999999999574958412069046680119051639753412378037565521359444170241346845522403274,
        0.9999999933226234193375324943920160947158239076786103108097456617750134812033362048,
        -0.9999999188923242461073033481053037468263536806742737922476636768006622772762168467,
        0.9999992195143483674402853783549420883055129680082932629160081128947764415749728967,
        -0.999993935137206712830997921913316971472227199741857386575097250553105958772041501,
        0.99996135597690552745362392866517133091672395614263398912807169603795088421057688716,
        -0.99979556366513946026406788969630293820987757758641211293079784585126692672425362469,
        0.999092789629617100153486251423850590051366661947344315423226082520411961968929483,
        -0.996593837411918202119308620432614600338157335862888580671450938858935084316004769854,
        0.98910017138386127038463510314625339359073956513420458166238478926511821146316469589567,
        -0.970078558040693314521331982203762771512160168582494513347846407314584943870399016019,
        0.92911438683263187495758525500033707204091967947532160289872782771388170647150321633673,
        -0.8542058695956156057286980736842905011429254735181323743367879525470479126968822863,
        0.73796526033030091233118357742803709382964420335559408722681794195743240930748630755,
        -0.58523469882837394570128599003785154144164680587615878645171632791404210655891158,
        0.415997776145676306165661663581868460503874205343014196580122174949645271353372263,
        -0.2588210875241943574388730510317252236407805082485246378222935376279663808416534365,
        0.1375535825163892648504646951500265585055789019410617565727090346559210218472356689,
        -0.0607952766325955730493900985022020434830339794955745989150270485056436844239206648,
        0.0216337683299871528059836483840390514275488679530797294557060229266785853764115,
        -0.00593405693455186729876995814181203900550014220428843483927218267309209471516256,
        0.0011743414818332946510474576182739210553333860106811865963485870668929503649964142,
        -1.489155613350368934073453260689881330166342484405529981510694514036264969925132e-4,
        9.072354320794357587710929507988814669454281514268844884841547607134260303118208e-6,
    ];

    let a2 = a * a;
    let hs = h * h;
    let y = hs.recip();

    let mut ii = 1.0;
    let mut vi = a * exp(-0.5 * (ah * ah)) * ONE_DIV_ROOT_TWO_PI;
    let mut zi = znorm1_ah / h;
    let mut val = 0.0;

    let mut i = 0;
    loop {
        debug_assert!(i < 31);
        val += zi * C2[i];
        if M <= i {
            val *= exp(-0.5 * hs) * ONE_DIV_ROOT_TWO_PI;
            return val;
        }
        zi = y * (ii * zi - vi);
        vi *= a2;
        ii += 2.0;
        i += 1;
    } // while( true )
} // RealType owens_t_T3(const RealType h, const RealType a, const RealType ah)

// compute the value of Owen's T function with method T4 from the reference paper
//template<typename RealType>
//inline RealType owens_t_T4(const RealType h, const RealType a, const unsigned short m)
pub(crate) fn owens_t_T4(h: f64, a: f64, m: u16) -> f64 {
    let maxii = m + m + 1;
    let hs = h * h;
    let a2 = -a * a;

    let mut ii: u16 = 1;
    let mut ai = a * exp(-0.5 * hs * (1.0 - a2)) * ONE_DIV_TWO_PI;
    let mut yi = 1.0;
    let mut val = 0.0;

    loop {
        val += ai * yi;
        if maxii <= ii {
            return val;
        }
        ii += 2;
        yi = (1.0 - hs * yi) / (ii as f64);
        ai *= a2;
    } // while( true )
} // RealType owens_t_T4(const RealType h, const RealType a, const unsigned short m)

// compute the value of Owen's T function with method T5 from the reference paper
// template<typename RealType>
//inline RealType owens_t_T5_imp(const RealType h, const RealType a, const std::integral_constant<int, 64>&)
pub(crate) fn owens_t_T5(h: f64, a: f64) -> f64 {
    /*
      NOTICE:
      - The pts[] array contains the squares (!) of the abscissas, i.e. the roots of the Legendre
      polynomial P_n(x), instead of the plain roots as required in Gauss-Legendre
      quadrature, because T5(h,a,m) contains only x^2 terms.
      - The wts[] array contains the weights for Gauss-Legendre quadrature scaled with a factor
      of 1/(2*pi) according to T5(h,a,m).
    */

    const M: usize = 19;
    const PTS: [f64; M] = [
        0.0016634282895983227941,
        0.014904509242697054183,
        0.04103478879005817919,
        0.079359853513391511008,
        0.1288612130237615133,
        0.18822336642448518856,
        0.25586876186122962384,
        0.32999972011807857222,
        0.40864620815774761438,
        0.48971819306044782365,
        0.57106118513245543894,
        0.6505134942981533829,
        0.72596367859928091618,
        0.79540665919549865924,
        0.85699701386308739244,
        0.90909804422384697594,
        0.95032536436570154409,
        0.97958418733152273717,
        0.99610366384229088321,
    ];
    const WTS: [f64; M] = [
        0.012975111395684900835,
        0.012888764187499150078,
        0.012716644398857307844,
        0.012459897461364705691,
        0.012120231988292330388,
        0.011699908404856841158,
        0.011201723906897224448,
        0.010628993848522759853,
        0.0099855296835573320047,
        0.0092756136096132857933,
        0.0085039700881139589055,
        0.0076757344408814561254,
        0.0067964187616556459109,
        0.005871875456524750363,
        0.0049082589542498110071,
        0.0039119870792519721409,
        0.0028897090921170700834,
        0.0018483371329504443947,
        0.00079623320100438873578,
    ];

    let a2 = a * a;
    let hs = -0.5 * h * h;

    let mut val = 0.0;
    for i in 0..M {
        let r = 1.0 + a2 * PTS[i];
        val += WTS[i] * exp(hs * r) / r;
    } // for(unsigned short i = 0; i < m; ++i)

    val * a
} // RealType owens_t_T5(const RealType h, const RealType a)

// compute the value of Owen's T function with method T6 from the reference paper
// template<typename RealType, class Policy>
// inline RealType owens_t_T6(const RealType h, const RealType a, const Policy& pol)
pub(crate) fn owens_t_T6(h: f64, a: f64) -> f64 {
    let normh = owens_t_znorm2(h);
    let y = 1.0 - a;
    let r = atan2(y, 1.0 + a);

    let mut val = normh * (1.0 - normh) * 0.5;

    if r != 0.0 {
        val -= r * exp(-0.5 * y * (h * h) / r) * ONE_DIV_TWO_PI;
    }

    val
} // RealType owens_t_T6(const RealType h, const RealType a, const unsigned short m)

// Auxiliary function, it computes an array key that is used to determine
// the specific computation method for Owen's T and the order thereof
// used in owens_t_dispatch.
// template<typename RealType>
// inline unsigned short owens_t_compute_code(const RealType h, const RealType a)
pub(crate) fn owens_t_compute_code(h: f64, a: f64) -> u16 {
    const HRANGE: [f64; 14] = [
        0.02, 0.06, 0.09, 0.125, 0.26, 0.4, 0.6, 1.6, 1.7, 2.33, 2.4, 3.36, 3.4, 4.8,
    ];

    const ARANGE: [f64; 7] = [0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999];
    /*
    original select array from paper:
    1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9
    1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9
    2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10
    2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10
    2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11
    2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12
    2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12
    2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12
    */
    // subtract one because the array is written in FORTRAN in mind - in C arrays start @ zero
    const SELECT: [u16; 120] = [
        0, 0, 1, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 8, 0, 1, 1, 2, 2, 4, 4, 13, 13, 14,
        14, 15, 15, 15, 8, 1, 1, 2, 2, 2, 4, 4, 14, 14, 14, 14, 15, 15, 15, 9, 1, 1, 2, 4, 4, 4, 4,
        6, 6, 15, 15, 15, 15, 15, 9, 1, 2, 2, 4, 4, 5, 5, 7, 7, 16, 16, 16, 11, 11, 10, 1, 2, 4, 4,
        4, 5, 5, 7, 7, 16, 16, 16, 11, 11, 11, 1, 2, 3, 3, 5, 5, 7, 7, 16, 16, 16, 16, 16, 11, 11,
        1, 2, 3, 3, 5, 5, 17, 17, 17, 17, 16, 16, 16, 11, 11,
    ];

    let ihint = HRANGE.iter().position(|val| h <= *val).unwrap_or(14);
    let iaint = ARANGE.iter().position(|val| a <= *val).unwrap_or(7);

    // interpret select array as 8x15 matrix
    SELECT[iaint * 15 + ihint]
}

//template<typename RealType>
//  inline unsigned short owens_t_get_order_imp(const unsigned short icode, RealType, const std::integral_constant<int, 64>&)
pub(crate) fn owens_t_get_method_and_order(icode: u16) -> (u16, u16) {
    const METHOD: [u16; 18] = [
        //
        1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6,
    ];
    const ORD: [u16; 18] = [
        3, 4, 5, 6, 8, 11, 13, 19, 10, 20, 30, 0, 7, 10, 11, 23, 0, 0,
    ];

    debug_assert!(icode < 18);
    (METHOD[icode as usize], ORD[icode as usize])
} // unsigned short owens_t_get_order(const unsigned short icode, RealType, std::integral_constant<int, 64> const&)

#[cfg(test)]
mod tests {
    use super::*;
    use assert_within::assert_within;

    // This T7 routine is described in Patefield-Tandy for validating the other routines
    // C2 array is chebyshev coefficients from Patefield-Tandy fortran sources.
    // TODO: Should this actually be replaced with the C2 from T3 routine? I'm not sure.
    const C2: [f64; 21] = [
        0.99999999999999987510,
        -0.99999999999988796462,
        0.99999999998290743652,
        -0.99999999896282500134,
        0.99999996660459362918,
        -0.99999933986272476760,
        0.99999125611136965852,
        -0.99991777624463387686,
        0.99942835555870132569,
        -0.99697311720723000295,
        0.98751448037275303682,
        -0.95915857980572882813,
        0.89246305511006708555,
        -0.76893425990463999675,
        0.58893528468484693250,
        -0.38380345160440256652,
        0.20317601701045299653,
        -0.82813631607004984866e-01,
        0.24167984735759576523e-01,
        -0.44676566663971825242e-02,
        0.39141169402373836468e-03,
    ];
    fn owens_t_T7(h: f64, a: f64) -> f64 {
        let h2 = h * h;
        let a2 = a * a;

        let mut u = 0.5 * FRAC_1_PI * a * exp(-0.5 * h2 * (1.0 + a2));
        let mut v = C2[0];

        let mut result = 0.0;
        let mut k: u16 = 0;
        while (u * v).abs() > 1e-20 {
            result += u * v;

            u *= a2;
            k += 1;

            v *= h2;
            if k < C2.len() as u16 {
                v += C2[k as usize];
            }
            v /= (2 * k + 1) as f64;
        }
        result
    }

    #[test]
    fn compare_with_T7() {
        for h in [
            0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0
        ] {
            for a in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95] {
                assert_within!(+1e-16, owens_t(h, a), owens_t_T7(h, a), "h = {h}, a = {a}");
            }
        }
    }

    #[test]
    fn compare_with_T7_extreme_a() {
        for h in [
            0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5
        ] {
            for a in [0.005, 0.01, 0.99, 0.995] {
                assert_within!(+1e-16, owens_t(h, a), owens_t_T7(h, a), "h = {h}, a = {a}");
            }
        }

        for h in [
            1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0
        ] {
            // TODO: Figure out why these have less agreement, maybe using T7 wrong...
            // or maybe T7 just isn't accurate here, we appear to agree with wolfram alpha on these points
            {
                let a = 0.005;
                assert_within!(+1e-12, owens_t(h, a), owens_t_T7(h, a), "h = {h}, a = {a}");
            }

            {
                let a = 0.01;
                assert_within!(+1e-10, owens_t(h, a), owens_t_T7(h, a), "h = {h}, a = {a}");
            }

            for a in [0.99, 0.995, 0.999] {
                assert_within!(+1e-16, owens_t(h, a), owens_t_T7(h, a), "h = {h}, a = {a}");
            }

        }
    }

    #[test]
    fn spot_check_owen_t() {
        let eps = 1e-16;
        // These values from wolframalpha OwenT function
        assert_within!(+eps, owens_t(4.0, 1.0), 0.00001583511938278);
        assert_within!(+eps, owens_t(3.0, 1.0), 0.00067403790346715);
        assert_within!(+eps, owens_t(2.0, 1.0), 0.01111628172225982);
        assert_within!(+eps, owens_t(1.0, 1.0), 0.06674188216570097);
        assert_within!(+eps, owens_t(0.0, 1.0), 0.125);
        assert_within!(+eps, owens_t(-1.0, 1.0), 0.06674188216570097);
        assert_within!(+eps, owens_t(-2.0, 1.0), 0.01111628172225982);
        assert_within!(+eps, owens_t(-3.0, 1.0), 0.00067403790346715);
        assert_within!(+eps, owens_t(-4.0, 1.0), 0.00001583511938278);
        assert_within!(+eps, owens_t(0.25, 0.25), 0.03776579074789887);
        assert_within!(+eps, owens_t(0.35, 0.25), 0.03662713104215306);
        assert_within!(+eps, owens_t(0.45, 0.25), 0.03516215929701228);
        assert_within!(+eps, owens_t(0.55, 0.25), 0.03341309193133329);
        assert_within!(+eps, owens_t(0.65, 0.25), 0.03142870674717572);
        assert_within!(+eps, owens_t(0.75, 0.25), 0.02926208986224979);
        assert_within!(+eps, owens_t(0.85, 0.25), 0.02696829336798952);
        assert_within!(+eps, owens_t(0.95, 0.25), 0.02460204670946997);

        // These are values from Wolfram alpha at locations where at some point
        // disagreed with T7
        assert_within!(+eps, owens_t(1.0, 0.05), 0.00482059414583312869);
        assert_within!(+eps, owens_t(1.0, 0.01), 0.00096527526326129856);
        assert_within!(+eps, owens_t(1.0, 0.005), 0.00048265572997626906);
        assert_within!(+eps, owens_t(2.0, 0.05), 0.00107427827101219616);
        assert_within!(+eps, owens_t(2.0, 0.01), 0.00021537125589291752);
        assert_within!(+eps, owens_t(2.0, 0.005), 0.00010769370416663914);

        let eps = 1e-14;
        assert_within!(~eps, owens_t(5.0, 0.05), 2.93255054168372620347e-8);
        assert_within!(~eps, owens_t(5.0, 0.01), 5.92848480305363926081e-9);
        assert_within!(~eps, owens_t(5.0, 0.005), 2.96524277424805823103e-9);
        assert_within!(~eps, owens_t(10.0, 0.05), 1.4720412827780475392e-24);
        assert_within!(~eps, owens_t(10.0, 0.01), 3.06449020934741106603e-25);
        assert_within!(~eps, owens_t(10.0, 0.005), 1.5341982995816032882e-25);
    }

    #[test]
    fn spot_check_owen_t_large_a() {
        let eps = 1e-16;
        assert_within!(+eps, owens_t(1.0, 10.0), 0.079327626965728525707);
        assert_within!(+eps, owens_t(1.1, 10.0), 0.067833030473191337587);
        assert_within!(+eps, owens_t(1.2, 10.0), 0.057534835110854134011);
        assert_within!(+eps, owens_t(1.5, 10.0), 0.033403600634429033002);
        assert_within!(+eps, owens_t(2.0, 10.0), 0.011375065974089603600);
        let eps = 1e-14;
        assert_within!(~eps, owens_t(5.0, 10.0), 1.433257859395969558369e-7);
        assert_within!(~eps, owens_t(10.0, 10.0), 3.809926512080263032987e-24);
    }
}
