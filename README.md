owens-t
=======

[![Crates.io](https://img.shields.io/crates/v/owens-t?style=flat-square)](https://crates.io/crates/owens-t)
[![Crates.io](https://img.shields.io/crates/d/owens-t?style=flat-square)](https://crates.io/crates/owens-t)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue?style=flat-square)](LICENSE-APACHE)
[![License](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](LICENSE-MIT)

[API Docs](https://docs.rs/owens-t/latest)

This crate provides an implementation of:

 * [Owen's T function](https://en.wikipedia.org/wiki/Owen%27s_T_function)
 * The [bivariate normal CDF](https://en.wikipedia.org/wiki/Multivariate_normal_distribution) (using Owen's T)
 
The Owen's T implementaiton is accurate to about 16 decimals.
The bivariate normal CDF is accurate to about 9-10 decimals in my testing.

On my x86-64 hardware, the cost of Owen's T is typically about 3x the cost
of `libm::erfc`. The cost of bivariate normal cdf computation is two calls to `erfc` and two calls to `owens_t`.

Owen's T can also be used to compute the CDF of the [skew-normal distribution](https://en.wikipedia.org/wiki/Skew_normal_distribution).

Quick-start
-----------

To evaluate `T(h,a) = 1/2π · ∫₀^a exp(½h²·(1+x²)) / (1+x²) dx`:

```rust
owens_t::owens_t(h,a)
```

To evaluate `Pr[X > a, Y > b]` where `X` and `Y` are mean 0 variance 1 Gaussians of correlation coefficient `ρ`:

```rust
owens_t::biv_norm(a, b, rho)
```

If you will use the same arguments with `biv_norm` repeatedly, you can save very modestly by using a slightly different API that caches some of the shared computation.

```rust
let a_vec = vec![0.1, 0.2, 0.3].into_iter().map(BivNormArg::from).collect::<Vec<_>>();
let b_vec = vec![0.4, 0.5, 0.6].into_iter().map(BivNormArg::from).collect::<Vec<_>>();
let rho = BivNormRho::from(0.7);

let mut answers = vec![];
for a in a_vec {
    for b in b_vec {
        answers.push(owens_t::biv_norm_inner(a, b, rho));
    }
}
```

Algorithm Notes
---------------

This is a partial port of the [C++ implementation](https://live.boost.org/doc/libs/1_81_0/boost/math/special_functions/owens_t.hpp) in [`boost::math`](https://live.boost.org/doc/libs/1_81_0/libs/math/doc/html/math_toolkit/owens_t.html).

* Owen's T was introduced [in 1954 by Owen](http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.aoms/1177728074) and used in the resolution of many Gaussian-type integrals and statistical applications. Owen also showed that it can be used to compute the bivariate normal CDF, which has scientific applications.
* The Patefield-Tandy algorithm and Fortran code was [published in 2000](https://www.jstatsoft.org/article/view/v005i05), which is close to state of the art, and achieved significantly higher precision for all inputs than all previous implementations.

The `boost::math` implementation is generic over an integer type, and uses SFINAE to change algorithms based on the precision of the integer.

> The Patefield-Tandy algorithm provides six methods of evaluation (`T1` to `T6`); the best method is selected according to the values of `a` and `h`.
> The Patefield-Tandy algorithm is accurate to approximately 20 decimal places, so for types with greater precision we use:
>
>  A modified version of T1 which folds the calculation of atan(h) into the T1 series (to avoid subtracting two values similar in magnitude), and then accelerates the resulting alternating series using method 1 from H. Cohen, F. Rodriguez Villegas, D. Zagier, "Convergence acceleration of alternating series", Bonn, (1991). The result is valid everywhere, but doesn't always converge, or may become too divergent in the first few terms to sum accurately. This is used for `ah < 1`.
> * A modified version of T2 which is accelerated in the same manner as T1. This is used for h > 1.
> * A version of T4 only when both T1 and T2 have failed to produce an accurate answer.
> * Fallback to the Patefiled Tandy algorithm when all the above methods fail: this happens not at all for our test data at 100 decimal digits precision. However, there is a difficult area when a is very close to 1 and the precision increases which may cause this to happen in very exceptional circumstances.
>
> Using the above algorithm and a 100-decimal digit type, results accurate to 80 decimal places were obtained in the difficult area where a is close to 1, and greater than 95 decimal places elsewhere.

I had no particular need to get 80 decimals of accuracy, and my goals were to get about 10 decimals of accuracy very rapidly. Also I'm a bit spooked by these series that may not converge everywhere -- the boost implementation attempts to throw exceptions when they don't converge, and I'd rather not have that either.

So, I did not implement their modifications to T1, T2, T4. I only implemented classic Patefield-Tandy, and I only did it for `f64`. Then I tested carefully, and measured performance, which was satisfactory for my needs.

Additionally, I found that `biv_norm` computes values of `Φ(x)` that are often recomputed by Patefield-Tandy, so I created an `owens_t_inner` that can optionally forward those values. This modestly accelerates `biv_norm`.

Future Directions
-----------------

* It would be interesting to have a trait-based API and do it for `f32` also, especially if it can be faster. Or to make the implementation generic over e.g. `num_traits::FloatCore`. Ideally we'd be able to do precision-based specialization similar to `boost::math`.
* It would be interesting if SIMD instructions can be used to evaluate Owen's T at multiple positions in parallel. (It seems difficult to use SIMD productively for just a single point.) The main stumbling block is, if those positions are not classified to use the same algorithm by the Patefield-Tandy selector, then most likely it won't work. However, if it works often when the inputs are close, it may still be worthwhile. As long as the points are classified to use the same series, it will work, even if the algorithm calls for a larger order for some points than for others. This is one argument in favor of using the `boost::math` "accelerated" series at least some of the time, since `boost::math` is using fewer series on a wider range of inputs.

Licensing and Distribution
--------------------------

MIT or Apache 2 at your option.
