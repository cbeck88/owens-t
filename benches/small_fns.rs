use criterion::{Criterion, black_box, criterion_group, criterion_main};
use owens_t::owens_t;
use rand::Rng;
use rand_pcg::{Pcg64Mcg, rand_core::SeedableRng};

fn div_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("f64::div", |b| {
        b.iter(|| {
            let arg: f64 = black_box(rng.random::<f64>());
            let arg2: f64 = black_box(rng.random::<f64>());
            arg / arg2
        })
    });
}

fn sqrt_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("f64::sqrt", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            arg.sqrt()
        })
    });
}

fn exp_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("f64::exp", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            arg.exp()
        })
    });
}

fn atan_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::atan", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            libm::atan(arg)
        })
    });
}

fn atan2_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::atan2", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            let arg2: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            libm::atan2(arg, arg2)
        })
    });
}

fn atan2f_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::atan2f", |b| {
        b.iter(|| {
            let arg: f32 = black_box(4.0 * rng.random::<f32>() - 2.0);
            let arg2: f32 = black_box(4.0 * rng.random::<f32>() - 2.0);
            libm::atan2f(arg, arg2)
        })
    });
}

fn erfc_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::erfc", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            libm::erfc(arg)
        })
    });
}

fn erfcf_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::erfcf", |b| {
        b.iter(|| {
            let arg: f32 = black_box(4.0 * rng.random::<f32>() - 2.0);
            libm::erfcf(arg)
        })
    });
}

fn expm1_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::expm1", |b| {
        b.iter(|| {
            let arg: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            libm::expm1(arg)
        })
    });
}

fn expm1f_bench(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("libm::expm1f", |b| {
        b.iter(|| {
            let arg: f32 = black_box(4.0 * rng.random::<f32>() - 2.0);
            libm::expm1f(arg)
        })
    });
}

fn owens_t_bench1(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("owens_t", |b| {
        b.iter(|| {
            let arg1: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            let arg2: f64 = black_box(4.0 * rng.random::<f64>() - 2.0);
            owens_t(arg1, arg2)
        })
    });
}

fn owens_t_bench2(c: &mut Criterion) {
    let mut rng = Pcg64Mcg::seed_from_u64(9);
    c.bench_function("owens_t", |b| {
        b.iter(|| {
            let arg1: f64 = black_box(rng.random::<f64>());
            let arg2: f64 = black_box(rng.random::<f64>());
            let arg3: f64 = black_box(rng.random::<f64>());
            owens_t(arg1 / arg2, arg3)
        })
    });
}

criterion_group!(
    benches,
    div_bench,
    sqrt_bench,
    exp_bench,
    atan_bench,
    atan2_bench,
    atan2f_bench,
    erfc_bench,
    erfcf_bench,
    expm1_bench,
    expm1f_bench,
    owens_t_bench1,
    owens_t_bench2
);
criterion_main!(benches);
