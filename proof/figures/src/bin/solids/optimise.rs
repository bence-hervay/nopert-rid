//! Nelder–Mead minimisation with random restarts, for finding passages to
//! draw.

/// A small deterministic generator (xorshift64*), uniform in `[0, 1)`.
pub struct Random(u64);

impl Random {
    pub fn new(seed: u64) -> Self {
        Random(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }

    pub fn next(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A local minimum of `f` near `x0`, with its value.
pub fn nelder_mead(f: &impl Fn(&[f64]) -> f64, x0: &[f64], step: f64, iterations: usize) -> (Vec<f64>, f64) {
    let n = x0.len();
    let mut simplex: Vec<(Vec<f64>, f64)> = (0..=n)
        .map(|i| {
            let mut x = x0.to_vec();
            if i > 0 {
                x[i - 1] += step;
            }
            let v = f(&x);
            (x, v)
        })
        .collect();
    let combine = |a: &[f64], b: &[f64], t: f64| a.iter().zip(b).map(|(x, y)| x + t * (y - x)).collect::<Vec<f64>>();
    for _ in 0..iterations {
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        let mut centroid = vec![0.0; n];
        for (x, _) in &simplex[..n] {
            for k in 0..n {
                centroid[k] += x[k] / n as f64;
            }
        }
        let worst = simplex[n].clone();
        let reflected = combine(&centroid, &worst.0, -1.0);
        let fr = f(&reflected);
        if fr < simplex[0].1 {
            let expanded = combine(&centroid, &worst.0, -2.0);
            let fe = f(&expanded);
            simplex[n] = if fe < fr { (expanded, fe) } else { (reflected, fr) };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (reflected, fr);
        } else {
            let contracted = combine(&centroid, &worst.0, 0.5);
            let fc = f(&contracted);
            if fc < worst.1 {
                simplex[n] = (contracted, fc);
            } else {
                let best = simplex[0].0.clone();
                for entry in simplex.iter_mut().skip(1) {
                    let x = combine(&best, &entry.0, 0.5);
                    *entry = (x.clone(), f(&x));
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
    simplex.swap_remove(0)
}

/// The best of `starts` local minimisations from points drawn by `start`,
/// each polished by a second, finer run.
pub fn minimise(
    f: &impl Fn(&[f64]) -> f64,
    starts: usize,
    seed: u64,
    start: impl Fn(&mut Random) -> Vec<f64>,
) -> (Vec<f64>, f64) {
    let mut random = Random::new(seed);
    let mut best: Option<(Vec<f64>, f64)> = None;
    for _ in 0..starts {
        let x0 = start(&mut random);
        let (x, _) = nelder_mead(f, &x0, 0.3, 400);
        let (x, _) = nelder_mead(f, &x, 0.02, 400);
        let (x, v) = nelder_mead(f, &x, 0.001, 600);
        if best.as_ref().map_or(true, |b| v < b.1) {
            best = Some((x, v));
        }
    }
    best.expect("at least one start")
}
