use std::collections::HashSet;
use std::{iter::once, usize};

use plotters::style::full_palette;
use rand::Rng;

use crate::gray_flips::GrayFlips;
use crate::lattice::Lattice;
use crate::plots::Series;

pub fn run_ising_simulation<Sample, Accumulator>(
    sampler: impl Iterator<Item = Sample>,
    mut trackers: Vec<&mut dyn ProgressTracker<Sample>>,
    aggregator: impl SampleAggregator<Sample, Accumulator>,
) -> Vec<Series> {
    let acc = (0..)
        .zip(sampler)
        .inspect(|x| {
            trackers
                .iter_mut()
                .for_each(|tracker| tracker.track(x.0, &x.1));
        })
        .fold(aggregator.init(), |mut acc, sample| {
            aggregator.handle(&mut acc, sample.1);
            acc
        });

    return aggregator.series(acc);
}

// Helpers

#[inline(always)]
fn energy_from_idx(idx: usize, lattice: Lattice) -> i32 {
    idx as i32 - lattice.edges as i32
}

#[inline(always)]
fn boltzmann_weight(idx: usize, lattice: Lattice, beta: f64) -> f64 {
    f64::exp(-beta * energy_from_idx(idx, lattice) as f64)
}

#[inline(always)]
fn magnet_from_idx(idx: usize, lattice: Lattice) -> i32 {
    idx as i32 - (lattice.vertices as i32)
}

// Samplers

pub type IsingSample = (usize, usize);

/// This will sample all lattice configurations states exhaustively
/// (which means 2^(n^2) possibilities -- really slow for n > 4)
pub fn gray_sampler(lattice: Lattice) -> impl Iterator<Item = IsingSample> {
    let (init_energy, init_magn) = (0_usize, 2 * lattice.vertices);
    once((init_energy, init_magn)).chain(GrayFlips::new(lattice.vertices).scan(
        (vec![true; lattice.vertices], init_energy, init_magn),
        move |(spins, energy, magn), idx| {
            let diff = if spins[idx] { -2 } else { 2 };
            spins[idx] = !spins[idx];
            *energy = (*energy as i32 - lattice.neighbouring(idx, spins) * diff) as usize;
            *magn = (*magn as i32 + diff) as usize;
            return Some((*energy, *magn));
        },
    ))
}

/// Naive probabilistic sampler based on a Markov chain
/// which inverts one spin at a time
pub fn mcmc_local_sampler(
    iter: usize,
    lattice: Lattice,
    beta: f64,
) -> impl Iterator<Item = IsingSample> {
    let (init_energy, init_magn) = (0_usize, 2 * lattice.vertices);
    let mut rng1 = rand::thread_rng();
    let mut rng2 = rand::thread_rng();
    (0..iter)
        .map(move |_| rng1.gen_range(0..lattice.vertices))
        .scan(
            (vec![true; lattice.vertices], init_energy, init_magn),
            move |(spins, energy, magn), idx| {
                let diff = if spins[idx] { -2 } else { 2 };
                let ngb_spins_diff = lattice.neighbouring(idx, spins) * diff;
                if rng2.gen::<f64>() < f64::exp(beta * ngb_spins_diff as f64) {
                    spins[idx] = !spins[idx];
                    *energy = (*energy as i32 - ngb_spins_diff) as usize;
                    *magn = (*magn as i32 + diff) as usize;
                }
                return Some((*energy, *magn));
            },
        )
}

/// Cluster algorithm based on Wolff algorithm as described
/// in <https://arxiv.org/abs/cond-mat/0311623>
pub fn mcmc_cluster_sampler<'a>(
    iter: usize,
    lattice: Lattice,
    beta: f64,
    final_inspect: &'a dyn Fn(&Vec<bool>),
) -> impl Iterator<Item = IsingSample> + 'a {
    let p = 1. - f64::exp(-2. * beta);
    let (init_energy, init_magn) = (0_usize, 2 * lattice.vertices);
    let mut rng1 = rand::thread_rng();
    let mut rng2 = rand::thread_rng();
    (0..iter)
        .map(move |it| (it, rng1.gen_range(0..lattice.vertices)))
        .scan(
            (vec![true; lattice.vertices], init_energy, init_magn),
            move |(spins, energy, magn), (it, idx)| {
                let mut cluster = HashSet::from([idx]);
                let mut pocket = vec![idx];
                let cluster_spin = spins[idx];
                while let Some(k) = pocket.pop() {
                    for l in lattice.neighbours(k) {
                        if !cluster.contains(&l)
                            && spins[l] == cluster_spin
                            && rng2.gen::<f64>() < p
                        {
                            pocket.push(l);
                            cluster.insert(l);
                        }
                    }
                }
                let diff = if cluster_spin { -2 } else { 2 };
                for k in cluster {
                    let ngb_spins_diff = lattice.neighbouring(k, spins) * diff;
                    spins[k] = !spins[k];
                    *energy = (*energy as i32 - ngb_spins_diff) as usize;
                    *magn = (*magn as i32 + diff) as usize;
                }
                if it + 1 == iter {
                    final_inspect(spins);
                }
                return Some((*energy, *magn));
            },
        )
}

// Progress trackers

pub trait ProgressTracker<Sample> {
    fn track(&mut self, it: usize, sample: &Sample);
}

pub struct AvgTracker {
    period: usize,
    avg_energy: f64,
    avg_magn: f64,
}

impl AvgTracker {
    pub fn new(period: usize) -> AvgTracker {
        AvgTracker {
            period: period,
            avg_energy: 0.,
            avg_magn: 0.,
        }
    }
}

#[inline(always)]
fn online_avg(it: usize, avg: f64, current: usize) -> f64 {
    (avg * (it as f64) + current as f64) / (it as f64 + 1.)
}

impl ProgressTracker<IsingSample> for AvgTracker {
    fn track(&mut self, it: usize, sample: &IsingSample) {
        self.avg_energy = online_avg(it, self.avg_energy, sample.0);
        self.avg_magn = online_avg(it, self.avg_magn, sample.1);
        if (it + 1) % self.period == 0 {
            println!(
                "{}: {:?} {:?}",
                it,
                sample,
                (self.avg_energy, self.avg_magn)
            );
        }
    }
}

// Aggregators

pub trait SampleAggregator<Sample, Accumulator> {
    fn init(&self) -> Accumulator;
    fn handle(&self, acc: &mut Accumulator, sample: Sample);
    fn series(&self, acc: Accumulator) -> Vec<Series>;
}

/// Aggregator compatible with exhaustive sampling strategy,
/// computes histograms and Boltzmann distributions
pub struct ExhaustiveAggregator {
    lattice: Lattice,
    beta: f64,
}

impl ExhaustiveAggregator {
    pub fn new(lattice: Lattice, beta: f64) -> ExhaustiveAggregator {
        ExhaustiveAggregator {
            lattice: lattice,
            beta: beta,
        }
    }
}

pub type ExhaustiveAccumulator = (f64, Vec<i32>, Vec<f64>, Vec<i32>, Vec<f64>);

impl SampleAggregator<IsingSample, ExhaustiveAccumulator> for ExhaustiveAggregator {
    fn init(&self) -> ExhaustiveAccumulator {
        (
            0.,
            vec![0; 2 * self.lattice.edges + 1],
            vec![0.; 2 * self.lattice.edges + 1],
            vec![0; 2 * self.lattice.vertices + 1],
            vec![0.; 2 * self.lattice.vertices + 1],
        )
    }
    fn handle(&self, acc: &mut ExhaustiveAccumulator, sample: IsingSample) {
        let weight = boltzmann_weight(sample.0, self.lattice, self.beta);
        acc.0 += weight;
        acc.1[sample.0] += 1;
        acc.2[sample.0] += weight;
        acc.3[sample.1] += 1;
        acc.4[sample.1] += weight;
    }
    fn series(&self, acc: ExhaustiveAccumulator) -> Vec<Series> {
        return [
            (
                "Energy histogram",
                full_palette::RED,
                acc.1
                    .iter()
                    .zip(0..)
                    .map(|(w, idx)| (energy_from_idx(idx, self.lattice), *w))
                    .filter(|(_, w)| *w > 0)
                    .collect(),
            ),
            (
                "Energy probabilities",
                full_palette::RED,
                acc.2
                    .iter()
                    .zip(0..)
                    .map(|(w, idx)| (energy_from_idx(idx, self.lattice), *w / acc.0))
                    .filter(|(_, w)| *w > f64::EPSILON)
                    .map(|(x, w)| (x, (w / 0.01) as i32))
                    .collect(),
            ),
            (
                "Magnetization histogram",
                full_palette::TEAL,
                acc.3
                    .iter()
                    .zip(0..)
                    .map(|(w, idx)| (magnet_from_idx(idx, self.lattice), *w))
                    .filter(|(_, w)| *w > 0)
                    .collect(),
            ),
            (
                "Magnetization probabilities",
                full_palette::TEAL,
                acc.4
                    .iter()
                    .zip(0..)
                    .map(|(w, idx)| (magnet_from_idx(idx, self.lattice), *w / acc.0))
                    .filter(|(_, w)| *w > f64::EPSILON)
                    .map(|(x, w)| (x, (w / 0.01) as i32))
                    .collect(),
            ),
        ]
        .map(|(title, color, data)| Series {
            name: format!(
                "[d={}, β={}, n={}, {}] {}",
                self.lattice.d, self.beta, self.lattice.n, self.lattice.boundary_conditions, title
            ),
            color: color,
            data: data,
        })
        .to_vec();
    }
}

/// Aggregator compatible with MCMC sampling strategies
pub struct MCMCAggregator {
    iterations: usize,
    lattice: Lattice,
    beta: f64,
}

impl MCMCAggregator {
    pub fn new(iterations: usize, lattice: Lattice, beta: f64) -> MCMCAggregator {
        MCMCAggregator {
            iterations: iterations,
            lattice: lattice,
            beta: beta,
        }
    }
}

type MCMCAccumulator = (Vec<i32>, Vec<i32>);

fn bucketize(
    data: Vec<i32>,
    n_buckets: usize,
    iterations: usize,
    value_from_idx: &dyn Fn(usize) -> i32,
) -> Vec<(i32, i32)> {
    let bucket_width = usize::div_ceil(data.len(), n_buckets);
    let mut buckets = vec![0; n_buckets];
    for (idx, w) in data.iter().enumerate() {
        buckets[idx / bucket_width] += w;
    }
    return buckets
        .into_iter()
        .zip(0..)
        .map(|(w, idx)| (value_from_idx(idx * bucket_width + bucket_width / 2), w))
        .map(|(x, w)| (x, ((w as f64 / iterations as f64) / 0.01) as i32))
        .collect();
}

impl SampleAggregator<IsingSample, MCMCAccumulator> for MCMCAggregator {
    fn init(&self) -> MCMCAccumulator {
        (
            vec![0; 2 * self.lattice.edges + 1],
            vec![0; 2 * self.lattice.vertices + 1],
        )
    }
    fn handle(&self, acc: &mut MCMCAccumulator, sample: IsingSample) {
        acc.0[sample.0] += 1;
        acc.1[sample.1] += 1;
    }
    fn series(&self, acc: MCMCAccumulator) -> Vec<Series> {
        let bucket_size = 100;
        return [
            (
                "Energy probabilities",
                full_palette::RED,
                bucketize(acc.0, bucket_size, self.iterations, &|idx| {
                    energy_from_idx(idx, self.lattice)
                }),
            ),
            (
                "Magnetization probabilities",
                full_palette::TEAL,
                bucketize(acc.1, bucket_size, self.iterations, &|idx| {
                    magnet_from_idx(idx, self.lattice)
                }),
            ),
        ]
        .map(|(title, color, data)| Series {
            name: format!(
                "[d={}, β={}, n={}] {}",
                self.lattice.d, self.beta, self.lattice.n, title
            ),
            color: color,
            data: data,
        })
        .to_vec();
    }
}
