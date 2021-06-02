#![allow(dead_code)]

mod registration;

// Help: https://docs.rs/libR-sys, https://github.com/hadley/r-internals

use dahl_salso::clustering::Clusterings;
use dahl_salso::optimize::{minimize_by_salso, SALSOParameters};
use dahl_salso::{LabelType, LossFunction, PartitionDistributionInformation};
use roxido::*;
use epa::epa::{sample, EpaParameters, SquareMatrixBorrower};
use epa::perm::Permutation;
use rand::SeedableRng;
use rand::Rng;
use rand_pcg::Pcg64Mcg;
use std::convert::TryFrom;

fn sample_epa_engine<T: Rng>(
    n_samples: usize,
    n_items: usize,
    similarity: &[f64],
    mass: f64,
    discount: f64,
    n_cores: usize,
    rng: &mut T,
) -> (Vec<LabelType>, Vec<LabelType>) {
    let n_cores = if n_cores == 0 {
        num_cpus::get()
    } else {
        n_cores
    };
    let n_samples = n_samples.max(1);
    let n_samples_per_core = 1 + (n_samples - 1) / n_cores;
    let chunk_size = n_samples_per_core * n_items;
    let mut samples: Vec<LabelType> = vec![0; n_cores * chunk_size];
    let mut n_clusters: Vec<LabelType> = vec![0; n_cores * n_samples_per_core];

    let _result = crossbeam::scope(|s| {
        let mut stick1 = &mut samples[..];
        let mut stick2 = &mut n_clusters[..];
        let mut plan = Vec::with_capacity(n_cores);
        for _ in 0..n_cores - 1 {
            let (left1, right1) = stick1.split_at_mut(chunk_size);
            let (left2, right2) = stick2.split_at_mut(n_samples_per_core);
            plan.push((left1, left2, rng.gen::<u128>()));
            stick1 = right1;
            stick2 = right2;
        }
        plan.push((stick1, stick2, rng.gen()));
        let sim = SquareMatrixBorrower::from_slice(similarity, n_items);
        plan.into_iter().for_each(|p| {
            s.spawn(move |_| {
                let mut rng = Pcg64Mcg::new(p.2);
                let mut params =
                    EpaParameters::new(sim, Permutation::natural(n_items), mass, discount).unwrap();
                for i in 0..n_samples_per_core {
                    params.shuffle_permutation(&mut rng);
                    let clustering = sample(&params, &mut rng);
                    let zero: LabelType = 0;
                    clustering.relabel_into_slice(zero, &mut p.0[i * n_items..(i + 1) * n_items]);
                    p.1[i] = LabelType::try_from(clustering.max_label() + 1).unwrap();
                }
            });
        });
    });
    (samples, n_clusters)
}

#[no_mangle]
pub extern "C" fn sample_epa(
    n_samples: SEXP,
    similarity: SEXP,
    mass: SEXP,
    discount: SEXP,
    n_cores: SEXP,
) -> SEXP {
    let mut rng = Pcg64Mcg::from_seed(r::random_bytes_from_r::<16>());
    let n_samples = n_samples.as_usize();
    let n_items = similarity.nrow_usize();
    let (samples, _) = sample_epa_engine(
        n_samples,
        n_items,
        similarity.as_double_slice(),
        mass.as_double(),
        discount.as_double(),
        usize::try_from(n_cores.as_integer()).unwrap(),
        &mut rng,
    );
    let n_samples = samples.len() / n_items;
    let result = r::integer_matrix(i32::try_from(n_samples).unwrap(), i32::try_from(n_items).unwrap()).protect();
    let result_slice = result.as_integer_slice_mut();
    for i in 0..n_items {
        for j in 0..n_samples {
            result_slice[i * n_samples + j] = i32::from(samples[j * n_items + i] + 1);
        }
    }
    r::unprotect(1);
    result
}


#[no_mangle]
pub extern "C" fn caviarpd_n_clusters(
    n_samples: SEXP,
    similarity: SEXP,
    mass: SEXP,
    discount: SEXP,
    use_vi: SEXP,
    n_runs: SEXP,
    max_size: SEXP,
    n_cores: SEXP,
) -> SEXP {
    let mut rng = Pcg64Mcg::from_seed(r::random_bytes_from_r::<16>());
    let n_samples = n_samples.as_usize();
    let n_items = similarity.nrow_usize();
    let (samples, n_clusters) = sample_epa_engine(
        n_samples,
        n_items,
        similarity.as_double_slice(),
        mass.as_double(),
        discount.as_double(),
        n_cores.as_usize(),
        &mut rng,
    );
    let n_samples = samples.len() / n_items;
    let clusterings = Clusterings::unvalidated(n_samples, n_items, samples, n_clusters);
    let pdi = PartitionDistributionInformation::Draws(&clusterings);
    let a = 1.0;
    let loss_function = if use_vi.as_bool() {
        LossFunction::VI(a)
    } else {
        LossFunction::BinderDraws(a)
    };
    let p = SALSOParameters {
        n_items,
        max_size: LabelType::try_from(max_size.as_integer()).unwrap(),
        max_size_as_rf: false,
        max_scans: u32::MAX,
        max_zealous_updates: 10,
        n_runs: u32::try_from(n_runs.as_integer()).unwrap(),
        prob_sequential_allocation: 0.5,
        prob_singletons_initialization: 0.0,
    };
    let fit = minimize_by_salso(
        pdi,
        loss_function,
        &p,
        f64::INFINITY,
        u32::try_from(n_cores.as_integer()).unwrap(),
        &mut rng,
    );
    r::integer((fit.clustering.into_iter().max().unwrap() + 1) as i32)
}
