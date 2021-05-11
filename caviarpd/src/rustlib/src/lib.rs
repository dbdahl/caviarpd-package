use dahl_salso::clustering::Clusterings;
use dahl_salso::optimize::{minimize_by_salso, SALSOParameters};
use dahl_salso::{LabelType, LossFunction, PartitionDistributionInformation};
use epa::epa::{sample, EpaParameters, SquareMatrixBorrower};
use epa::perm::Permutation;
use extendr_api::prelude::*;
use rand::thread_rng;

fn sample_epa_engine(
    n_samples: i32,
    similarity: &Robj,
    mass: f64,
    discount: f64,
    n_cores: i32,
) -> (Vec<LabelType>, Vec<LabelType>) {
    let n_samples = n_samples as usize;
    let n_items = similarity.nrows();
    let n_cores = if n_cores == 0 {
        num_cpus::get()
    } else {
        n_cores as usize
    };
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
            plan.push((left1, left2));
            stick1 = right1;
            stick2 = right2;
        }
        plan.push((stick1, stick2));
        let sim = SquareMatrixBorrower::from_slice(similarity.as_real_slice().unwrap(), n_items);
        plan.into_iter().for_each(|p| {
            s.spawn(move |_| {
                let mut rng = thread_rng();
                let mut params =
                    EpaParameters::new(sim, Permutation::natural(n_items), mass, discount).unwrap();
                for i in 0..n_samples_per_core {
                    params.shuffle_permutation(&mut rng);
                    let clustering = sample(&params, &mut rng);
                    let zero: LabelType = 0;
                    clustering.relabel_into_slice(zero, &mut p.0[i * n_items..(i + 1) * n_items]);
                    p.1[i] = (clustering.max_label() + 1) as LabelType;
                }
            });
        });
    });
    (samples, n_clusters)
}

#[extendr]
fn sample_epa(
    n_samples: i32,
    similarity: Robj,
    mass: f64,
    discount: f64,
    n_cores: i32,
) -> RMatrix<Vec<i32>> {
    let (samples, _) = sample_epa_engine(n_samples, &similarity, mass, discount, n_cores);
    let n_samples = n_samples as usize;
    let n_items = similarity.nrows();
    let mut result = Vec::with_capacity(n_samples * n_items);
    for i in 0..n_items {
        for j in 0..n_samples {
            result.push((samples[j * n_items + i] + 1) as i32);
        }
    }
    RMatrix::new(result, n_samples, n_items)
}

#[extendr]
fn caviarpd(
    n_samples: i32,
    similarity: Robj,
    mass: f64,
    discount: f64,
    use_vi: bool,
    n_runs: i32,
    max_size: i32,
    n_cores: i32,
) -> Vec<i32> {
    let (samples, n_clusters) = sample_epa_engine(n_samples, &similarity, mass, discount, n_cores);
    let n_items = similarity.nrows();
    let n_samples = samples.len() / n_items;
    let mut rng = thread_rng();
    let clusterings = Clusterings::unvalidated(n_samples, n_items, samples, n_clusters);
    let pdi = PartitionDistributionInformation::Draws(&clusterings);
    let a = 1.0;
    let loss_function = if use_vi {
        LossFunction::VI(a)
    } else {
        LossFunction::BinderDraws(a)
    };
    let p = SALSOParameters {
        n_items,
        max_size: max_size as LabelType,
        max_size_as_rf: false,
        max_scans: u32::MAX,
        max_zealous_updates: 10,
        n_runs: n_runs as u32,
        prob_sequential_allocation: 0.5,
        prob_singletons_initialization: 0.0,
    };
    let fit = minimize_by_salso(
        pdi,
        loss_function,
        &p,
        f64::INFINITY,
        n_cores as u32,
        &mut rng,
    );
    fit.clustering.into_iter().map(|x| (x+1) as i32).collect()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod caviarpd;
    fn sample_epa;
    fn caviarpd;
}
