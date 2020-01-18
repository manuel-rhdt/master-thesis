mod configuration;
mod gillespie;
mod likelihood;

use std::env;
use std::fs::{self, File};
use std::io::prelude::*;

use configuration::{calculate_hash, create_dir_if_not_exists, Config};
use gillespie::{SimulationCoordinator, Trajectory, TrajectoryIterator};
use likelihood::log_likelihood;

use ndarray::{Array, Array1, Array2, ArrayView1, Axis};
use ndarray_npy::WriteNpyExt;
use rayon;
use rayon::prelude::*;
use toml;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

fn conditional_likelihood(
    traj_lengths: &[f64],
    num_res: usize,
    coordinator: &mut SimulationCoordinator<impl rand::Rng>,
) -> Array1<f64> {
    let res_network = coordinator.res_network.clone();
    let sig = coordinator.generate_signal().collect();
    let mut result = Array2::zeros((num_res, traj_lengths.len()));
    for mut out in result.outer_iter_mut() {
        let res = coordinator.generate_response(sig.iter());
        for (ll, out) in log_likelihood(traj_lengths, sig.iter(), res, &res_network).zip(&mut out) {
            *out = ll;
        }
    }
    result.mean_axis(Axis(0)).unwrap()
}

fn marginal_likelihood(
    traj_lengths: &[f64],
    signals_pre: &[Trajectory<Vec<f64>, Vec<f64>, Vec<u32>>],
    coordinator: &mut SimulationCoordinator<impl rand::Rng>,
) -> Array1<f64> {
    let sig = coordinator.generate_signal().collect();
    let res = coordinator.generate_response(sig.iter()).collect();

    let mut result = Array2::zeros((signals_pre.len(), traj_lengths.len()));
    for (sig, mut out) in signals_pre.iter().zip(result.outer_iter_mut()) {
        for (ll, out) in log_likelihood(
            traj_lengths,
            sig.iter(),
            res.iter(),
            &coordinator.res_network,
        )
        .zip(&mut out)
        {
            *out = ll;
        }
    }

    result.map_axis(Axis(0), log_mean_exp)
}

pub fn log_mean_exp(values: ArrayView1<f64>) -> f64 {
    use std::ops::Div;

    let max = values.fold(std::f64::NEG_INFINITY, |a, &b| a.max(b));
    values
        .iter()
        .map(|&x| (x - max).exp())
        .sum::<f64>()
        .div(values.len() as f64)
        .ln()
        + max
}

enum EntropyType {
    Conditional,
    Marginal,
}

fn print_help() -> ! {
    println!("\
Usage:
    gillespie [conf]

Arguments:
    conf - The path to the configuration file. The default is 
           'configuration.toml'.
");
    std::process::exit(0)
}

fn main() -> std::io::Result<()> {
    use std::os::unix::fs::OpenOptionsExt;

    let args: Vec<String> = env::args().collect();
    let configuration_filename = match args.len() {
        // no arguments passed
        1 => "configuration.toml",
        // one argument passed
        2 => &args[1],
        _ => {
            print_help()
        }
    };

    let mut config_file = File::open(configuration_filename)?;
    let mut contents = String::new();
    config_file.read_to_string(&mut contents)?;
    let conf: Config = toml::from_str(&contents)?;

    create_dir_if_not_exists(&conf.output)?;
    let info_toml_path = &conf.output.join("info.toml");

    let worker_name = std::env::var("GILLESPIE_WORKER_ID").ok();

    match fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .mode(0o444)
        .open(info_toml_path)
    {
        Ok(mut config_file_copy) => write!(config_file_copy, "{}", contents)?,
        Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => {
            if worker_name.is_none() {
                panic!(
                    "{:?} already exists!\n\
                     to run multiple jobs in parallel set the environment \
                     variable GILLESPIE_WORKER_ID.",
                    info_toml_path
                );
            }

            // sleep a little so the other workers have time to write the info.toml
            std::thread::sleep(std::time::Duration::from_secs(1));

            let mut info_toml = File::open(info_toml_path)?;
            let mut contents = String::new();
            info_toml.read_to_string(&mut contents)?;
            let other_conf: Config = toml::from_str(&contents)?;

            if other_conf == conf {
                println!("configuration.toml matches existing info.toml");
            } else {
                panic!("{:?} does not match configuration.toml", info_toml_path);
            }
        }
        Err(other) => return Err(other),
    }

    let worker_dir = &conf
        .output
        .join(&worker_name.as_deref().unwrap_or("default"));
    fs::create_dir(worker_dir)?;

    let seed_base = conf.get_relevant_hash() ^ calculate_hash(&worker_name);

    let traj_lengths = Array::linspace(0.0, conf.length, 2000);

    let ce_chunks = (0..conf.conditional_entropy.num_signals)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0xabcd_abcd;
            let mut coordinator = conf.create_coordinator(seed);
            -conditional_likelihood(
                traj_lengths.as_slice().unwrap(),
                conf.conditional_entropy.responses_per_signal,
                &mut coordinator,
            )
        })
        .chunks(conf.batch_size)
        .enumerate();

    let mut coordinator = conf.create_coordinator(seed_base);
    let mut signals_pre = Vec::with_capacity(conf.marginal_entropy.num_signals);
    for _ in 0..conf.marginal_entropy.num_signals {
        signals_pre.push(coordinator.generate_signal().collect())
    }

    let me_chunks = (0..conf.marginal_entropy.num_responses)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0x1234_1234;
            let mut coordinator = conf.create_coordinator(seed);
            -marginal_likelihood(
                traj_lengths.as_slice().unwrap(),
                &signals_pre,
                &mut coordinator,
            )
        })
        .chunks(conf.batch_size)
        .enumerate();

    ce_chunks
        .map(|(chunk_num, log_lh)| (EntropyType::Conditional, chunk_num, log_lh))
        .chain(me_chunks.map(|(chunk_num, log_lh)| (EntropyType::Marginal, chunk_num, log_lh)))
        .for_each(|(entropy_type, chunk_num, log_lh)| {
            let mut array = Array2::zeros((conf.batch_size, traj_lengths.len()));
            for (mut row, ll) in array.outer_iter_mut().zip(log_lh.iter()) {
                row.assign(ll);
            }

            let filename = match entropy_type {
                EntropyType::Conditional => format!("ce-{:03}.npy", chunk_num),
                EntropyType::Marginal => format!("me-{:03}.npy", chunk_num),
            };

            if let Ok(out_file) = fs::OpenOptions::new()
                .create_new(true)
                .write(true)
                .mode(0o444)
                .open(worker_dir.join(&filename))
            {
                array.write_npy(out_file).expect("could not write npy");
            } else {
                panic!("could not write chunk {:?}", &filename);
            }
        });

    Ok(())
}
