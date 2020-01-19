mod configuration;
mod gillespie;
mod likelihood;

use std::env;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    Mutex,
};

use configuration::{calculate_hash, create_dir_if_not_exists, Config};
use gillespie::{SimulationCoordinator, Trajectory, TrajectoryIterator};
use likelihood::log_likelihood;

use ndarray::{Array, Array1, Array2, ArrayView1, Axis};
use rayon;
use rayon::prelude::*;
use serde::Serialize;
use toml;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

static UNCAUGHT_SIGNAL: AtomicBool = AtomicBool::new(false);

macro_rules! check_abort_signal {
    () => {
        if UNCAUGHT_SIGNAL.compare_and_swap(true, false, Ordering::SeqCst) {
            panic!("Process received SIGINT or SIGTERM!");
        }
    };
}

fn conditional_likelihood(
    traj_lengths: &[f64],
    num_res: usize,
    coordinator: &mut SimulationCoordinator<impl rand::Rng>,
) -> Array1<f64> {
    let res_network = coordinator.res_network.clone();
    let sig = coordinator.generate_signal().collect();
    let mut result = Array2::zeros((num_res, traj_lengths.len()));
    for mut out in result.outer_iter_mut() {
        check_abort_signal!();
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
        check_abort_signal!();
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
    println!(
        "\
Usage:
    gillespie [conf]

Arguments:
    conf - The path to the configuration file. The default is 
           'configuration.toml'.
"
    );
    std::process::exit(0)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    ctrlc::set_handler(move || {
        UNCAUGHT_SIGNAL.store(true, Ordering::SeqCst);
    })
    .expect("Error setting Ctrl-C handler");

    let args: Vec<String> = env::args().collect();
    let configuration_filename = match args.len() {
        // no arguments passed
        1 => "configuration.toml",
        // one argument passed
        2 => &args[1],
        _ => print_help(),
    };

    let mut config_file = File::open(configuration_filename)?;
    let mut contents = String::new();
    config_file.read_to_string(&mut contents)?;
    let conf: Config = toml::from_str(&contents)?;

    create_dir_if_not_exists(&conf.output)?;
    let file_marginal = Mutex::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(conf.output.join("marginal_entropy.txt"))?,
    );
    let file_conditional = Mutex::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(conf.output.join("conditional_entropy.txt"))?,
    );

    let traj_lengths = Array::linspace(0.0, conf.length, conf.num_trajectory_lengths);

    let worker_name = env::var("GILLESPIE_WORKER_ID").ok();

    let worker_info = WorkerInfo {
        hostname: env::var("HOSTNAME").ok(),
        worker_id: worker_name.clone(),
        start_time: chrono::Local::now()
            .to_rfc3339()
            .parse()
            .expect("could not parse current time"),
        version: VersionInfo {
            build_time: env!("VERGEN_BUILD_TIMESTAMP").parse().ok(),
            commit_sha: env!("VERGEN_SHA"),
            commit_date: env!("VERGEN_COMMIT_DATE").parse().ok(),
            version: env!("VERGEN_SEMVER"),
        },
    };
    write!(
        File::create(conf.output.join("worker.toml"))?,
        "{}",
        toml::to_string_pretty(&worker_info)?
    )?;

    let seed_base = conf.get_relevant_hash() ^ calculate_hash(&worker_name);

    let ce_chunks = (0..conf.conditional_entropy.num_signals)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0xabcd_abcd;
            let mut coordinator = conf.create_coordinator(seed);
            (
                num,
                -conditional_likelihood(
                    traj_lengths.as_slice().unwrap(),
                    conf.conditional_entropy.responses_per_signal,
                    &mut coordinator,
                ),
            )
        });

    check_abort_signal!();
    let mut coordinator = conf.create_coordinator(seed_base);
    let mut signals_pre = Vec::with_capacity(conf.marginal_entropy.num_signals);
    for _ in 0..conf.marginal_entropy.num_signals {
        signals_pre.push(coordinator.generate_signal().collect())
    }
    check_abort_signal!();

    let me_chunks = (0..conf.marginal_entropy.num_responses)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0x1234_1234;
            let mut coordinator = conf.create_coordinator(seed);
            (
                num,
                -marginal_likelihood(
                    traj_lengths.as_slice().unwrap(),
                    &signals_pre,
                    &mut coordinator,
                ),
            )
        });

    ce_chunks
        .map(|(row, log_lh)| (EntropyType::Conditional, row, log_lh))
        .chain(me_chunks.map(|(row, log_lh)| (EntropyType::Marginal, row, log_lh)))
        .panic_fuse()
        .try_for_each(|(entropy_type, _, log_lh)| -> std::io::Result<()> {
            let mut buffer = ryu::Buffer::new();
            let mut line = String::with_capacity(log_lh.dim() * 22);
            for &val in &log_lh {
                line += buffer.format(val);
                line += " ";
            }
            let mut file = match entropy_type {
                EntropyType::Conditional => file_conditional.lock().unwrap(),
                EntropyType::Marginal => file_marginal.lock().unwrap(),
            };
            writeln!(file, "{}", line)?;
            Ok(())
        })?;

    Ok(())
}

#[derive(Debug, Clone, Serialize)]
struct WorkerInfo {
    hostname: Option<String>,
    worker_id: Option<String>,
    start_time: toml::value::Datetime,
    version: VersionInfo,
}

#[derive(Debug, Clone, Serialize)]
struct VersionInfo {
    commit_sha: &'static str,
    commit_date: Option<toml::value::Datetime>,
    version: &'static str,
    build_time: Option<toml::value::Datetime>,
}
