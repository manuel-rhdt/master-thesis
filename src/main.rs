mod configuration;
mod gillespie;
mod kde;
mod likelihood;

use std::env;
use std::fs::{File, OpenOptions, Permissions};
use std::io::prelude::*;
use std::os::unix::fs::{OpenOptionsExt, PermissionsExt};
use std::sync::{
    atomic::{AtomicBool, Ordering},
    Mutex,
};

use configuration::{calculate_hash, create_dir_if_not_exists, parse_configuration, Config};
use gillespie::{SimulationCoordinator, Trajectory, TrajectoryIterator};
use likelihood::log_likelihood;

use ndarray::{array, Array, Array1, Array2, ArrayView1, Axis};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use rayon;
use rayon::prelude::*;
use serde::Serialize;
use toml;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

static UNCAUGHT_SIGNAL: AtomicBool = AtomicBool::new(false);

macro_rules! check_abort_signal {
    () => {
        check_abort_signal!(Err(
            std::io::Error::from(std::io::ErrorKind::Interrupted).into()
        ));
    };

    ($default:expr) => {
        if UNCAUGHT_SIGNAL.load(Ordering::SeqCst) {
            return $default;
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

    let kde = coordinator.equilibrate_respones_dist(sig.iter().components(), 1000);

    let mut result = Array2::zeros((num_res, traj_lengths.len()));
    for mut out in result.outer_iter_mut() {
        check_abort_signal!(array![]);
        let res = coordinator.generate_response(sig.iter());
        let log_p0 = kde.pdf(res.components()[0]).ln();
        for (ll, out) in log_likelihood(traj_lengths, sig.iter(), res, &res_network).zip(&mut out) {
            *out = log_p0 + ll;
        }
    }
    result.mean_axis(Axis(0)).unwrap()
}

fn marginal_likelihood(
    traj_lengths: &[f64],
    signals_pre: &[(
        Trajectory<Vec<f64>, Vec<f64>, Vec<u32>>,
        kde::NormalKernelDensityEstimate,
    )],
    coordinator: &mut SimulationCoordinator<impl rand::Rng>,
) -> Array1<f64> {
    let sig = coordinator.generate_signal().collect();
    let res = coordinator.generate_response(sig.iter()).collect();

    let mut result = Array2::zeros((signals_pre.len(), traj_lengths.len()));
    for ((sig, kde), ref mut out) in signals_pre.iter().zip(result.outer_iter_mut()) {
        check_abort_signal!(array![]);
        let logp = kde.pdf(res.iter().components()[0]).ln();
        for (ll, out) in log_likelihood(
            traj_lengths,
            sig.iter(),
            res.iter(),
            &coordinator.res_network,
        )
        .zip(out)
        {
            *out = logp + ll;
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

    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let configuration_filename = match args.len() {
        // no arguments passed
        1 => "configuration.toml",
        // one argument passed
        2 => &args[1],
        _ => print_help(),
    };
    let conf = parse_configuration(configuration_filename)?;

    create_dir_if_not_exists(&conf.output)?;
    let file_marginal = Mutex::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .mode(0o444)
            .open(conf.output.join("marginal_entropy.txt"))?,
    );
    let file_conditional = Mutex::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .mode(0o444)
            .open(conf.output.join("conditional_entropy.txt"))?,
    );

    let traj_lengths = Array::linspace(0.0, conf.length, conf.num_trajectory_lengths);

    let worker_name = env::var("GILLESPIE_WORKER_ID").ok();

    let mut worker_info = WorkerInfo {
        hostname: env::var("HOSTNAME").ok(),
        worker_id: worker_name.clone(),
        start_time: chrono::Local::now()
            .to_rfc3339()
            .parse()
            .expect("could not parse current time"),
        end_time: None,
        error: None,
        version: VersionInfo {
            build_time: env!("VERGEN_BUILD_TIMESTAMP").parse().ok(),
            commit_sha: env!("VERGEN_SHA"),
            commit_date: env!("VERGEN_COMMIT_DATE").parse().ok(),
            version: env!("VERGEN_SEMVER"),
        },
        configuration: conf.clone(),
    };
    write!(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(conf.output.join("worker.toml"))?,
        "{}",
        toml::to_string_pretty(&worker_info)?
    )?;

    let seed_base = conf.get_relevant_hash() ^ calculate_hash(&worker_name);

    let ce_chunks = (0..conf.conditional_entropy.num_signals)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0xabcd_abcd;
            let mut coordinator = conf.create_coordinator(seed);
            let before = std::time::Instant::now();
            let ce = -conditional_likelihood(
                traj_lengths.as_slice().unwrap(),
                conf.conditional_entropy.responses_per_signal,
                &mut coordinator,
            );
            log::info!(
                "timing: conditional_entropy: {:.2} s",
                (std::time::Instant::now() - before).as_secs_f64()
            );
            ce
        })
        .map(|log_lh| (EntropyType::Conditional, log_lh));

    let coordinator = conf.create_coordinator(seed_base);
    let before = std::time::Instant::now();
    let signals_pre = (0..conf.marginal_entropy.num_signals)
        .into_par_iter()
        .map_with(coordinator, |coordinator, i| {
            check_abort_signal!();
            coordinator.rng = Pcg64Mcg::seed_from_u64(seed_base ^ i as u64);
            let sig = coordinator.generate_signal().collect();
            let kde = coordinator.equilibrate_respones_dist(sig.iter().components(), 1_000);
            Ok((sig, kde))
        })
        .collect::<Result<Vec<_>, std::io::Error>>()?;
    log::info!(
        "timing: generate_signals: {:.2} s",
        (std::time::Instant::now() - before).as_secs_f64()
    );

    let me_chunks = (0..conf.marginal_entropy.num_responses)
        .into_par_iter()
        .map(|num| {
            let seed = num as u64 ^ seed_base ^ 0x1234_1234;
            let mut coordinator = conf.create_coordinator(seed);
            let before = std::time::Instant::now();
            let me = -marginal_likelihood(
                traj_lengths.as_slice().unwrap(),
                &signals_pre,
                &mut coordinator,
            );
            log::info!(
                "timing: marginal_entropy: {:.2} s",
                (std::time::Instant::now() - before).as_secs_f64()
            );
            me
        })
        .map(|log_lh| (EntropyType::Marginal, log_lh));

    match ce_chunks
        .chain(me_chunks)
        .try_for_each(|(entropy_type, log_lh)| -> std::io::Result<()> {
            check_abort_signal!();
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
        }) {
        Ok(()) => {}
        Err(error) => {
            log::error!("Error: {:?}", error);
            worker_info.error = Some(format!("{}", error))
        }
    }

    worker_info.end_time = Some(
        chrono::Local::now()
            .to_rfc3339()
            .parse()
            .expect("could not parse current time"),
    );
    let mut worker_toml = File::create(conf.output.join("worker.toml"))?;
    write!(worker_toml, "{}", toml::to_string_pretty(&worker_info)?)?;
    worker_toml.set_permissions(Permissions::from_mode(0o444))?;

    Ok(())
}

#[derive(Debug, Clone, Serialize)]
struct WorkerInfo {
    hostname: Option<String>,
    worker_id: Option<String>,
    start_time: toml::value::Datetime,
    end_time: Option<toml::value::Datetime>,
    error: Option<String>,
    version: VersionInfo,
    configuration: Config,
}

#[derive(Debug, Clone, Serialize)]
struct VersionInfo {
    commit_sha: &'static str,
    commit_date: Option<toml::value::Datetime>,
    version: &'static str,
    build_time: Option<toml::value::Datetime>,
}
