mod gillespie;
mod likelihood;

use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::path::PathBuf;

use gillespie::{ReactionNetwork, SimulationCoordinator, Trajectory, TrajectoryIterator};
use likelihood::log_likelihood;

use ndarray::{Array, Array1, Array2, ArrayView1, Axis};
use ndarray_npy::WriteNpyExt;
use rand::{self, SeedableRng};
use rand_chacha::ChaChaRng;
use rayon;
use rayon::prelude::*;
use serde::Deserialize;
use toml;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

fn calculate_hash<T: Hash + ?Sized>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
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
        let res = coordinator.generate_response(sig.iter());
        for (ll, out) in
            log_likelihood(traj_lengths, sig.iter(), res, &res_network).zip(out.iter_mut())
        {
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
    let res_network = coordinator.res_network.clone();
    let sig = coordinator.generate_signal().collect();
    let res = coordinator.generate_response(sig.iter()).collect();

    let mut result = Array2::zeros((signals_pre.len(), traj_lengths.len()));
    for (sig, mut out) in signals_pre.iter().zip(result.outer_iter_mut()) {
        for (ll, out) in
            log_likelihood(traj_lengths, sig.iter(), res.iter(), &res_network).zip(out.iter_mut())
        {
            *out = ll;
        }
    }

    result.map_axis(Axis(0), log_mean_exp)
}

pub fn log_mean_exp(values: ArrayView1<f64>) -> f64 {
    use std::ops::Div;

    let max = values.fold(std::f64::NEG_INFINITY, |a, &b| a.max(b));
    let count_non_nan = values.iter().filter(|x| !x.is_nan()).count() as f64;
    values
        .iter()
        .filter(|x| !x.is_nan())
        .map(|&x| (x - max).exp())
        .sum::<f64>()
        .div(count_non_nan)
        .ln()
        + max
}

#[derive(Deserialize, Debug, Clone, PartialEq)]
struct Config {
    output: PathBuf,
    batch_size: usize,
    length: f64,
    conditional_entropy: ConfigConditionalEntropy,
    marginal_entropy: ConfigMarginalEntropy,
    signal: ConfigReactionNetwork,
    response: ConfigReactionNetwork,
}

impl Config {
    pub fn hash_relevant<H: Hasher>(&self, hasher: &mut H) {
        self.conditional_entropy.hash(hasher);
        self.marginal_entropy.hash(hasher);
        self.length.to_bits().hash(hasher);
        self.signal.hash(hasher);
        self.response.hash(hasher);
    }

    pub fn get_relevant_hash(&self) -> u64 {
        let mut s = DefaultHasher::new();
        self.hash_relevant(&mut s);
        s.finish()
    }

    pub fn create_coordinator(&self, seed: u64) -> SimulationCoordinator<ChaChaRng> {
        let rng = ChaChaRng::seed_from_u64(seed);
        SimulationCoordinator {
            trajectory_len: self.length,

            sig_network: self.signal.to_reaction_network(),
            res_network: self.response.to_reaction_network(),

            rng,
        }
    }
}

#[derive(Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
struct ConfigConditionalEntropy {
    num_signals: usize,
    responses_per_signal: usize,
}

#[derive(Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
struct ConfigMarginalEntropy {
    num_signals: usize,
    num_responses: usize,
}

#[derive(Deserialize, Debug, Clone)]
struct ConfigReactionNetwork {
    initial: f64,
    components: Vec<String>,
    reactions: Vec<Reaction>,
}

impl PartialEq for ConfigReactionNetwork {
    fn eq(&self, right: &ConfigReactionNetwork) -> bool {
        self.initial.to_bits() == right.initial.to_bits()
            && self.components == right.components
            && self.reactions == right.reactions
    }
}

impl Hash for ConfigReactionNetwork {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.initial.to_bits().hash(state);
        self.components.hash(state);
        self.reactions.hash(state);
    }
}

impl ConfigReactionNetwork {
    fn to_reaction_network(&self) -> ReactionNetwork {
        let components: HashSet<_> = self.components.iter().collect();
        let mut external_components = HashSet::new();
        for reaction in &self.reactions {
            for name in reaction.reactants.iter().chain(reaction.products.iter()) {
                if !components.contains(name) {
                    external_components.insert(name);
                }
            }
        }
        let num_ext_components = external_components.len();

        let mut name_table = HashMap::new();
        name_table.extend(
            external_components
                .iter()
                .enumerate()
                .map(|(value, &key)| (key, value)),
        );
        name_table.extend(
            components
                .iter()
                .enumerate()
                .map(|(value, &key)| (key, value + num_ext_components)),
        );

        let mut k = Vec::with_capacity(self.reactions.len());
        let mut reactants = Vec::with_capacity(self.reactions.len());
        let mut products = Vec::with_capacity(self.reactions.len());

        for reaction in &self.reactions {
            k.push(reaction.k);
            let r = reaction
                .reactants
                .iter()
                .map(|reactant| name_table[reactant] as u32)
                .collect();
            reactants.push(r);
            let p = reaction
                .products
                .iter()
                .map(|product| name_table[product] as u32)
                .collect();
            products.push(p);
        }

        ReactionNetwork {
            k,
            reactants,
            products,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
struct Reaction {
    k: f64,
    reactants: Vec<String>,
    products: Vec<String>,
}

impl PartialEq for Reaction {
    fn eq(&self, right: &Reaction) -> bool {
        self.k.to_bits() == right.k.to_bits()
            && self.reactants == right.reactants
            && self.products == right.products
    }
}

impl Hash for Reaction {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.k.to_bits().hash(state);
        self.reactants.hash(state);
        self.products.hash(state);
    }
}

enum EntropyType {
    Conditional,
    Marginal,
}

fn create_dir_if_not_exists<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<()> {
    match fs::create_dir(path) {
        Ok(()) => Ok(()),
        Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => Ok(()),
        other => other,
    }
}

fn main() -> std::io::Result<()> {
    use std::os::unix::fs::OpenOptionsExt;

    let mut config_file = File::open("configuration.toml")?;
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