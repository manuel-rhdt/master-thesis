use super::gillespie::{Event, ReactionNetwork, TrajectoryIterator};

fn sum_of_reaction_propensities(components: &[f64], reactions: &ReactionNetwork) -> f64 {
    let mut acc = 0.0;
    for j_reaction in 0..reactions.len() {
        let k = reactions.k[j_reaction];
        let reactants = &reactions.reactants[j_reaction];

        let mut tmp = k;
        for &reactant in reactants.iter() {
            tmp *= components[reactant as usize];
        }

        acc += tmp
    }
    acc
}

fn propensity_of_event(
    components: &[f64],
    reaction_event: u32,
    reactions: &ReactionNetwork,
) -> f64 {
    let k = reactions.k[reaction_event as usize];
    let reactants = &reactions.reactants[reaction_event as usize];

    let mut result = k;
    for &reactant in reactants.iter() {
        result *= components[reactant as usize];
    }
    result
}

struct TrajectoryIter<'a> {
    index: usize,
    timestamps: &'a [f64],
    values: &'a [f64],
}

impl<'a> Iterator for TrajectoryIter<'a> {
    type Item = (f64, (f64, f64));

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.timestamps.len() || self.index >= self.values.len() {
            let last_val = self.values[self.values.len() - 1];
            return Some((std::f64::INFINITY, (last_val, last_val)));
        }
        if self.index == 0 {
            self.index = 1;
        }

        let result = (
            self.timestamps[self.index],
            (self.values[self.index - 1], self.values[self.index]),
        );
        self.index += 1;
        Some(result)
    }
}

struct ComponentIter<'r, Signal, Response> {
    signal: Signal,
    response: Response,
    reactions: &'r ReactionNetwork,
    remaining_constant_time: f64,
    components: Vec<f64>,
}

impl<'r, Signal, Response> ComponentIter<'r, Signal, Response>
where
    Signal: TrajectoryIterator,
    Response: TrajectoryIterator<Ex = u32>,
{
    fn new(mut signal: Signal, response: Response, reactions: &'r ReactionNetwork) -> Self {
        let remaining_constant_time = signal
            .advance()
            .map(|event| event.time)
            .unwrap_or(std::f64::INFINITY);
        let num_sig_comp = signal.components().len();
        let num_res_comp = response.components().len();
        ComponentIter {
            signal,
            response,
            reactions,
            remaining_constant_time,
            components: vec![0.0; num_sig_comp + num_res_comp],
        }
    }
}

impl<'r, Signal, Response> Iterator for ComponentIter<'r, Signal, Response>
where
    Signal: TrajectoryIterator,
    Response: TrajectoryIterator<Ex = u32>,
{
    type Item = (f64, f64, f64);

    /// Average of `trajectory` with `old_timestamp` in the time-intervals specified by
    /// `new_timestamps`.
    ///
    /// Note: This function assumes that both `old_timestamps` and `new_timestamps` are
    /// ordered.
    ///
    /// # Discussion
    ///
    /// This function is used for the calculation of the mean propensity over a timestep in
    /// the response trajectory. Since the propensity for a reaction is variable, the
    /// survival probability depends on the time-integrated value of the total propensity.
    /// Since the signal trajectory is assumed to be piecewise constant we can calculate the
    /// average signal between every pair of response timestamps. If we use this average
    /// signal to compute the total propensity, we don't require the calculation of the time
    /// integral anymore.
    ///
    /// This function computes the average of a trajectory given by the values `trajectory`
    /// at the timestamps `old_timestamps` where the averaging happens for the intervals
    /// between pairs of timestamps in `new_timestamps`.
    ///
    /// Returns a list of averages of size `len(new_timestamps) - 1`.
    ///
    /// ```
    ///             |                                    |
    ///             |        +---------------------------| <-- trajectory[k + 1]
    ///             |========|===========================| <== average[i]
    ///             |        |                           |
    ///          ---|--------+                           | <-- trajectory[k]
    ///             | old_timestamps[k + 1]              |
    ///             |                                    |
    ///             +------------------------------------+---> time
    ///     new_timestamps[i]                  new_timestamps[i+1]
    /// ```
    ///
    fn next(&mut self) -> Option<(f64, f64, f64)> {
        let num_ext_comp = self.signal.components().len();
        if let Some(Event { time: delta_t }) = self.response.advance() {
            let mut rem_time = delta_t;
            let mut integrated_propensity = 0.0;
            loop {
                self.components[..num_ext_comp].copy_from_slice(self.signal.components());
                self.components[num_ext_comp..].copy_from_slice(self.response.components());
                if rem_time < self.remaining_constant_time {
                    let reaction_event = self.response.update();
                    let event_prop =
                        propensity_of_event(&self.components, reaction_event, self.reactions);
                    integrated_propensity +=
                        rem_time * sum_of_reaction_propensities(&self.components, self.reactions);
                    break Some((delta_t, event_prop, integrated_propensity / delta_t));
                } else {
                    rem_time -= self.remaining_constant_time;
                    integrated_propensity += self.remaining_constant_time
                        * sum_of_reaction_propensities(&self.components, self.reactions);
                    self.signal.update();
                    if let Some(Event { time }) = self.signal.advance() {
                        self.remaining_constant_time = time;
                    } else {
                        self.remaining_constant_time = std::f64::INFINITY;
                    }
                }
            }
        } else {
            None
        }
    }
}

/// Returns an Iterator over pairs of timestamps and log_likelihoood values.
///
fn log_likelihood_inner<'a>(
    signal: impl 'a + TrajectoryIterator<Ex = u32>,
    response: impl 'a + TrajectoryIterator<Ex = u32>,
    reactions: &'a ReactionNetwork,
) -> impl 'a + Iterator<Item = (f64, f64)> {
    let component_iter = ComponentIter::new(signal, response, reactions);
    component_iter.scan(
        (0.0, 0.0),
        |(current_time, ll), (dt, inst_rate, avrg_rate)| {
            *current_time += dt;
            *ll += inst_rate.ln() - avrg_rate * dt;
            Some((*current_time, *ll))
        },
    )
}

pub fn log_likelihood<'a>(
    traj_lengths: &'a [f64],
    signal: impl 'a + TrajectoryIterator<Ex = u32>,
    response: impl 'a + TrajectoryIterator<Ex = u32>,
    reactions: &'a ReactionNetwork,
) -> impl 'a + Iterator<Item = f64> {
    let mut ll_iter = log_likelihood_inner(signal, response, reactions).peekable();

    let bin_iter = traj_lengths.iter();

    bin_iter.scan(0.0, move |acc, &bin_edge| {
        while let Some(&(t, lh)) = ll_iter.peek() {
            if t < bin_edge {
                *acc = lh;
                ll_iter.next();
            } else {
                return Some(*acc);
            }
        }
        None
    })
}
