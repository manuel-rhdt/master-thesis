---
bibliography: ["library.bib"]
autoEqnLabels: true
link-citations: true
linkReferences: true
cref: true
---

# Monte Carlo Estimate of the Mutual Information

Equipped with analytical formulae for the computation of trajectory probabilities and with methods for efficient stochastic simulation, we can start to develop estimates for the mutual information. The basis for our method is @eq:mi_form2; specifically we separate the mutual information into two parts that are computed independently, the _marginal entropy_ $\mathrm H(\mathcal X)$ and the _conditional entropy_ $\mathrm H(\mathcal X|\mathcal S)$. Since signals and responses are time-varying, both entropies involve integrals over spaces of trajectories which are high-dimensional. At high dimensionality, direct numerical integration is not viable and instead, we use Monte-Carlo approaches based on random sampling of signals and responses.

While Monte-Carlo methods comprise a wide variety of approaches to stochastically evaluate integrals or sums the common idea is easily stated. We have a state space $U$ and a probability distribution $p_U$ on that state space. The problem is to numerically evaluate integrals of the form
$$
F = \langle f(u) \rangle = \int\limits_{U} \mathrm du\ p_U(u)\; f(u)
$$
where $f: U\rightarrow\mathbb R$ is some function of interest. If $U$ is high-dimensional it is very time-consuming to estimate it by direct numerical integration. Instead, we generate random samples $u_1,u_2,\ldots$ from the probability distribution $p_U$ such that by the _law of large numbers_ we have the equality
$$
F = \lim_{N\rightarrow\infty} \frac{\sum^N_{i=1} f(u_i)}{N}
$$ {#eq:inf_mc_estimate}
i.e. the sample average of random samples converges towards their mean. In practice, we can't evaluate the limit in @eq:inf_mc_estimate and we approximate $F$ using a finite number $N$ of random samples
$$
\hat{F} = \frac{\sum^N_{i=1} f(u_i)}{N}\,.
$$
The variance of Monte-Carlo estimates typically decreases much faster than TODO. Monte-Carlo methods depend on the efficient and correct sampling of the corresponding probability distribution. In this thesis ...

In this chapter we will show how to use Monte-Carlo integration to compute the mutual information between stochastic trajectories. The computation requires the estimation of two quantities, marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$. We will show how to perform Monte-Carlo estimates for both entropies, however as we will describe the marginal entropy turns out to be computationally much more difficult to estimate. We conclude the chapter with a thorough analysis of the challenges that arise by ...

## Monte-Carlo Estimate for the Marginal Entropy

We start by considering an abstract system, consisting of two random variables $\mathcal S$ and $\mathcal X$ representing the possible signals and responses, respectively. We will assume the availability of efficient computational methods to generate random samples from either of these random variables. Later in this chapter we will describe some methods for random sample generation. In this subsection we show how to use Monte-Carlo techniques to estimate the marginal entropy $\mathrm H(\mathcal X)$ while in the following subsection we describe a similar method to estimate the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$. 

We intend to compute the marginal entropy $\mathrm H(\mathcal X)$ using Monte-Carlo (MC) sampling to evaluate the necessary integrals. First, given a number of random samples $(\mathbf x_i)_{i=1,\ldots,N_x}$ taken from $\mathcal X$ we propose as an estimate for the entropy
$$
\mathrm H(\mathcal X) = -\int\mathrm d\mathbf x\ \mathrm P(\mathbf x)\ln \mathrm P(\mathbf x) \approx - \frac{\sum\limits_{i=1}^{N_x} \ln\mathrm P(\mathbf x_i)}{N_x} \,.
$$ {#eq:mc_entropy}
I.e. using numerous response samples, we compute the sample average of their individual log-probabilities $\ln\mathrm P(\mathbf x_i)$. As explained in the following section, from a biochemical of a cell we are not able to directly compute $\mathrm P(\mathbf x_i)$. We do however have an efficient way of _generating_ the responses $\mathbf x_1,\mathbf x_2,\ldots$ using the stochastic simulation algorithm.

Nonetheless, we see from @eq:mc_entropy that we _do_ have to evaluate $\mathrm P(\mathbf x_i)$ for every generated sample. To solve this problem we choose to evaluate $\mathrm P(\mathbf x_i)$ by doing another nested Monte-Carlo integration using signal samples $(\mathbf s_j)_{j=1,\ldots,N_s}$ from $\mathcal S$ to write
$$
\mathrm P(\mathbf x_i) = \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\ \mathrm P(\mathbf x_i|\mathbf s) \approx \frac{\sum\limits_{j=1}^{N_s} \mathrm P(\mathbf x_i | \mathbf s_j)}{N_s} \,.
$$ {#eq:mc_marginal}
Hence, for every response sample $\mathbf x_i$ we have to perform a sample average over the conditional probabilities $\mathrm P(\mathbf x_i|\mathbf s_j)$. Again, anticipating results from the following section, we make use of the fact that a stochastic model for a biochemical network allows us to efficiently compute $\mathrm P(\mathbf x_i|\mathbf s_j)$. Therefore it seems viable to use @eq:mc_marginal to compute the marginal probability for a given response.

While for a low-dimensional signal space it is feasible to instead compute the marginalization integral @eq:mc_marginal using direct evaluation [@2019.Cepeda-Humerez] we choose to use a nested MC simulation to also be able to handle high-dimensional signal spaces. This is crucial since time-varying signals are described by high-dimensional trajectories.

We can summarize the estimation procedure for the marginal entropy using the equation
$$
\mathrm H(\mathcal X) = -\left\langle \ln \left\langle \mathrm P(\mathbf x | \mathbf s) \right\rangle_{\mathrm P(\mathbf s)} \right\rangle_{\mathrm P(\mathbf x)}
$$ {#eq:mc_entropy_notation}
where we use the notation $\langle f(x)\rangle_{g(x)}$ for the expected value of $f(x)$ when $x$ is distributed according to the probability density given by $g(x)$. Thus when thinking in mathematical terms we have the shorthand $\langle f(x)\rangle_{g(x)} \equiv\int \mathrm dx\ g(x) f(x)$. We can also easily translate this notation into a Monte-Carlo estimate, i.e. $\langle f(x)\rangle_{g(x)} = \lim\limits_{N\rightarrow\infty}\frac{\sum_{i=1}^N f(x_i)}{N}$ where $x_1, x_2,\ldots$ are independent samples of the probability distribution given by $g(x)$.

## Estimating the Conditional Entropy

We can also estimate the _conditional entropy_ using MC averages over trajectories. We express the conditional entropy using the notation introduced in @eq:mc_entropy_notation
$$
\mathrm H(\mathcal X|\mathcal S) = -\iint \mathrm d\mathbf s\mathrm d\mathbf x\ \mathrm P(\mathbf s)\mathrm P(\mathbf x | \mathbf s) \ln\mathrm P(\mathbf x|\mathbf s) = -\left\langle\langle\ln\mathrm P(\mathbf x | \mathbf s)\rangle_{\mathrm P(\mathbf x | \mathbf s)} \right\rangle_{\mathrm P(\mathbf s)}
$$
to show that we require nested Monte Carlo integrations to evaluate the integral. We first generate signal samples $\mathbf s_1, \ldots, \mathbf s_{N_s}$ from the density $\mathrm P(\mathbf s)$. Let $\mathbf x_i^1,\ldots,\mathbf x_i^{N_x}$ be response samples generated from $\mathrm P(\mathbf x | \mathbf s_i)$. The Monte Carlo estimate for the conditional entropy then reads
$$
\mathrm H(\mathcal X|\mathcal S) \approx - \frac1{N_s N_x} \sum\limits_{i=1}^{N_s} \sum\limits_{j=1}^{N_x} \ln\mathrm P(\mathbf x_i^j | \mathbf s_i)\,.
$$ {#eq:conditional_entropy_estimate}

Using Monte-Carlo computations we can in principle compute the mutual information between arbitrary random variables.

## Monte-Carlo Simulations for Trajectories

- Idea: Use SSA to generate response trajectories for given signals
- The signals themselves are taken to be realizations of a given stochastic process
- Using a set of predefined signals, we compute many responses
- from a set of generated responses we can estimate the conditional entropy
- and with some additional work the marginal entropy

### General Stochastic Dynamics of Signals

In the previous section we found that both, the Monte-Carlo computation of the marginal entropy and of the conditional entropy require the generation of stochastic samples of the signal. This is not surprising since the stochastic dynamics of the signal reflect the amount of information that it contains. For example, if the signal is governed by mostly noise, then even the most efficient biochemical network could not extract a lot of meaningful information out of it. Therefore, to compute the mutual information between the signal and the responses, apart from modeling the biochemical network that processes the signal we also have to model the signal-generating process.

Generally, any kind of noisy signal can be modeled by a _stochastic process_. In many physically relevant cases, these processes arise as the solutions of SDEs that describe the underlying deterministic physics of the signal together with a _noise term_ that gives rise to the stochastic nature of the signal trajectories. An example for such a process is the _Ornstein-Uhlenbeck process_ that describes the random motion of a diffusing particle in a harmonic potential well @2009.Gardiner. It represents a combination of a deterministic harmonic oscillator together with fluctuations arising from the diffusion process. For the kinds of stochastic processes that are described by a SDE there exist many numerical methods to generate approximate signals trajectories. A simple, yet effective method to generate stochastic realizations of SDEs is the _Euler–Maruyama method_ that is a generalization of the Euler method for integrating ODEs @1992.Kloeden. Using such methods we can in principle choose any integrable SDE as the signal-generating process for the computation of the mutual information.

If we want to investigate the simple reaction network for gene expression given in @eq:simple_reaction_network we find that the signal dynamics are themselves described by a reaction network. In that case it is possible to generate _exact_ stochastic realizations of signal trajectories using a stochastic simulation algorithm as discussed in @sec:ssa.

### Generating Responses for Time-Varying Signals

Both, the estimate for the marginal entropy $\mathrm H(\mathcal X)$ and the estimate for the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$ require the generation of appropriate responses to a given time-varying signal. To find stochastically correct responses for a given signal we have to understand how the biochemical network interacts with the signal. In many cases the signal itself is a molecular species that participates in the reaction network. In this scenario, an abundance of signal directly leads to an increased rate of all reactions in which the signal molecule acts as a reactant. In other cases the signal may be a thermodynamic quantity as for example the environmental temperature. Changes in temperature usually effect the rates of _all_ reactions in a chemical reaction network. Thus, in general we think of a signal as a quantity whose value affects some or all of the reaction rates of the biochemical network at a given instant. Since we model the signal as a time-varying process, for a given realization of the signal we can describe the corresponding responses using a master equation with time-dependent reaction rates.

We describe the coupling of the biochemical network to the signal through an explicit dependence of the reaction rates from @eq:general_master_eq on the signal level $s(t)$. Thus we write the master equation for the responses as
$$
\frac{\partial \mathrm P(x, t|x_0, t_0, \mathbf s)}{\partial t} = \sum\limits_{x^\prime\in\mathcal U} \left[
w_t(x, x^\prime, s(t))\ \mathrm P(x^\prime, t|x_0, t_0, \mathbf s)
- w_t(x^\prime, x, s(t))\ \mathrm P(x, t|x_0, t_0, \mathbf s)
\right]
$$ {#eq:general_master_eq2}
which shows the explicit dependence of the stochastic dynamics on the given signal trajectory $\mathbf s$. Since we assume that the reaction rates are the _only_ way in which the signal trajectory affects the response we conclude that stochastic realizations of the process described by @eq:general_master_eq2 are distributed according to the conditional distribution $\mathrm P(\mathbf x|\mathbf s)$.

The formulation of the master equation with an explicit dependence on the signals makes it possible to perform the Monte-Carlo estimates as described above. Specifically, by making use of results derived in the previous chapter we can perform three crucial operations needed for the computational estimate of the mutual information.

#### Generation of response trajectories according to  $\mathrm P(\mathbf x|\mathbf s)$

As already described above, for any signal trajectory $\mathbf s$ we can formulate the master equation as written in @eq:general_master_eq2 and use the stochastic simulation algorithm for time-dependent reaction rates to generate correctly distributed response trajectories.

#### Generation of response trajectories according to  $\mathrm P(\mathbf x)$

Similarly we can generate independent draws from the distribution $\mathrm P(\mathbf x)$ by first generating a stochastic realization $\mathbf s$ of the signal from the corresponding stochastic process and then generating exactly one response according to $\mathrm P(\mathbf x|\mathbf s)$.

#### Compuation of $\mathrm P(\mathbf x|\mathbf s)$

Finally, for any signal-response pair of trajectories, we can compute the conditional probability $\mathrm P(\mathbf x|\mathbf s)$ using the formulae for the trajectory probability derived before in [@eq:trajectory_probability_product;@eq:transition_probability] and inserting the reaction rates for the given signal trajectory $\mathbf s$.

- Note that the definition of the signal and responses as done in this section implies a strict separation of signal and response (make figure?). We will discuss this issue at the end of this thesis in section TODO



<!-- ## Practical Concerns

In the previous sections we explained a) how to use Monte-Carlo integration to compute the mutual information from stochastic realizations of signals and responses and b) how to use the biochemical network model to compute the (conditional) probabilities of the respective trajectories. In this section we discuss some small but important details for the computation of the mutual information that were not mentioned so far. -->

### Probability Distribution of the Initial State {#sec:initial_condition}

We have shown that we can use the stochastic simulation algorithm (SSA) to simulate response trajectories for a given signal. Note however, that the SSA needs us to specify the initial copy numbers of all species in the biochemical network. From this initial state we can then evolve the stochastic dynamics of the reactions in time. Hence, we have to think about the correct choice of the initial condition. In some cases it might be reasonable to have the same predefined initial condition for all responses that are simulated, e.g. in the beginning there is no response $x(t_0) = 0$. In many cases however, the initial states $x_0$ will be distributed according to some probability density $\mathrm P(x_0,t_0)$. A common example are systems that start out in a _steady state_, i.e. at $t_0$ both the signal and the response are distributed according to their joint stationary distribution $p_{S,X}(s,x)$.

In that case, to set up the initial conditions correctly, we start by generating one signal trajectory whose initial point $s_0$ is drawn from the signal's stationary distribution $p_S(s_0)$. Thus, to ensure that the initial condition is drawn from $p_{S,X}(s,x)$, we need to pick the initial condition of the response $x_0$ from the conditional distribution $p_{X|S}(x_0|s_0)$. For our simulations, we estimate that conditional distribution in a separate computation for every signal that is generated.

If the signal dynamics in steady state are time-reversible, we can estimate the distribution function $p_{X|S}(x_0|s_0)$ for a given $s_0$ by generating time-reversed signal trajectories. The idea is to generate $N_\text{past}$ signal trajectories of duration $t_\text{past}$, _backwards in time_, starting from $s_0$ at time $t_0$ that represent possible "pasts" of the system.  For every such trajectory we then generate a response trajectory of the same duration _forward in time_. That is, the response starts at time $t_0-t_\text{past}$ with some arbitrary initial state $x_\text{past}$ and continues until time $t_0$. Since we generate $N_\text{past}$ signals into the past, we end up with $N_\text{past}$ response trajectories, whose end points at $t_0$ are distributed according to $p_{X|S}(x_0|s_0)$.

From these points $x^{(1)}_0, \ldots, x^{(N_\text{past})}_0$ we estimate the distribution function $p_{X|S}(x_0|s_0)$ using Gaussian kernel density estimation (KDE) [@1956.Rosenblatt;@1962.Parzen]. The KDE yields a smooth estimate of the conditional distribution $p_{X|S}(x_0|s_0)$ from the sample points. This estimated distribution is used in two different ways
1. TODO
2. TODO


In summary, as long as we can formulate a stochastic process that describes the signal dynamics _and_ we understand how the signal level affects the reaction rates of the biochemical network, we can in principle use the Monte-Carlo procedure described in the previous section to estimate the mutual information. From the signal's SDE we generate stochastic realizations using appropriate numerical integrators for SDEs such as the _Euler–Maruyama method_ or the _stochastic simulation algorithm_. For every signal realization we have shown how to use the SSA to generate corresponding response trajectories. In the following section we will use these methods on a simple model of gene expression to see how the method is implemented in practice.

<!-- ### Computation in Logarithmic Space

- define likelihood
- there can arise numerical issues for the computation of the likelihoods
- since the computation of the likelihood is a large product (#factors grow with trajectory length), on average the likelihood scales exponentially with trajectory length
- very quickly, we find trajectory lengths were the likelihood is outside the representable range of floating point numbers
- given the exponential nature of the likelihood we find, it is easier to instead compute the log-likelihood $\ln\mathrm P(\mathbf x|\mathbf s)$

For increasingly long trajectories this quantity will in many physically relevant cases either grow or decay exponentially (*TODO: explain why*). Thus sufficiently long trajectories, the numerical values of the likelihood will not be directly representable by conventional floating-point numbers.

This problem can be avoided if we compute the *log-likelihood* $\ell(X|S) \equiv \ln\mathrm P(X|S)$ instead. We can easily rephrase the equation for the likelihood:

$$
\ell(X|S) = \ln\left[ \mathrm P(x_0,t_0 | S) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S) \right] = \ln \mathrm P(x_0,t_0 | S) +\sum\limits^{N-1}_{n=1} \ln \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S)\,.
$$

In practice (due to limited precision of floating-point arithmetic) it is only possible to evaluate the log-likelihood $\ell(X|S) \equiv \ln\mathrm P(X|S)$. This means that the calculation of the averaged likelihood involves the quantity

$$
\ln \sum^{N_S}_{i=1} \mathrm P(X|S^{(i)}) = \ln \sum\limits^{N_S}_{i=1} \exp \ell(X|S^{(i)}) \equiv \mathrm{LSE}\left( \ell(X|S^{(1)}),\ldots, \ell(X|S^{(N_S)})\right)
$$

where $\mathrm{LSE} : \mathbb{R}^n \rightarrow \mathbb{R}$ is called log-sum-exp \cite{blanchard-2019}. An interesting property of $\mathrm{LSE}$ is that it's a smooth approximation to the $\max$ function. This means that for finite sample sizes the monte-carlo estimate of the averaged likelihood will always be too small!


### Estimating the marginal probability of response trajectories

To calculate the mutual information between trajectories we need to have a good estimate for $\ln\left\langle \mathrm P(X | S) \right\rangle_\mathcal{S}$. We calculate this average by sampling of trajectories $(S^{(i)})_{i=1\ldots N_S}$ from the probability distribution of $\mathcal{S}$:

$$
\ln\left\langle\mathrm P(X | S) \right\rangle_\mathcal{S} \approx \ln \frac{\sum^{N_S}_{i=1} \mathrm P(X|S^{(i)})}{N_S} = \ln \sum^{N_S}_{i=1} \mathrm P(X|S^{(i)}) - \ln N_S
$$

Thus we find that it is enough to be able to compute the likelihood between trajectories to estimate the marginal distribution of trajectories.

In practice (due to limited precision of floating-point arithmetic) it is only possible to evaluate the log-likelihood $\ell(X|S) \equiv \ln\mathrm P(X|S)$. This means that the calculation of the averaged likelihood involves the quantity

$$
\ln \sum^{N_S}_{i=1} \mathrm P(X|S^{(i)}) = \ln \sum\limits^{N_S}_{i=1} \exp \ell(X|S^{(i)}) \equiv \mathrm{LSE}\left( \ell(X|S^{(1)}),\ldots, \ell(X|S^{(N_S)})\right)
$$

where $\mathrm{LSE} : \mathbb{R}^n \rightarrow \mathbb{R}$ is called log-sum-exp \cite{blanchard-2019}. An interesting property of $\mathrm{LSE}$ is that it's a smooth approximation to the $\max$ function. This means that for finite sample sizes the monte-carlo estimate of the averaged likelihood will always be too small!

We approximate the mutual information between trajectories as

$$
\mathrm{I}(\mathcal{X}; \mathcal{S}) = \left\langle \ln \frac{\mathrm{P} ( X |  S)}{\left\langle\mathrm P(X | S) \right\rangle_\mathcal{S}} \right\rangle_{\mathcal{X},\mathcal{S}} = \left\langle \ell ( X |  S) - \mathrm{LSE}\left( \ell(X|S^{(1)}),\ldots, \ell(X|S^{(N_S)})\right) + \ln N_S\right\rangle_{\mathcal{X},\mathcal{S}}
$$

which means that for finite amount of signal samples we will _systematically over-estimate_ the mutual information. **TODO: Is that really true? What about $\ln N_S$?** Even worse: the longer the trajectories the bigger the error becomes since the dimensionality of the space of possible signals is growing.

Another way to phrase this insight is that to get a good approximation for the logarithmic average likelihood, our set of signals that we use for monte-carlo sampling should contain many signals that produce a high likelihood. **Therefore it probably is necessary to come up with a scheme to specifically sample signal trajectories for which the likelihood of a particular trajectory is high**. On the other hand the results do not seem to get significantly better when averaging over more trajectories. -->


## Returning to the Simple Model of Gene Expression

To evaluate the estimation procedure described so far in this chapter, we used it to compute the mutual information for the simple biochemical network introduced in @sec:simple_model. As the initial condition we assumed that the system is in the steady state such that we were able to compare our Monte Carlo estimates with the analytical results from @eq:analytical_mi. We generated all signals and responses using a stochastic simulation algorithm (SSA) for the reaction network in @eq:simple_reaction_network. The reaction rates were chosen as specified in @tbl:k and the computed mutual information is always presented in units of _nats_, where $1\ \text{nat} = (1/\ln 2)\, \text{bits}$.

### Simulating Signals and Responses

In the model, the signal dynamics arise from the reactions $\emptyset\xrightarrow{\kappa} S, S\xrightarrow{\lambda}\emptyset$ that describe a simple birth-death process. Therefore, to generate signal trajectories we used the SSA with the initial copy numbers distributed according to the stationary distribution @2009.Gardiner $$P_S(s_0) = e^{-\kappa/\lambda} \frac{(\kappa/\lambda)^{s_0}}{s_0!}$$
for $s_0\in\mathbb{N}_0$. The trajectories generated by the SSA consist of piecewise constant segments in between the transition times $t_0<\ldots<t_{N}$. We can thus describe a trajectory $\mathbf s$ as a function $s: [t_0, t_N]\rightarrow\mathbb{N}_0$ given by
$$
s(t) = \sum^{N-1}_{i=0} s_i\ \theta(t-t_{i}) \theta(t_{i+1}-t) 
$$ {#eq:piecewise_const}
where $s_1,\ldots,s_{N-1}$ are the sequence of copy numbers of $S$ and $\theta$ is the Heaviside function. When we generate response trajectories corresponding to a given signal trajectory, the transition rate $w_t(x+1, x)$ for the reaction $S\xrightarrow{\rho} S+X$ is proportional to the function in @eq:piecewise_const, i.e. 
$$w_t(x+1, x) = \rho\ s(t)\,.$$ {#eq:transition_rate}
@Eq:transition_rate thus precisely specifies how a given signal trajectory interacts with the response. We use [@eq:piecewise_const;@eq:transition_rate;@eq:inverse_function_method] to compute the stochastic waiting times inside the SSA for the response trajectories.

To use the SSA for the generation of response trajectories, we need to solve @eq:inverse_function_method which includes performing a time integral of the transition rates.
The expression in @eq:piecewise_const allows us to integrate a signal trajectory with respect to time to arrive at
$$
\int^t_{t_0} \mathrm d\tau\ s(\tau) = \sum^{N-1}_{i=0} s_i \left[ (t-t_i) \ \theta(t-t_i) - (t-t_{i+1}) \ \theta(t-t_{i+1}) \right]
$$ {#eq:piecewise_linear}
which naturally is a piecewise linear function of $t$. Using @eq:piecewise_linear it is then straightforward to solve @eq:inverse_function_method for $u$ and use the SSA as described in @sec:ssa to generate response trajectories. For response trajectories, we set the initial copy number $x_0$ of $X$ according to the empirically generated distribution $\mathrm P(x_0|s_0)$.

Similarly, for the computation of the likelihood $\mathrm P(\mathbf x|\mathbf s)$ we use  [@eq:transition_probability;@eq:trajectory_probability_product] which also depend on time integrals of transition rates. Hence using the results in this section, we can generate stochastic realizations of signal trajectories and corresponding response trajectories. In other words, we can draw samples from the distributions $\mathrm P(\mathbf s)$ and $\mathrm P(\mathbf x| \mathbf s)$ which we will use to estimate the mutual information.

### Consistent Bias in Comparisons with Analytic Approximations

The Monte Carlo estimate of the mutual information (MI) was performed using two independent computations, a) an average over the _likelihoods_ $\mathrm P(\mathbf x|\mathbf s)$ of random responses with respect to their corresponding signals in @eq:conditional_entropy_estimate that yields the conditional entropy and b) an average over the _marginal probability densities_ $\mathrm P(\mathbf x)$ of sampled responses in @eq:mc_entropy for the marginal entropy. The difference of these two quantities is our estimate of the MI $\mathrm I(\mathcal S, \mathcal X) = \mathrm H(\mathcal X) - \mathrm H(\mathcal X|\mathcal S)$.

It is crucial to note that the two independent computations involved in the MI are not of equal difficulty. The evaluation of $\mathrm P(\mathbf x)$ requires a nested Monte Carlo simulation @eq:mc_marginal for every generated response. The inner compuation of $\mathrm P(\mathbf x)$ can suffer especially high Monte Carlo variance since in @eq:mc_marginal the sampling distribution $\mathrm P(\mathbf s)$ may have only very little overlap with the integrand $\mathrm P(\mathbf x|\mathbf s)$. Specifically, it might be the case that most of the generated signals $\mathbf s_1,\ldots,\mathbf s_{N_S}\sim\mathrm P(\mathcal S)$ contribute very weakly to the integral since their likelihoods $\mathrm P(\mathbf x|\mathbf s_i)\approx 0$ for almost all $i$. Consequently, only a very small number of signals may actually determine the estimate of $\mathrm P(\mathbf x)$ which leads to a high variance of the estimate. In conclusion, we expect the marginal entropy $\mathrm H(\mathcal X)$ to be more difficult to estimate than the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$. This expectation is confirmed by the results in this chapter.

In order to understand how the inner Monte Carlo compuation of $\mathrm P(\mathbf x)$ influences the estimate of the MI we performed an array of simulations with different values of $N_S$. The number $N_S$ specifies how many signals we generate for every nested Monte Carlo computation of $\mathrm P(\mathbf x)$. We expect that increasing $N_S$ leads to a better estimate of the marginal density for every response and therefore improves the estimate of the MI. We compare the simulation results with a reference value for the MI that was obtained by a Gaussian approximation. Using @eq:analytical_mi we can analytically compute the MI $\mathrm I(\mathcal S_T, \mathcal X_T)$ between trajectories of duration $T$ and using @eq:analytical_rate the information rate $I_R$ between $\mathcal S$ and $\mathcal X$. The analytically obtained values for the MI and the information rate are shown as red dotted lines in [@fig:a;@fig:b]. From a single simulation of the MI with a trajectory duration of $T_\text{max}$ we can similarly compute the MI $\mathrm I(\mathcal S_T, \mathcal X_T)$ for any $T\in[0, T_\text{max}]$ by simply truncating all trajectories to duration $T$. An estimate for the information rate can be obtained by approximating the limit in @eq:analytical_rate
$$
\hat I_R = \frac{\mathrm I(\mathcal S_{T_\text{max}}, \mathcal X_{T_\text{max}})}{T_\text{max}}
$$
assuming $T_\text{max} \gg \tau$ where $\tau$ is the longest timescale in the system.

![Estimates of the mutual information $\mathrm I(\mathcal S_T, \mathcal X_T)$ for varying trajectory duration $T$. Every line shows an estimate of $\mathrm I(\mathcal S_T, \mathcal X_T) = \mathrm H(\mathcal X_T) - \mathrm H(\mathcal X_T|\mathcal S_T)$ where the marginal and the conditional entropy were computed in separate simulations. For the marginal entropy we generated $2\times 10^5$ independent responses per simulation and for every response $\mathbf x$ computed the marginal density $\mathrm P(\mathbf x)$ using $N_S$ signal trajectories. Hence, this figure shows how the nested estimate of the marginal density influences the overall result. While increasing $N_S$ helps, all simulations consistently over-estimate the mutual information when compared to the analytic reference result computed using @eq:analytical_mi. The estimation error appears to increase linearly with trajectory length and therefore leads to a biased estimate of the information rate. For the conditional entropy we generated $2\times 10^5$ signals trajectories and for every signal we generated $10^3$ response trajectories for the Monte Carlo average in @eq:conditional_entropy_estimate. ](figures/mutual_information_traj.svg){#fig:a}

In @fig:a we show the the results of the different simulations to compute the MI between trajectories of varying durations $T$. We start by describing the features that all simulations have in common. In all cases we see a linear increase of the mutual information with trajectory length. The linearity is to be expected for trajectory durations $T\gg\tau$ where $\tau$ is the correlation time of the system @2010.Tostevin. If the signal is much longer than the correlation time of the system it is impossible for any future response to still gain information about the beginning of the signal. Hence, for short trajectories, the increase in MI per unit time may be smaller than the asymptotic information rate since the signal has not yet reached the maximal integration length of the system. Thus, while there exist situations where for small $T$ the MI increases faster than linear, analytic computations using the Gaussian approximation confirm a fully linear increase of the MI from the start for the simple reaction network under consideration. In @fig:a, we show the result of @eq:analytical_mi for the Gaussian approximation as the dotted "reference" line.

For trajectory length $T = 0$ @fig:a shows that all estimates—including the analytical result—coincide at a small non-zero value for their mutual information. This value is the mutual information between the initial conditions of the signal and the response. The agreement of all simulations shows that using the approach detailed in @sec:initial_condition we were able correctly estimate the conditional distribution $\mathrm P(x_0|s_0)$. For every response needed during our simulations we first generated 1000 signals backwards in time for $T_\text{past} = 500$ to choose a correct starting point for the response and to correctly compute the likelihood $\mathrm P(x_0|s_0)$ of the initial condition.

![Estimated information rate $I_R=\mathrm I(\mathcal S_T, \mathcal X_T)/T$ from the same data that was used in @fig:a. The reference line represents the analytical result $I_R = 0.00183$ from @eq:analytical_rate. For a correct Monte Carlo simulation we expect the estimated information rate to converge towards the reference line as $T\rightarrow\infty$. However, while we see that the estimates converge towards a stable value for large $T$, all simulations considerably over-estimate the information rate. Even though for increasing $N_S$ the estimates _do_ improve we can not acquire reasonable results for the simulations shown.](figures/information_rate_traj.svg){#fig:b}

Thus, we see from @fig:a that the estimated MI grows linearly with $T$ and that all estimates agree on the MI for $T=0$. However, we see that each of the MI curves has a different slope. All Monte Carlo simulations significantly over-estimate the slope of the analytical reference line. As $N_S$ increases we see a slight decrease of the MI slope, however even for largest simulation with $N_S=1000$ the slope differs by more than 50% from the analytical reference. Further increasing $N_S$ seems to yield diminishing returns. Hence, from the simulations we can't seem to accurately compute the information rate since it strongly depends on estimating the MI slopes correctly.

Asymptotically for $T\rightarrow\infty$ the slope of the MI with respect to the trajectory duration $T$ converges towards the information rate $I_R$. In @fig:b we see that all simulations show the convergence of the estimated information rate $I_R\approx \mathrm I(\mathcal S_T, \mathcal X_T) / T$ to a stable value already at $T\approx 150$. We again find that we over-estimate the information rate for simulations with $N_S\in\{10, 20, 50, 100, 200, 500, 1000\}$ even though we see slight improvement for larger $N_S$. We propose that even further increasing $N_S$ and thus spending more CPU time per simulation is not worthwhile to improve the estimates.

![Estimated information rate $I_R=\mathrm I(\mathcal S_T, \mathcal X_T)/T$ for three estimate with very large $N_S$. The reference line represents the analytical result $I_R = 0.00183$ from @eq:analytical_rate. This plot highlights that even strong increases of $N_S$ beyond $10^3$ do not yield accurate results for the information rate. Our conclusion from this is that the simple Monte Carlo approach is too naïve and we have to find another way to improve accuracy, instead of spending more CPU time.](figures/information_rate_traj_long.svg){#fig:c}

To support that claim, in @fig:c we show the results of two additional simulation with very high $N_S$ and compare them to the analytical reference and the MI for $N_S=1000$. Even for these very resource-intensive simulations we don't get better approximations for the information rate than for $N_S=1000$. From these results we conclude that our approach is plagued by some more fundamental problem that needs to be solved with other means than merely using more computational resources.

While we expect the main cuprit to be the nested Monte Carlo estimates of the marginal density $\mathrm P(\mathbf x)$, it is difficult to obtain convincing evidence for that claim using the results shown so far. For an indivdual response trajectory $\mathbf x$ we don't have a way to judge the accuracy of our estimate of $\mathrm P(\mathbf x)$. The only analytical result that we can use to check our results is the MI from the Gaussian approximation which is shown as the reference line in [@fig:a;@fig:b;@fig:c]. Therefore we decided to use the same Monte Carlo procedure to compute the MI with a simpler stochastic description of the signal and response trajectories which allows us to analytically validate all intermediate results. Thus, instead of using the SSA to fully simulate the signal and response dynamics we make the Gaussian approximation as introduced in @sec:simple_model. Using that approximation, generating a signal or response amounts to taking a draw from a multivariate normal distributions with known covariance matrix.

## Gaussian Approximation

Clearly, there is an issue with the computation of the mutual information (MI) for trajectories that prevents the accurate estimation using the simple Monte Carlo strategy shown above. In [@fig:a;@fig:b;@fig:c] we saw that we consistently over-estimate the information rate regardless of the amount of computing time that was invested. To allow a deeper analysis of the underlying issue we repeated the estimates of the MI using a simpler, Gaussian model for signals and responses. Instead of simulating trajectories exactly in continuous time using the SSA, we choose a small time-interval $\Delta t$ which we use to discretize the system. Conversely, while the SSA simulates individual reactions that lead to a discrete jumps in the trajectories, in the Gaussian approximation the copy numbers are approximated by continuous variables. Hence, we switch from a _discrete states and continuous time_ description to a _continuous states and discrete time_ approximation of the biochemical network. 


### Time Discretization

We discretize the time axis using a set of equidistant points $t_n = n\Delta t$ for $n=1,\ldots,d$. A signal trajectory is described a $d$-dimensional vector $\mathbf s = (s(t_1), \ldots, s(t_d))^T$ where $s(t_i)$ describes the deviation of the signal from its mean value $\bar s$ at time $t_i$. The null vector therefore represents a trajectory that remains constant at $\bar s$ for its entire duration whereas a non-zero vector indicates fluctuations. The Gaussian approximation for a biochemical network uses the correlation functions to derive a multivariate normal distribution such that the discretized system dynamics approximate the true stochastic dynamics of the network. Given a time-discretization we can use @eq:covariance_from_corr compute the covariance matrices $C_{ss}, C_{sx}, C_{xs}, C_{xx}$ from the correlation functions in @eq:correlation_functions. From these covariance matrices we can derive all stochastic properties of the Gaussian approximation.

The signal dynamics are given by the probability distribution
$$
\mathrm P(\mathbf{s}) = \frac{1}{\sqrt{\left( 2\pi  \right)^{d} \det C_{ss}}} \;\exp\left[-\frac12\ \mathbf s^T\ C_{ss}^{-1}\ \mathbf s\right]\,,
$$
and similarly the response dynamics are determined by the marginal distribution
$$
\mathrm P(\mathbf{x}) = \frac{1}{\sqrt{\left( 2\pi  \right)^{d} \det C_{xx}}} \;\exp\left[-\frac12\ \mathbf x^T\ C_{xx}^{-1}\ \mathbf x\right]\,.
$$ {#eq:analytical_marginal}
Notably, @eq:analytical_marginal allows us to analytically compute the marginal density $\mathrm P(\mathbf x)$ without performing a Monte Carlo estimate of the integral $\int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\,\mathrm P(\mathbf x|\mathbf s)$.
The joint distribution of $\mathcal S$ and $\mathcal X$ is given in @eq:joint_multivariate which is a $2d$-dimensional multivariate normal density function with covariance matrix
$$
Z =  \begin{pmatrix}
C_{ss} & C_{xs} \\
C_{sx} & C_{xx}
\end{pmatrix}\,.
$$ {#eq:corr_z2}

![Matrix plots of the full correlation matrix $Z$ from @eq:corr_z2 for different values of dimensionality $d$ and $\Delta t$. Brighter colors indicate higher matrix element values. We can clearly observe the block structure of $Z$ in every matrix plot. For every matrix plot, the element with coordinates $(i,j)$ in the top left quadrant shows the correlations $\langle s(i\Delta t) s(j \Delta t)\rangle$. In the top right quadrant we see the correlations $\langle s(i\Delta t) x(j \Delta t)\rangle$ and in the lower quadrants we see $\langle x(i\Delta t) s(j \Delta t)\rangle$ and $\langle x(i\Delta t) x(j \Delta t)\rangle$ on the left and right side respectively. The product $T = d\Delta t$ is the duration of the signal and response trajectories. The quantity $d\Delta t$ also serves as a rough measure of the sparsity of the correlation matrices (i.e. the fraction of matrix elements lower than some cutoff). Along diagonals from bottom left to top right, the product $T = d \Delta t$ is constant. Indeed, we see that along these diagonals the sparsity is roughly constant. Yet, as we move along such a diagonal of equal sparsity $d\Delta t$ but increasing dimensionality $d$ we see correlation matrices that display the same features in a gradually more refined and smooth way.](matrix_plots.png){#fig:corr}

Using this discretization we have two parameters left to tune. We can freely choose the number $d$ and offsets $\Delta t$ of our time samples. The duration of the trajectories $\mathbf s$ and $\mathbf x$ is given by the product $T=d\Delta t$. In @fig:corr we show matrix plots of the joint covariance matrix $Z$ for different values of $d$ and $\Delta t$. We can also observe that $d$ determines the dimensionality of the problem while the product $d \Delta t$ serves as a measure for the sparsity of the correlation matrices. 

Note that the choice of $\Delta t$ affects how well the discretized trajectories approximate physical continuous-time trajectories. Thus, usually we would have to choose $\Delta t$ to be relatively small such that the discretized process has similar dynamics compared to the continuous Gaussian process. For our purposes, this consideration is of secondary importance. We are merely comparing the Monte-Carlo estimates for the mutual information to the analytic results, obtained by computing the determinants of the covariance matrices. For a consistent Monte-Carlo procedure we should expect convergence of the estimates to the analytic results given enough computing time. That is, while the choices of $d$ and $\Delta t$ should not affect the correctness of our procedure, it very well may strongly affect the computational resources needed for efficient estimation.

### Direct Importance Sampling

We want to use this fully Gaussian model to understand how the sample sizes of the different Monte Carlo steps affect the estimate and whether there exists a bias in the approximation. We calculate the marginal entropy as a Monte Carlo average over the logarithms of the marginal distribution densities of $N_x$ sampled responses as shown in @eq:mc_entropy. The evaluation of the marginal density itself requires a Monte Carlo average over $N_s$ sampled signals (@eq:mc_marginal). Hence to evaluate the marginal density we need to perform nested averaging as shown in @eq:mc_entropy_notation. We performed this procedure for various values of $N_s$ and $N_x$ and compared the estimate with reference results using the analytical expression for the entropy of a multivariate Gaussian distribution.

Both, increase of $N_x$ and increase of $N_s$ should lead to an improved estimate of $\mathrm H(\mathcal X)$. To understand the accuracy of an estimate with a given $N_s$ and $N_x$ we repeat the estimation procedure multiple times and compute the mean and the standard deviation of the individual estimation results.

![Top: relative error for the marginal entropy as a function of $1/N_x$. Bottom: empirical variance of ensembles of 144 estimates. The solid lines show a linear extrapolation of the data points for $N_x \rightarrow\infty$. All estimates were performed using a constant number of signal samples $N_s = 400$ and for $d = 200$. The linear extrapolation in the bottom plot indicates that we do predict the variance of the results to vanish in the limit of infinite sampling. This behavior is generally expected for Monte Carlo estimates. Strikingly however, we find that there is a consistent offset of the average estimate from the correct result, even in the limit $N_x \rightarrow\infty$. We see that the bias scales with the sparsity of the covariance matrices. The relative error is computed using $\mathrm H_\text{estimate}/\mathrm H_\text{analytical} - 1$ where $\mathrm H_\text{estimate}=\sum^{M}_{i=1} \hat{\mathrm H}_i/M$ is the average over the results of the $M=144$ marginal entropy estimates that were performed using @eq:mc_entropy_notation and $\mathrm H_\text{analytical}$ is the value resulting from an analytical computation of the marginal entropy. The empirical variance shown is $\sum^{M}_{i=1} (\hat{\mathrm H}_i - \mathrm H_\text{estimate})^2/M$.](relative_error_responses.svg){#fig:rel_err_responses}

In @fig:rel_err_responses we see how the relative error of our estimate varies with the number of simulated responses $N_x$. Here use the same number of signals per response $N_s$ for all estimates. While---as expected---the variance of the estimate decreases when we increase $N_x$ we find that especially for very sparse covariance matrices we consistently over-estimate the marginal entropy. Indeed, we find that the systematic bias in our results seems to be independent of $N_x$.

We found that the bias is stronger for correlation matrices with higher sparsities $d\Delta t$. Since the sparsity grows with trajectory duration we can expect an increasingly strong over-estimation for longer trajectories. The sparsity can be increased either by decreasing the time-resolution or by increasing the dimensions of the covariance matrix. To understand how these parameters relate to each other we tested how the estimation error changes when we increase the dimensionality of the correlation matrices while the sparsity remains constant.

@Fig:sparsity shows how large the estimation error for the marginal entropy $\mathrm H(\mathcal X)$ is on average for different levels of sparsity. We see that in all cases that increasing the sparsity leads to larger errors in the estimates. Additionally we find that for a given sparsity value, the estimates with high-dimensional covariance matrices are slightly worse. As we keep increasing the number of dimensions $d$ at constant sparsity $d\Delta t$, thus decreasing $\Delta t$, the matrices gradually become a more faithful representation of the continuous correlation functions of the system (see @fig:corr). Extrapolating the lines in @fig:sparsity we project that for very large covariance matrices, the sparsity is the only determining factor of the estimation bias.

![Absolute error of marginal entropy estimates for different values of the sparsity $d\Delta t$ of the correlation matrices. We see that for high dimensionality the lines of constant sparsity become increasingly flat. This indicates that for high-dimensional systems the sparsity of the covariance matrix is a good measure for the difficulty of correct estimation. We therefore claim that the bias of the entropy estimate for the Gaussian system primarily depends on the sparsity of the covariance matrix. Note that for lower numbers of dimensions the covariance matrices of along the diagonals of equal sparsity look more blocky (see @fig:corr). That may be an indicator why the estimation error is not constant for a given sparsity at lower dimensions.](sparsity.svg){#fig:sparsity}

![Relative error $\mathrm H_\text{estimate}/\mathrm H_\text{analytical} - 1$ as a function of $1/N_s$. We can see that the relative error in the marginal entropy estimate increases with the sparsity $d\Delta t$ (i.e. with trajectory duration). The linear extrapolating lines emphasize that there is a noticeable but very slight decrease in error as $N_s\to\infty$. This seems puzzling since for infinite sampling we should expect the error to vanish. Apparently for high-sparsity covariance matrices we need extraordinarily many signal samples to achieve unbiased estimates.](error_grid.svg){#fig:error_regression}

As a next step we investigated how changes in the sampling for the marginal density $\mathrm P(\mathbf x_i)$ affect the estimation bias. Thus in @fig:error_regression we show how the increase of simulated signals per response $N_s$ improves the estimate of the marginal entropy. Here we again see that for high number of dimensions we over-estimate the marginal entropy. An increase of $N_s$ does lead to slightly less over-estimation but the linear extrapolation indicates that even if we choose enormously high values for $N_s$ we can not expect to reduce the bias substantially.

For a given number of samples, the fraction of the trajectory space probed by the Monte Carlo scheme is lower for longer durations. Therefore, for a given number of Monte Carlo samples we expect the estimate to become worse for longer trajectories, i.e. when the sparsity of the covariance matrix is high. This is confirmed by our results. Furthermore and more surprisingly we find that we consistently over-estimate the marginal entropy and while increasing $N_s$ _does_ reduce the bias slightly it appears to require an astronomically high sampling in signal trajectory space to reach arbitrary low errors. An increase in $N_x$ however reduces the variance of the results but does not influence the bias at all.

Thus we are lead to believe that the main difficulty in estimating the marginal entropy is the Monte-Carlo marginalization of the probability density function. To estimate $\mathrm P(\mathbf x)$ we sample signals from the marginal distribution $\mathrm P(\mathbf s)$ and average over the likelihoods $\mathrm P(\mathbf x | \mathbf s)$. However as the space of signals becomes increasingly vast for longer trajectories, it becomes more and more unlikely to sample a signal $\mathbf s$ where $\mathbf x$ has a non-vanishing likelihood of occurring. Hence the duration of the trajectories strongly influences the bias of the marginal entropy computation and we are well advised to keep the trajectories as short as possible. To capture the essential system dynamics however, the trajectories must be at least as long as the longest timescale $\tau$ in the system. If we are interested in the marginal entropy for timescales longer than the longest timescale in the system we can then use the fact that trajectory pieces of duration $> \tau$ become independent of each other and we can add the individual entropies of these pieces to get the full entropy for such a long trajectory. For example, if we were interested in the marginal entropy of trajectories of duration $T=\ell\tau$ where $\ell\in\mathbb N^+$ we could approximate it by estimating the entropy $\hat{H}_\tau$ of trajectories of duration $\tau$ and then take $\ell \hat{H}_\tau$ as an approximation for $\hat{H}_{\ell\tau}$. While performing the estimate $\hat{H}_\tau$ is computationally much cheaper than directly estimating $\hat{H}_{\ell\tau}$ and therefore for a given amount of CPU-time we get a lower systematic bias of the estimate $\hat{H}_\tau$, the bias of $\ell \hat{H}_\tau$ of course still scales linearly with $\ell$. Therefore it is still in our best interest to reduce the systematic bias of the marginal entropy as much as possible.

### Umbrella Sampling {#sec:umbrella}

To get a better estimate of $\mathrm P(\mathbf x)$ we decided to use _umbrella sampling_, i.e. to bias our sampling strategy towards signals that we expect to have a high likelihood for the given response. Since Monte-Carlo estimates depend on the distribution of the chosen samples we must correct our estimate by re-weighing the samples accordingly.

More specifically, given a sampling distribution $w(\mathbf s)$ (which is normalized like a probability density function) we can write
$$
\mathrm P(\mathbf x) = \int \mathrm d\mathbf s\ w(\mathbf s)\frac{\mathrm P(\mathbf s)\mathrm P(\mathbf x | \mathbf s)}{w(\mathbf s)} = \left\langle \frac{\mathrm P(\mathbf s)\mathrm P(\mathbf x | \mathbf s)}{w(\mathbf s)} \right\rangle_{w(\mathbf s)}
$$
which shows us how to compute $\mathrm P(\mathbf x_i)$ from signal trajectories $\mathbf s_1^w, \ldots, \mathbf s_{N_s}^w$ which are distributed according to the sampling distribution given by $w$:
$$
\mathrm P(\mathbf x_i) \approx \frac1{N_s} \sum\limits_{j=1}^{N_s} \frac{\mathrm P(\mathbf s_j^w)\mathrm P(\mathbf x_i | \mathbf s_j^w)}{w(\mathbf s_j^w)} \equiv P_{\mathbf x_i, N_s}^w \,.
$$

The choice of the sampling distribution $w$ has a direct impact on the variance of an ensemble of estimates $P_{\mathbf x_i, N_s}^w$. Indeed, for any given $\mathbf x_i$ there is an optimal choice for the sampling distribution $w_\text{opt}$ such that the variance of the estimates vanishes. This optimal choice is given by $w_\text{opt}(\mathbf s) = \mathrm P(\mathbf s | \mathbf x_i)$ which is easily confirmed by the calculation
$$
P_{\mathbf x_i, N_s}^{w_\text{opt}} = \frac1{N_s}\sum\limits_{j=1}^{N_s} \frac{\mathrm P(\mathbf s_j^w)\mathrm P(\mathbf x_i | \mathbf s_j^w)}{\mathrm P(\mathbf s_j^w | \mathbf x_i)} = \frac1{N_s}\sum\limits_{j=1}^{N_s} \mathrm P(\mathbf x_i)
$$ {#eq:opt_sampling}
where in the last step we applied Bayes' rule. Since the expression above is completely independent of the chosen signal samples the result is deterministic and thus has zero variance. We also see from @eq:opt_sampling that in practice we can't directly use $w_\text{opt}(\mathbf s) = \mathrm P(\mathbf s | \mathbf x_i)$ as our sampling distribution since the evaluation of $\mathrm P(\mathbf s | \mathbf x_i) = \frac{\mathrm P(\mathbf x_i|\mathbf s) \mathrm P(\mathbf s)}{\mathrm P(\mathbf x_i)}$ itself depends on $\mathrm P(\mathbf x_i)$ which is precisely the quantity we are interested in estimating.

![Relative error as a function of the dimensionality $d$. The solid lines show the results using non-optimized sampling while the dashed lines show the results when using a sampling distribution close to the optimal distribution $\mathrm P(\mathbf s|\mathbf x)$. We see that with optimized sampling there is no consistent over-estimation anymore. All estimated were done using $d = 200$ dimensional covariance matrices.](sampling2.svg){#fig:rel_err_opt}

Instead, we can try to obtain a sampling distribution that is as close as possible to $w_\text{opt}(\mathbf s)$. A known approach involves using random samples from $\mathrm P(\mathbf s | \mathbf x_i)$ to pick the most optimal sampling distribution from a family of candidate distributions [@2011.Chan]. Generating so-called _posterior samples_ from $\mathrm P(\mathbf s | \mathbf x_i) \sim \mathrm P(\mathbf x_i|\mathbf s) \mathrm P(\mathbf s)$ is generally possible without knowledge of the normalization factor $\mathrm P(\mathbf x_i)$ e.g. by using Metropolis-Sampling [@1991.Müller;@1994.Tierney]. To test within the Gaussian framework whether such an approach to importance sampling could work in principle, we generate 400 posterior samples by directly sampling from the analytically known posterior distribution $\mathrm P(\mathbf s | \mathbf x_i)$. We compute the empirical mean $\bar{\mathbf s}$ and the empirical covariance $\bar C_{\mathbf s|\mathbf x_i}$ of these samples as parameter estimates for a multivariate Gaussian $\mathcal N(\bar{\mathbf s}, \bar C_{\mathbf s|\mathbf x_i})$ and use the latter as an optimized sampling distribution.

In @fig:rel_err_opt we show that using optimized sampling we can strongly reduce the systematic bias in marginal entropy estimation. As expected, importance sampling is especially useful when the sparsity is very high, i.e. the trajectories are long. It is clear that for longer trajectories we expect $\mathrm P(\mathbf s | \mathbf x_i)$ to be a much more narrow sampling distribution than $\mathrm P(\mathbf s)$ whenever the $\mathcal S$ and $\mathcal X$ are not completely independent. Consequently, it becomes more and more unlikely to obtain a sample $s^\prime$ from $\mathrm P(\mathbf s)$ such that $\mathrm P(s^\prime | \mathbf x_i) > \epsilon$ for any $\epsilon > 0$ and therefore more difficult to accurately estimate $\mathrm P(\mathbf x_i)$ using an unbiased sampling distribution.

<!-- To choose a sensible sampling distribution we use the fact that we generate the responses $\mathbf x_i$ by first picking a signal $\mathbf s_i^\star$ and subsequently picking from the likelihood distribution $\mathrm P(\mathbf x|\mathbf s_i^\star)$. As a matter of fact, $\mathbf s_i^\star$ can be regarded as a single random sample from $\mathrm P(\mathbf s | \mathbf x_i)$, i.e. the “optimal” sampling distribution. It is not unreasonable to assume in this context that a multivariate Gaussian distribution centered at $\mathbf s_i^\star$ could turn out to be a good weighing function. -->

### Estimating the Conditional Entropy


![Comparison of the relative error of conditional entropy estimates versus marginal error estimates. The relative errors are shown on a logarithmic scale as a function of the sparsity. We can see that the relative error for the estimate of the conditional entropy is a few orders of magnitude smaller than the estimates of the marginal entropy. All estimates were performed with $N_x=25600$ and $N_s=1000$.](conditional.svg){#fig:conditional}

For both, marginal entropy and conditional entropy we have to evaluate the likelihood $\mathrm P(\mathbf x| \mathbf s)$ a total of $N_s N_x$ times. To compare the accuracy of we performed estimates of the marginal entropy with and without optimized sampling together with estimates of the conditional entropy for $N_s = 1000$ and $N_x = 25600$. In @fig:conditional we show the relative error of both, marginal and conditional entropy estimates as a function of the sparsity. We find that the estimate of the conditional entropy is very accurate regardless of sampling size. Even with optimized sampling the marginal entropy estimate is roughly two orders of magnitude worse than a comparable conditional entropy estimate.

## Discussion

Since both the estimate of the marginal entropy and of the conditional entropy require the computation of two nested Monte-Carlo averages one could expect the results of both estimates to be of similar accuracy. Yet we find that computing the marginal entropy is much more challenging than computing the conditional entropy. While analyzing the estimation procedure for the marginal entropy we found the main source of error to arise from the computation of the marginal probability density $\mathrm P(\mathbf x)$. While the computation of this density suffers from high Monte-Carlo variance when we are not carefully optimizing our trajectory sampling procedure, the real issue arises in the next step when we have to _compute the logarithm of the marginal probability density_ to estimate the marginal entropy $\mathrm H(\mathcal X) \approx \sum^{N_x}_{i=1}-\ln\mathrm P(\mathbf x_i) / N_x$ from our sampled response trajectories (see also @eq:mc_entropy_notation). For us, computing $\ln\mathrm P(\mathbf x_i)$ means computing the logarithm of an average. Taking the logarithm---which is concave function---of a Monte-Carlo average leads to a consistent bias. We see the existence of this bias in our results since we consistently over-estimate the marginal entropy. Indeed, the reason why it is much easier to get a good estimate of the conditional entropy is that in the latter case we have to average the logarithm of a quantity (see @eq:conditional_entropy) rather than taking the logarithm of an average, as we need for the marginal entropy.

We can thus conclude that the main difficulty of obtaining a good estimate for the mutual information between trajectories lies in the efficient and accurate computation of the marginal entropy. A viable approach for this seems to be to use importance sampling in the signal space in the computation of the marginal probability density. Our results indicate that such an approach could also work for trajectories generated using a fully stochastic model of a biochemical network. Another direction to pursue might be to use the replica trick, based on the mathematical identity that $\ln Z = \lim_{n\to 0} (Z^n-1) / n$. This may allow us to eliminate the systematic bias by circumventing the need to take the logarithm of an estimate.
