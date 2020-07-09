---
bibliography: ["library.bib"]
autoEqnLabels: true
link-citations: true
linkReferences: true
cref: true
---

# Monte-Carlo Estimate of the Mutual Information

Equipped with analytical formulae for the computation of trajectory probabilities and with methods for efficient stochastic simulation, we can start to develop estimates for the mutual information. The basis for our method is @eq:mi_form2; specifically we separate the mutual information into two parts that are computed independently, the _marginal entropy_ $\mathrm H(\mathcal X)$ and the _conditional entropy_ $\mathrm H(\mathcal X|\mathcal S)$. Since signals and responses are time-varying, both entropies involve integrals over spaces of trajectories which are high-dimensional. At high dimensionality, direct numerical integration is not viable and instead, we use Monte-Carlo approaches based on random sampling of signals and responses. 

While Monte-Carlo methods comprise a wide variety of approaches to stochastically evaluate integrals or sums the common idea is easily stated. We have a state space $U$ and a probability distribution $p_U$ on that state space. The problem is to numerically evaluate integrals of the form
$$
F = \langle f(u) \rangle = \int\limits_{U} \mathrm du\ p_U(u)\; f(u)
$$
where $f: U\rightarrow\mathbb R$ is some function of interest. If $U$ is high-dimensional it is very time-consuming to estimate it by direct numerical integration. Instead, we generate random samples $u_1,u_2,\ldots$ from the probability distribution $p_U$ such that by the _law of large numbers_ we have the equality
$$
F = \lim_{N\rightarrow\infty} \frac{\sum^N_{i=1} f(u_i)}{N}
$$ {#eq:inf_mc_estimate}
i.e. the sample average of random samples converges towards their mean. In practice we can't evaluate the limit in @eq:inf_mc_estimate and we approximate $F$ using a finite number $N$ of random samples
$$
\hat{F} = \frac{\sum^N_{i=1} f(u_i)}{N}\,.
$$
The variance of Monte-Carlo estimates typically decreases much faster than TODO. Monte-Carlo methods depend on the efficient and correct sampling of the corresponding probability distribution. In this thesis ...

In this chapter we will show how to use Monte-Carlo integration to compute the mutual information between stochastic trajectories. The computation requires the estimation of two quantities, marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$. We will show how to perform Monte-Carlo estimates for both entropies, however as we will describe the marginal entropy turns out to be computationally much more difficult to estimate. We conclude the chapter with a thorough analysis of the challenges that arise by ...

## Monte-Carlo Estimate for the Marginal Entropy

We start by considering an abstract system, consisting of two random variables $\mathcal S$ and $\mathcal X$ representing the possible signals and responses, respectively. We'll assume to have efficient computational methods to generate random samples from either of these random variables. In this subsection we show how to use Monte-Carlo techniques to estimate the marginal entropy $\mathrm H(\mathcal X)$ while in the following subsection we describe a similar method to estimate the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$. 

We intend to compute the marginal entropy $\mathrm H(\mathcal X)$ using Monte-Carlo (MC) sampling to evaluate the necessary integrals. First, given number of random samples $(\mathbf x_i)_{i=1,\ldots,N_x}$ taken from $\mathcal X$ we propose as an estimate for the entropy
$$
\mathrm H(\mathcal X) = -\int\mathrm d\mathbf x\ \mathrm P(\mathbf x)\ln \mathrm P(\mathbf x) \approx - \frac{\sum\limits_{i=1}^{N_x} \ln\mathrm P(\mathbf x_i)}{N_x} \,.
$$ {#eq:mc_entropy}
I.e. using a large number of response samples, we compute the sample average of their individual log-probabilities $\ln\mathrm P(\mathbf x_i)$. As explained in the following section, from a biochemical of a cell we are not able to directly compute $\mathrm P(\mathbf x_i)$. We do however have an efficient way of _generating_ the responses $\mathbf x_1,\mathbf x_2,\ldots$ using the stochastic simulation algorithm.

Nonetheless we see from @eq:mc_entropy that we _do_ have to evaluate $\mathrm P(\mathbf x_i)$ for every generated sample. To solve this problem we choose to evaluate $\mathrm P(\mathbf x_i)$ by doing another nested Monte-Carlo integration using signal samples $(\mathbf s_j)_{j=1,\ldots,N_s}$ from $\mathcal S$ to write
$$
\mathrm P(\mathbf x_i) = \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\ \mathrm P(\mathbf x_i|\mathbf s) \approx \frac{\sum\limits_{j=1}^{N_s} \mathrm P(\mathbf x_i | \mathbf s_j)}{N_s} \,.
$$ {#eq:mc_marginal}
Hence, for every response sample $\mathbf x_i$ we have to perform a sample average over the conditional probabilities $\mathrm P(\mathbf x_i|\mathbf s_j)$. Again, anticipating results from the following section, we make use of the fact that a stochastic model for a biochemical network allows us to efficiently compute $\mathrm P(\mathbf x_i|\mathbf s_j)$.

While for a low-dimensional signal space it is feasible to instead compute the marginalization integral @eq:mc_marginal using direct evaluation [@2019:Cepeda-Humerez] we choose to use a nested MC simulation to also be able to handle high-dimensional signal spaces. This is crucial since time-varying signals are described by high-dimensional trajectories.

We can summarize the estimation procedure for the marginal entropy using the equation
$$
\mathrm H(\mathcal X) = -\left\langle \ln \left\langle \mathrm P(\mathbf x | \mathbf s) \right\rangle_{\mathrm P(\mathbf s)} \right\rangle_{\mathrm P(\mathbf x)}
$$ {#eq:mc_entropy_notation}
where we use the notation $\langle f(x)\rangle_{g(x)}$ for the expected value of $f(x)$ when $x$ is distributed according to the probability density given by $g(x)$. Thus when thinking in mathematical terms we have the shorthand $\langle f(x)\rangle_{g(x)} \equiv\int \mathrm dx\ g(x) f(x)$. We can also easily translate this notation into a Monte-Carlo estimate, i.e. $\langle f(x)\rangle_{g(x)} = \lim\limits_{N\rightarrow\infty}\frac{\sum_{i=1}^N f(x_i)}{N}$ where $x_1, x_2,\ldots$ are independent samples of the probability distribution given by $g(x)$.

## Estimating the Conditional Entropy

We can also estimate the _conditional entropy_ using MC averages over trajectories. We express the conditional entropy using the notation introduced above
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

### Generating Responses for Time-Varying Signals

- Signals are generated from a known stochastic process

Both, the estimate for the marginal entropy $\mathrm H(\mathcal X)$ and the estimate for the conditional entropy $\mathrm H(\mathcal X|\mathcal S)$ require the generation of appropriate responses to a given time-varying signal. To find stochastically correct responses for a given signal we have to understand how the biochemical network _couples_ to the signal. In many cases the signal itself is a molecular species that participates in the reaction network. In this scenario, an abundance of signal directly leads to an increased rate of all reactions in which the signal molecule acts as a reactant. In other cases the signal may be a thermodynamic quantity as for example the environmental temperature. Changes in temperature usually effect the rates of _all_ reactions in a chemical reaction network. Thus, in general we think of a signal as a quantity whose value affects some or all of the reaction rates of the biochemical network at a given instant. Since we model the signal as a time-varying process, for a given realization of the signal we can describe the corresponding responses using a master equation with time-dependent reaction rates.

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

In summary, as long as we can formulate a stochastic process that describes the signal dynamics _and_ we understand how the signal level affects the reaction rates of the biochemical network, we can in principle use the Monte-Carlo procedure described in the previous section to estimate the mutual information. From the signal's SDE we generate stochastic realizations using appropriate numerical integrators for SDEs such as the _Euler–Maruyama method_ or the _stochastic simulation algorithm_. 

## Practical Concerns

### Computing the probability distribution of the initial state

We first look at the term $P_0 = \mathrm P(x_0,t_0 | S)$. Since $S$ is a trajectory in time we can directly conclude from causality that

$$
\mathrm P(x_0, t_0 | S) = \mathrm P(x_0, t_0 | S_{t \leq t_0})
$$

where $S_{t \leq t_0}$ is the temporal piece of the signal up to $t_0$. We further suppose that the signal itself is markovian and therefore has no memory of its past. With this simplification we get

$$
\mathrm P(x_0, t_0 | S) = \mathrm P(x_0, t_0 | S_{t = t_0}) = \frac{\mathrm P((x_0, t_0), (s_0, t_0))}{\mathrm P(s_0, t_0)} \,.
$$

We estimate $P_0$ using gaussian kernel density estimation to approximate both, the joint distribution of $X_0, S_0$ and the marginal distribution of $S_0$.

Knowing the probabilities of the initial condition of both response and signal we can directly estimate the mutual information of $\mathcal{X}_{t=t_0}$ and $\mathcal{S}_{t=t_0}$:

$$
\mathrm I(\mathcal{X}_{t=t_0}, \mathcal{S}_{t=t_0}) = \int ds_0\int dx_0\; \mathrm{P}(x_0, s_0)\; \ln \frac{\mathrm{P} (x_0, s_0)}{\mathrm{P} (x_0) \mathrm{P} (s_0)}
$$

### Computing Probability Densities

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

### Computing the likelihood

The Probability density of a markovian trajectory can be expressed as

$$
\mathrm P(X) = \mathrm P(x_0,t_0;x_1,t_1;\ldots;x_{N-1},t_{N-1}) = \mathrm P(x_0,t_0 ) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}) \,.
$$

Therefore the problem of calculating the likelihood for a particular trajectory amounts to solving two independent problems:

1. estimating the probability density of the starting point $\mathrm P (x_0, t_0)$ of a response
2. calculating the transition probabilities $\mathrm P(x_n,t_n|x_{n-1},t_{n-1})$

For a given chemical reaction network we can write down the chemical master equation. The chemical master equation contains all the information needed to compute the individual terms $\mathrm P(x_n,t_n|x_{n-1},t_{n-1})$ for the entire system.

To calculate the mutual information between $\mathcal{S}$ and $\mathcal X$ we have to consider the entire reaction network containing the components both in $S$ and in $X$. The precise reaction dynamics of the response part of the chemical network crucially depend on the observed signal trajectory. Therefore the chemical master equation for the whole reaction network allows us to compute the likelihood of a response trajectory for a particular signal trajectory:

$$
\mathrm P(\mathcal X = X|\mathcal S = S) = \mathrm P(x_0,t_0;x_1,t_1;\ldots;x_{N-1},t_{N-1} | S) = \mathrm P(x_0,t_0 | S) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S) \,.
$$


### The probability density for the starting point of a trajectory



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

Another way to phrase this insight is that to get a good approximation for the logarithmic average likelihood, our set of signals that we use for monte-carlo sampling should contain many signals that produce a high likelihood. **Therefore it probably is necessary to come up with a scheme to specifically sample signal trajectories for which the likelihood of a particular trajectory is high**. On the other hand the results do not seem to get significantly better when averaging over more trajectories.


## Computation of the Mutual Information for a Simple Model of Gene Expression

- Show the reaction network
- Show example trajectories

### Consistent Bias in Comparisons with Analytic Approximations



<!-- ## Simulating a Biochemical Network Driven by an External Signal



## Monte-Carlo computation of the Mutual Information

While so-called Monte-Carlo methods comprise a wide variety of approaches to stochastically evaluate integrals or sums the common idea is easily stated. We have a state space $U$ and a probability distribution $p_U$ over that state space. The problem is to evaluate

$$
\langle f(u) \rangle \equiv \int\limits_{u \in U} \mathrm du\; f(u) p_U(u)
$$

where $f: U\rightarrow\mathbb R$ is some smooth function. If $U$ is high-dimensional it is very time-consuming to estimate it by direct numerical integration.

### Computing the likelihood

The Probability density of a markovian trajectory can be expressed as

$$
\mathrm P(X) = \mathrm P(x_0,t_0;x_1,t_1;\ldots;x_{N-1},t_{N-1}) = \mathrm P(x_0,t_0 ) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}) \,.
$$

Therefore the problem of calculating the likelihood for a particular trajectory amounts to solving two independent problems:

1. estimating the probability density of the starting point $\mathrm P (x_0, t_0)$ of a response
2. calculating the transition probabilities $\mathrm P(x_n,t_n|x_{n-1},t_{n-1})$

For a given chemical reaction network we can write down the chemical master equation. The chemical master equation contains all the information needed to compute the individual terms $\mathrm P(x_n,t_n|x_{n-1},t_{n-1})$ for the entire system.

To calculate the mutual information between $\mathcal{S}$ and $\mathcal X$ we have to consider the entire reaction network containing the components both in $S$ and in $X$. The precise reaction dynamics of the response part of the chemical network crucially depend on the observed signal trajectory. Therefore the chemical master equation for the whole reaction network allows us to compute the likelihood of a response trajectory for a particular signal trajectory:

$$
\mathrm P(\mathcal X = X|\mathcal S = S) = \mathrm P(x_0,t_0;x_1,t_1;\ldots;x_{N-1},t_{N-1} | S) = \mathrm P(x_0,t_0 | S) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S) \,.
$$

For increasingly long trajectories this quantity will in many physically relevant cases either grow or decay exponentially (*TODO: explain why*). Thus sufficiently long trajectories, the numerical values of the likelihood will not be directly representable by conventional floating-point numbers.

This problem can be avoided if we compute the *log-likelihood* $\ell(X|S) \equiv \ln\mathrm P(X|S)$ instead. We can easily rephrase the equation for the likelihood:

$$
\ell(X|S) = \ln\left[ \mathrm P(x_0,t_0 | S) \prod\limits^{N-1}_{n=1} \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S) \right] = \ln \mathrm P(x_0,t_0 | S) +\sum\limits^{N-1}_{n=1} \ln \mathrm P(x_n,t_n|x_{n-1},t_{n-1}, S)\,.
$$

### The probability density for the starting point of a trajectory

We first look at the term $P_0 = \mathrm P(x_0,t_0 | S)$. Since $S$ is a trajectory in time we can directly conclude from causality that

$$
\mathrm P(x_0, t_0 | S) = \mathrm P(x_0, t_0 | S_{t \leq t_0})
$$

where $S_{t \leq t_0}$ is the temporal piece of the signal up to $t_0$. We further suppose that the signal itself is markovian and therefore has no memory of its past. With this simplification we get

$$
\mathrm P(x_0, t_0 | S) = \mathrm P(x_0, t_0 | S_{t = t_0}) = \frac{\mathrm P((x_0, t_0), (s_0, t_0))}{\mathrm P(s_0, t_0)} \,.
$$

We estimate $P_0$ using gaussian kernel density estimation to approximate both, the joint distribution of $X_0, S_0$ and the marginal distribution of $S_0$.

Knowing the probabilities of the initial condition of both response and signal we can directly estimate the mutual information of $\mathcal{X}_{t=t_0}$ and $\mathcal{S}_{t=t_0}$:

$$
\mathrm I(\mathcal{X}_{t=t_0}, \mathcal{S}_{t=t_0}) = \int ds_0\int dx_0\; \mathrm{P}(x_0, s_0)\; \ln \frac{\mathrm{P} (x_0, s_0)}{\mathrm{P} (x_0) \mathrm{P} (s_0)}
$$

### Estimating the marginal probability of response trajectories

To calculate the mutual information between trajectories we need to have a good estimate for $\ln\left\langle \mathrm P(X | S) \right\rangle_\mathcal{S}$. We calculate this average by sampling of trajectories $(S^{(i)})_{i=1\ldots N_S}$ from the probability distribution of $\mathcal{S}$:

$$
\ln\left\langle\mathrm P(X | S) \right\rangle_\mathcal{S} \approx \ln \frac{\sum^{N_S}_{i=1} \mathrm P(X|S^{(i)})}{N_S} = \ln \sum^{N_S}_{i=1} \mathrm P(X|S^{(i)}) - \ln N_S
$$

Thus we find that it is enough to be able to compute the likelihood between trajectories to estimate the marginal distribution of trajectories.



We approximate the mutual information between trajectories as

$$
\mathrm{I}(\mathcal{X}; \mathcal{S}) = \left\langle \ln \frac{\mathrm{P} ( X |  S)}{\left\langle\mathrm P(X | S) \right\rangle_\mathcal{S}} \right\rangle_{\mathcal{X},\mathcal{S}} = \left\langle \ell ( X |  S) - \mathrm{LSE}\left( \ell(X|S^{(1)}),\ldots, \ell(X|S^{(N_S)})\right) + \ln N_S\right\rangle_{\mathcal{X},\mathcal{S}}
$$

which means that for finite amount of signal samples we will _systematically over-estimate_ the mutual information. **TODO: Is that really true? What about $\ln N_S$?** Even worse: the longer the trajectories the bigger the error becomes since the dimensionality of the space of possible signals is growing.

Another way to phrase this insight is that to get a good approximation for the logarithmic average likelihood, our set of signals that we use for monte-carlo sampling should contain many signals that produce a high likelihood. **Therefore it probably is necessary to come up with a scheme to specifically sample signal trajectories for which the likelihood of a particular trajectory is high**. On the other hand the results do not seem to get significantly better when averaging over more trajectories. -->

## References


