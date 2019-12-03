# Mutual Information between Trajectories

## Information theory for trajectories

A trajectory $X$ with $N$ steps is defined by a set of pairs $X=\{(t_i, \mathbf{x}_i)\; |\; i=0\ldots N-1 \}$ where $\mathbf{x}_i$ defines the trajectory value at time $t_i$. We can also have random variables over trajectories and therefore probability distributions over the space of all trajectories.

As a next step we can make sense of the entropy of a trajectory. Let $\mathcal{X}_N$ be a random variable over trajectories of length $N$. We call

$$
\mathrm H(\mathcal{X}_N) = - \int\limits_{X\in \sigma(\mathcal{X}_N)} dX\; \mathrm{P}(\mathcal{X}_N = X)\; \ln \mathrm{P} (\mathcal{X}_N = X)
$$

the entropy of $\mathcal{X}_N$ where $\mathrm{P}(\mathcal{X}_N = X)$ is the probability density function of a trajectory $X=\{(t_i, \mathbf{x}_i)\; |\; i=0\ldots N-1 \}$. We can also define the conditional entropy for trajectories

$$
\mathrm H(\mathcal{X}_N | \mathcal{S}_M) = -\int\limits_{S\in \sigma(\mathcal{S}_N)} dS\; \mathrm{P} (\mathcal{S}_N = S) \int\limits_{X\in \sigma(\mathcal{X}_N)} dX\; \mathrm{P}(\mathcal{X}_N = X | \mathcal{S}_N = S)\; \ln \mathrm{P} (\mathcal{X}_N = X | \mathcal{S}_N = S) \,.
$$

With these two quantities we can express the _mutual information_ between trajectories

$$
\mathrm{I}(\mathcal{X}_N; \mathcal{S}_M) = \mathrm H(\mathcal{X}_N) - \mathrm H(\mathcal{X}_N | \mathcal{S}_M) \,.
$$

The mutual information between two random variables quantifies by how much our certainty of the value of one variable increases if we know the other one.

To shorten the notation we write $\mathrm{P} (\mathcal{X}_N = X)$ as $\mathrm P_{\mathcal{X}_N}(X)$ and if the random variable is clear from the context we even drop the index and only write $\mathrm P(X)$. With the short notation we can rewrite the mutual information

$$
\mathrm{I}(\mathcal{X}_N; \mathcal{S}_M) = \int\limits_{S\in \sigma(\mathcal{S}_N)} dS \int\limits_{X\in \sigma(\mathcal{X}_N)} dX\; \mathrm{P}( X , S)\; \ln \frac{\mathrm{P} ( X |  S)}{\mathrm P(X)} \,.
$$

To evaluate $\mathrm P(X)$ we have to expand it as follows

$$
\mathrm P(X) = \int\limits_{S\in \sigma(\mathcal{S}_N)} dS\; \mathrm P(X, S) = \int\limits_{S\in \sigma(\mathcal{S}_N)} dS\; \mathrm P(X|S) \ \mathrm P (S) \equiv \left\langle P(X | S) \right\rangle_{\mathcal{S}_N} \,.
$$

These relations let us state the mutual information as nested averages over the likelihood $P(X|S)$:

$$
\mathrm{I}(\mathcal{X}; \mathcal{S}) = \left\langle \ln \frac{\mathrm{P} ( X |  S)}{\mathrm P(X)} \right\rangle_{\mathcal{X},\mathcal{S}} = \left\langle \ln \frac{\mathrm{P} ( X |  S)}{\left\langle\mathrm P(X | S) \right\rangle_\mathcal{S}} \right\rangle_{\mathcal{X},\mathcal{S}} \,.
$$

These averages are defined as integrals over the very high-dimensional space of trajectories and thus very hard to evaluate analytically or numerically in the general case. Our goal is use *Monte-Carlo sampling* in the trajectory space to evaluate the above averages. To do this we have to sample trajectories from their probably distribution and we need to evaluate the likelihood for a response given a signal.

## Monte-Carlo simulation

While so-called Monte-Carlo methods comprise a wide variety of approaches to stochastically evaluate integrals or sums the common idea is easily stated. We have a state space $U$ and a probability distribution $p_U$ over that state space. The problem is to evaluate

$$
\langle f(u) \rangle \equiv \int\limits_{u \in U} \mathrm du\; f(u) p_U(u)
$$

where $f: U\rightarrow\mathbb R$ is some smooth function. If $U$ is high-dimensional it is very time-consuming to estimate it by direct numerical integration.

## Estimating the likelihood

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

### The transition probabilities

**TODO:** describe how the transition rates follow from the master equation (probably follow the style of Cepeda-Humerez, et. al.)

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

## Simulating chemical networks

As a model for the biochemical processing that takes place inside a cell we suppose that all interactions can be described by a chemical networks composed of a number of components (different types of molecules) and a number of reactions between them. As a very basic example we consider a highly simplified model of gene expression consisting of two components and four reactions:

\begin{align}
\emptyset &\longrightarrow S \longrightarrow \emptyset \\
S &\longrightarrow S + X \\
X &\longrightarrow \emptyset \,.
\end{align}

We might interpret $S$ as some signal whose quantity varies stochastically. Then there is a certain chance that a signal molecule is registered by the cell which triggers the creation of an $X$. The amount of $X$ then also decays over time. We call the trajectory of $S$ the "signal" and the trajectory of $X$ the "response".
