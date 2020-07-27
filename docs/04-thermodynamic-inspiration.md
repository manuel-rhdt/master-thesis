---
author:
  - Manuel Reinhardt
  - Pieter Rein ten Wolde
title: Notes for Estimates Based on Statistical Physics
institute: AMOLF
bibliography: ["library.bib"]
link-citations: true
linkReferences: true
autoEqnLabels: true
cref: true
lang: en-US
---

# Directed Sampling in Trajectory Space

In the previous chapters we established a technique to compute the mutual information between time-variable signals and responses using Monte-Carlo integration together with stochastic simulations. In practice we found that this method often severely over-estimates the mutual information. In the previous chapter we came to the conclusion that the root of the biased estimates is the computation of the marginal probability density $\mathrm P(\mathbf x)$.

The fundamental identity for the prediction of model parameters $\mathbf s$ given some data or response $\mathbf x$ is Bayes' theorem
$$
\mathrm P(\mathbf s | \mathbf x) = \frac{\mathrm P(\mathbf x|\mathbf s)\ \mathrm P(\mathbf s)}{\mathrm P(\mathbf x)}
$$ {#eq:bayes_thm}
which describes how the observation of $\mathbf x$ changes the probability distribution for the signals from $\mathrm P(\mathbf s)$ to $\mathrm P(\mathbf s | \mathbf x)$. In statistics, one usually denotes $\mathrm P(\mathbf s)$ as the _prior_, $\mathrm P(\mathbf x|\mathbf s)$ as the _likelihood_, $\mathrm P(\mathbf x)$ as the _evidence_ (or _marginal likelihood_) and $\mathrm P(\mathbf s | \mathbf x)$ as the _posterior_ distribution. It is no accident that for the description of statistical inference we reuse the same letters $\mathbf x$ and $\mathbf s$ as in the previous chapters. Indeed, at its core, the biochemical network performs statistical inference of the signal from the data given by its response. Thus, for now, we are going to change the perspective and look at the problem through the eyes of a Bayesian statistician. 

From the statistician's point of view, we have a _model_ with a set of parameters $\mathbf s = (s_1, s_2,\ldots)^T$ that we want to fit to the observations or data $\mathbf x$. The Bayesian approach involves choosing a reasonable _prior_, i.e. making an educated guess for the distribution $\mathrm P(\mathbf s)$ using justifiable assumptions about the underlying system. Using @eq:bayes_thm, the statistician can compute the _posterior_ probability distribution of the model parameters $\mathrm P(\mathbf s | \mathbf x)$ from the observations $\mathbf x$. This involves computing the _likelihood_ $\mathrm P(\mathbf x|\mathbf s)$ and the _evidence_ $\mathrm P(\mathbf x)$ of the acquired data $\mathbf x$ for the given model. Usually, models are chosen that make it easy to compute the likelihood of some data given the model parameters. It is the computation of the _evidence_ $\mathrm P(\mathbf x)$ which usually requires more sophisticated computational methods @2014.Held. In most cases, the evidence is estimated using the identity $\mathrm P(\mathbf x)=\int\mathrm d\mathbf s\ \mathrm P(\mathbf s, \mathbf x)$ which has motivated the use of the term _marginal likelihood_. In the following section we show the evolution of methods to compute such integrals in statistical science.

## Previous Work on the Computation of the Marginal Likelihood

The estimation of the marginal density $\mathrm P(\mathbf x)=\int\mathrm d\mathbf s\ \mathrm P(\mathbf s, \mathbf x)$ is an important ingredient for Bayesian statistics, and is used for example in the computation of so-called _Bayes factors_ for hypothesis testing. 

Written as
$$
\mathrm P(\mathbf x) = \langle \mathrm P(\mathbf x|\mathbf s) \rangle_{\mathrm P(\mathbf s)}
$$ {#eq:marginal_as_mean}
the marginal density is an example of an expectation value with respect to a potentially complex distribution. For simple distributions, expectations can be approximated by a simple Monte-Carlo approach, taking the sample average $1/N \sum^N_{i=1} \mathrm P(\mathbf x|\mathbf s_i)$ of $N$ independent random samples $\mathbf s_1,\ldots,\mathbf s_N$ generated from $\mathrm P(\mathbf s)$. That is, we estimate the marginal density using samples from the _prior_ and averaging over their likelihoods. As demonstrated in the previous chapters, this method is not very efficient when the individual signals $\mathbf s_i$ are actually trajectories. Once the space of possible signals is sufficiently large it becomes extremely unlikely to sample a $\mathbf s^\prime$ such that $\mathrm P(\mathbf x|\mathbf s^\prime)$ presents a strong contribution to the mean. A well-known approach to circumvent this problem is _importance sampling_.

### Importance Sampling

In importance sampling, we choose a different distribution than the prior from which to generate samples. Let $\hat{\mathbf s}_1, \ldots,\hat{\mathbf s}_N$ be samples drawn from a distribution with a known density function $\hat q$. Then we can estimate @eq:marginal_as_mean as
$$
\hat{\mathrm P}^{(q)}(\mathbf x) = \frac{1}{N} \sum^N_{i=1} \frac{\mathrm P(\hat{\mathbf s}_i)\,\mathrm P(\mathbf x | \hat{\mathbf s}_i)}{\hat q(\hat{\mathbf s}_i)} \equiv 
\frac{1}{N} \sum^N_{i=1} \frac{q(\hat{\mathbf s}_i)}{\hat q(\hat{\mathbf s}_i)}
$$ {#eq:importance_sampling}
where by $q(\mathbf s)=\mathrm P(\mathbf s)\,\mathrm P(\mathbf x|\mathbf s)$ we denote the unnormalized posterior distribution (note that the posterior is $\mathrm P(\mathbf s|\mathbf x) = q(\mathbf s)/C$ with $C=\mathrm P(\mathbf x)$).
This is precisely the approach that was described in @sec:umbrella using the term _umbrella sampling_. For this section we will continue to use the statisticians' terminology of _importance sampling_. The choice of $\hat q$ requires crucial consideration for @eq:importance_sampling to yield a good estimate. To reduce the variance of the ratio $q(\hat{\mathbf s}_i) / \hat q(\hat{\mathbf s}_i)$, the sampling distribution $\hat q$ should be similar to the posterior and have tails no thinner than $q$ @1997.Diciccio. In practice it is often very difficult to find an adequate choice for $\hat q$, especially in high-dimensional spaces. Therefore, it has been suggested to choose an approproate sampling function using one or more _posterior samples_ i.e. random samples taken from the posterior distribution $\mathrm P(\mathbf s|\mathbf x)$ [@1997.Diciccio;@2011.Chan;@2014.Chan;@2014.Perrakis]. Even if it is possible to acquire a sufficient amount of posterior samples, these methods usually require the choice of an appropriate family of distributions from which the sampling distribution is picked. Since in our case, the individual $\mathbf s$ are entire trajectories it is unclear what a useful family of distributions would be.

Some authors have instead proposed methods to compute the marginal density _directly_ from the posterior samples without the need to pick a separate sampling distribution $\hat q$. This promises to be more efficient than using prior samples since we don't face the problem of generating mostly irrelevant signals for a given response.
However, the generation of random samples from $\mathrm P(\mathbf s|\mathbf x)$ is often a non-trivial problem. Thus, there have been proposed a variety of Markov chain Monte Carlo methods to generate approximately independent samples from the posterior distribution [@1994.Tierney;@1990.Gelfand;@1987.Tanner;@2001.Neal;@2018.Warne]. We will leave the discussion of efficient posterior sampling for later and assume for now to have a set of corresponding samples available.

### Estimating the Marginal Density from Posterior Samples

Let $\tilde{\mathbf s}_1,\ldots,\tilde{\mathbf s}_N$ be independent posterior draws. Newton, et. al. @1994.Newton propose the harmonic mean estimator
$$
\hat{\mathrm P}^{(\text{harm.})}(\mathbf x) = \left[ \frac{1}{N} \sum^N_{i=1} \mathrm P(\mathbf x|\tilde{\mathbf s}_i)^{-1} \right]^{-1}\,,
$$ {#eq:harmonic_rule}
i.e. the marginal likelihood is estimated by the harmonic mean of the likelihoods of a sample from the posterior distribution. It has been noted that this estimate is rather unstable because of the occasional occurence of a value $\tilde{\mathbf s}_i$ with a small likelihood and hence a large effect on the final result @1994.Newton. The estimate in @eq:harmonic_rule can therefore be improved using an arbitrary probability density $f$ as proposed by Gelfand, et. al. @1994.Gelfand
$$
\hat{\mathrm P}^{(\text{reciprocal})}(\mathbf x) = \left[ \frac{1}{N} \sum^N_{i=1} \frac{f(\tilde{\mathbf s}_i)}{q(\tilde{\mathbf s}_i)} \right]^{-1}\,.
$$ {#eq:reciprocal_importance}
If we choose $f(\mathbf s)=\mathrm P(\mathbf s)$ we recover the harmonic rule in @eq:harmonic_rule. If instead, we choose $f$ with tails thinner than $q$ this estimator is well-behaved @1997.Diciccio. This requirement for $f$ is the opposite as for importance sampling which lead to the estimate in @eq:reciprocal_importance being called _reciprocal importance sampling_. However, similarly as before it is often difficult to find a density $f$ that has sufficiently thin tails in all dimensions. Both, importance sampling and reciprocal importance sampling are special cases of a more general identity.

### Bridge Sampling and Beyond

The idea of bridge sampling was first given by Bennet @1976.Bennett for the simulation of free energy differences. Given two unnormalized probability densities, $q_1$ and $q_2$, Meng and Wong @1996.Meng state the central identity of bridge sampling as
$$
\frac{c_1}{c_2} = \frac{\langle q_1(\mathbf s)\,\alpha(\mathbf s) \rangle_2}{\langle q_2(\mathbf s)\,\alpha(\mathbf s) \rangle_1}
$$ {#eq:bridge_sampling}
where $c_1, c_2$ are the normalization constants corresponding to $q_1, q_2$, and $\alpha$ is an arbitrary function. $\langle\,\rangle_i$ denotes the expectation with respect to the distribution corresponding to $q_i$. If we choose $q_1(\mathbf s)=q(\mathbf s)$ and $q_2(\mathbf s)=\hat q(\mathbf s)$ then for $\alpha(\mathbf s)=1/\hat q(\mathbf s)$ we get importance sampling as defined in @eq:importance_sampling. If we choose $\alpha(\mathbf s) = [q_1(\mathbf s)\,q_2(\mathbf s)]^{-1}$ we get
$$
\frac{c_1}{c_2} = \frac{\langle q^{-1}_2(\mathbf s) \rangle_2}{\langle q^{-1}_1(\mathbf s) \rangle_1}
$$
which is a generalization of the harmonic rule in @eq:harmonic_rule.

For bridge sampling, Bennet, Meng and Wong [@1976.Bennett;@1996.Meng] discuss optimal choices for $\alpha$ such that the mean square error of the estimate is minimized. It turns out that a good choice of $\alpha$ acts as a _bridge_ between the densities $q_1$ and $q_2$, such that both the numerator and the denominator of @eq:bridge_sampling can be reliably computed. Hence, the name: bridge sampling.

It did not go unnoticed that techniques originally developed for the computation of free energy differences might be useful in Bayesian computation. Gelman and Weng @1998.Gelman introduce _path sampling_ which is a generalized form of _thermodynamic integration_. It takes the idea of bridge sampling one step further: instead of using one bridge to join two distributions $q_1$ and $q_2$ and compute the ratio of their normalization constants, thermodynamic integration constructs a continuous _path in distribution space_ between $q_1$ and $q_2$ to further increase the efficiency. Another related method is _annealed importance sampling_, proposed by Neal @2001.Neal. It is based on simulated annealing, another technique originating in statistical physics.

We have seen that many methods to compute the marginal likelihood draw inspiration from statistical physics. Specifically we can map the computation of the normalization constant of a probability distribution to the simulation of free energy differences. In the next subsection we will introduce notation and terminology reminiscent of Statistical physics to describe our problem. In this way we can understand why those methods are so helpful for Bayesian computation. Using that description we show that both, Thermodynamic Integration and the Wang-Landau algorithm are very promising candidates for the computation of $\mathrm P(\mathbf x)$.


## Borrowing Terminology from Statistical Physics

The mutual information then measures the difference in entropies between these two distributions or in other words: _how much uncertainty about $\mathbf s$ the observation of $\mathbf x$ was able to remove_.  For the computation of the _mutual information_ between $\mathcal S$ and $\mathcal X$ we need to estimate the marginal entropy $\mathrm H(\mathcal X)$ which requires averaging over the logarithm of the marginal densities $\ln\mathrm P(\mathbf x_i)$ for many samples $\mathbf x_i\sim\mathcal X$. Since we only have access to the _prior_ and the _likelihood_ we have to compute the marginal likelihood as $\mathrm P(\mathbf x) = \int \mathrm d\mathbf s\ \mathrm P(\mathbf s | \mathbf x) \mathrm P(\mathbf s)$. These kinds of integrals are very familiar to researchers in statistical physics albeit they use a different terminology to describe them. In the following we will describe how @eq:bayes_thm is analogous to the distribution of the canonical ensemble and how this analogy allows us to make use of efficient computational methods from statistical physics.

In the framework employed by statistical physics, Bayes' theorem corresponds to the canonical ensemble distribution of $\mathbf s$ (for $\beta=1$)
$$
\mathrm P(\mathbf s | \mathbf x) = \frac{1}{Z(\mathbf x)}\exp\left[-E(\mathbf s, \mathbf x)\right]
$$
where the _partition function_ is defined by $Z(\mathbf x) = \int \mathrm d\mathbf s\ \exp\left[-E(\mathbf s, \mathbf x)\right]$ and $E(\mathbf s, \mathbf x)$ denotes the total energy of the system at state $\mathbf s$. In this context $\mathbf x$ is considered a parameter vector for the specific model used to compute the energy. In classical problems of statistical physics (such as e.g. the _Ising model_) the state space spans the single particle states $\mathbf{s} = (\sigma_1,\ldots,\sigma_n)\in\Omega^n$ for all particles and the energy is given by the _Hamiltonian_ $\mathcal H(\sigma_1,\ldots,\sigma_n;\mathbf x)$ where $\mathbf x$ could contain parameters describing e.g. the interaction strength between neighboring spins. In our case however we define our energy function by comparison with @eq:bayes_thm as
$$E(\mathbf s, \mathbf x) = -\ln\mathrm P(\mathbf x|\mathbf s)-\ln\mathrm P(\mathbf s)\,.$$ 
From this point of view the marginal density $\mathrm P(\mathbf x) = Z(\mathbf x)$ _is_ the partition function of the canonical ensemble. In statistical physics the partition function is of central importance since its partial derivatives include all thermodynamic properties of a physical system. The free energy of the canonical ensemble is defined by $\mathcal F(\mathbf x) = -\ln Z(\mathbf x)$ (for $\beta=1$) such that using this terminology we can write the marginal entropy $\mathrm H(\mathcal X)$ as an average over the _"free energies of response trajectories"_
$$
\mathrm H(\mathcal X) = \int\mathrm d\mathbf x\ \mathrm P(\mathbf x)\ \mathcal F(\mathbf x) = \left\langle \mathcal F(\mathbf x) \right\rangle_{\mathrm P(\mathbf x)}\,.
$$

Since the computation of the partition function is central to the solution of many statistical problems there has been done considerable work on efficient estimation of the partition function, the free energy and other related quantities such as the _density of states_.

## Thermodynamic Integration

One well-established technique to estimate free energy (differences) is by thermodynamic integration (TI) @1998.Gelman. It allows the accurate computation of the ratio between the normalization constants of two different probability distributions using a continuous path in _distribution space_ that connects both. Since this strategy uses random samples taken from many different distributions along this path it is especially robust when the two distributions have very little overlap. For the computation of the marginal density $\mathrm P(\mathbf x)$ we can (for a given $\mathbf x$) define a suitable path in distribution space between $\mathrm P(\mathbf s)$ and $\mathrm P(\mathbf s, \mathbf x)$. The normalization constants of these distributions are $z_0 = 1$ and $z_1 = \mathrm P(\mathbf x)$, respectively such that the ratio $r=z_1/z_0$ of these normalization constants directly corresponds to the marginal density. Using TI we estimate this ratio using approximately independent samples from a _Markov chain Monte Carlo_ (MCMC) simulation.

In the following sections we will give a quick summary of TI followed by an explanation of the Markov chain Monte Carlo simulation and a discussion of the resulting accuracy of the estimates.

### Summary of the Technique

Let $q_0$ and $q_1$ be the unnormalized distribution functions and $z_0, z_1$ the corresponding normalization constants such that $z_i=\int\mathrm d\mathbf s\ q_i(\mathbf s)$. Next we construct a path between $q_0$ and $q_1$, parametrized by $\theta\in[0,1]$ such that $q_\theta$ smoothly connects the end points. We similarly define $z(\theta)$ as the normalization constant of $q_\theta$. A smooth path that can be constructed for any pair of distributions $(q_0, q_1)$ is the _geometric path_ given by $q_\theta=q^{1-\theta}_0\ q^\theta_1$. Note however that variance of the estimate depends on the chosen path and that the geometric path is not the optimal path in general. 

For the estimation of free energy differences we are interested in the ratio $r=z(1)/z(0)$. To find an estimate we differentiate the logarithm of $z(\theta)$ with respect to $\theta$ to arrive at
$$
\frac{\mathrm d\ln z(\theta)}{\mathrm d\theta} = \frac{1}{z(\theta)} \frac{\partial}{\partial\theta}  \int\mathrm d\mathbf s\ q_\theta(\mathbf s) = \int\mathrm d\mathbf s\ \frac{q_\theta(\mathbf s)}{z(\theta)} \frac{\partial}{\partial\theta} \ln q_\theta(\mathbf s) = \left\langle \frac{\partial}{\partial\theta} \ln q_\theta(\mathbf s) \right\rangle_{p_\theta(\mathbf s)}
$$
where $p_\theta(\mathbf s) = q_\theta(\mathbf s)/z(\theta)$ is the normalized probability distribution corresponding to $q_\theta$. By analogy to the potential in statistical physics we define
$$
U(\mathbf s, \theta) = -\frac{\partial}{\partial\theta} \ln q_\theta(\mathbf s)\,.
$$
Now we can express the log-ratio $\lambda=\ln r$ by the integral
$$
\lambda = \ln z(1) - \ln z(0) = -\int\limits^1_0 \mathrm d\theta\ \left\langle 
U(\mathbf s, \theta)
\right\rangle_{p_\theta(\mathbf s)}
$$ {#eq:path_sampling_int}
which forms the basis of all thermodynamic integration estimates. One advantage of the TI estimators is that we directly estimate the log-ratio $\lambda$, i.e. the free energy difference as opposed to the ratio of partition functions. Indeed, to eventually compute the marginal entropy $\mathrm H(\mathcal X) = -\langle\ln P(\mathbf x)\rangle$ we require the logarithm of the marginal density thus no further error is introduced by taking the logarithm of an estimated quantity.

Using the previous identities, one possible way to estimate $\lambda$ is to regard $\theta$ as a random variable with a density $p(\theta)$, allowing us to compute the integral in @eq:path_sampling_int using the Monte Carlo estimator
$$
\hat{\lambda}_\text{MC} = -\frac{1}{n} \sum\limits^n_{i=1}\frac{U(\mathbf s_i, \theta_i)}{p(\theta_i)}
$$ {#eq:lambda_mc}
with draws $(\mathbf s_1, \theta_1),\ldots,(\mathbf s_n, \theta_n)$ from the joint probability density $p(\mathbf s, \theta) = p_\theta(\mathbf s)\ p(\theta)$. Alternatively, we can perform numerical integration using the trapezoidal rule by evaluating the potential over values $\theta_1<\cdots<\theta_{n-1}$ between $\theta_0=0$ and $\theta_n=1$
$$
\hat{\lambda}_\text{NI} = -\sum\limits^n_{i=1}\frac{
  \langle U(\mathbf s, \theta_{i-1}) \rangle_{p_{\theta_{i-1}}(\mathbf s) }
  + \langle U(\mathbf s, \theta_{i}) \rangle_{p_{\theta_{i}}(\mathbf s) }
  }{2} (\theta_i - \theta_{i-1})
$$ {#eq:lambda_ni}
where each average over the potential is performed using a Monte Carlo simulation.

To use these estimators for the computation of the marginal density $\mathrm P(\mathbf x)$ at a given $\mathbf x$ we need to construct a path between the densities $q_0(\mathbf s) = \mathrm P(\mathbf s)$ and $q_1(\mathbf s) = \mathrm P(\mathbf s)\mathrm P(\mathbf x|\mathbf s)$. For simplicity and convenience we choose the geometric path $q_\theta(\mathbf s) = \mathrm P(\mathbf s)\ [\mathrm P(\mathbf x|\mathbf s)]^\theta$. Taking the logarithm of this density we get $\ln q_\theta(\mathbf s) = \ln P(\mathbf s) + \theta \ln \mathrm P(\mathbf x|\mathbf s)$ which prompts us to define the "energy" of a signal trajectory with respect to $\theta$ as
$$
E(\mathbf s, \theta) = -\ln q_\theta(\mathbf s) = -\ln P(\mathbf s) - \theta \ln \mathrm P(\mathbf x|\mathbf s)\,.
$$
For $\theta = 1$ this definition of the energy matches our previous definition by analogy with the canonical ensemble whereas for $\theta = 0$ this definition of the energy is equivalent to the energy of a signal trajectory for a system where $\mathcal S$ and $\mathcal X$ are completely independent and thus $\mathrm P(\mathbf s|\mathbf x) = \mathrm P(\mathbf s)$. Thus, this nomenclature also motivates the name "potential" for the quantity
$$
U(\mathbf s, \theta) = - \frac{\partial}{\partial\theta} \ln q_\theta(\mathbf s) = \frac{\partial}{\partial\theta} E_\theta(\mathbf s) = -\ln \mathrm P(\mathbf x|\mathbf s)
$$
i.e. $\theta$ acts as a _"knob"_ that allows us to gradually turn the potential on or off. The potential term itself characterizes the amount of dependence between the random variables $\mathcal S$ and $\mathcal X$. Note that all energetic quantities depend on the specific response $\mathbf x$ (except at $\theta=0$) even if this dependence is suppressed in the notation.

To use the TI estimators introduced in [@eq:lambda_mc;@eq:lambda_ni] we need to generate samples from arbitrary distributions along our chosen geometric path. Since we can compute the unnormalized densities of these distributions, we can use the Metropolis-Hastings algorithm as a very general method to sample from arbitrary distributions @1970.Hastings.

### Markov Chain Monte Carlo

To generate approximately independent samples from a distribution given by the unnormalized density $q_\theta$ we start from an (in principle arbitrary) initial signal $\mathbf s$. Next, a new signal $\mathbf s^\prime$ is proposed from the proposal distribution $\mathrm T(\mathbf s \rightarrow \mathbf s^\prime)$ which is typically chosen to yield a $\mathbf s^\prime$ close to $\mathbf s$. Then with some probability $A(\mathbf s^\prime, \mathbf s)$ we _accept_ the new configuration and our first generated sample is $\mathbf s_1 = \mathbf s^\prime$. Otherwise we _reject_ the new configuration and our first sample is equal to the initial signal $\mathbf s_1 = \mathbf s$. For the next iteration of the algorithm we then set our new initial signal to be $\mathbf s \leftarrow \mathbf s_1$ such that when we repeat this procedure many times we generate a sequence of signals $\mathbf s_1, \mathbf s_2, \ldots$ where each sample is a random value only directly dependent on the immediately preceding sample. Thus we have defined a Markov process that generates a _chain_ of signals with the transition probability given by $\mathrm P(\mathbf s \rightarrow \mathbf s^\prime) = T(\mathbf s \rightarrow \mathbf s^\prime)\,A(\mathbf s^\prime, \mathbf s)$. We want to choose the acceptance probability $A(\mathbf s^\prime, \mathbf s)$ such that the stationary distribution of this Markov process is precisely $q_\theta$. It can be shown that the _Metropolis choice_
$$
A(\mathbf s, \mathbf s^\prime) = \min\left( 1, \frac{q_\theta(\mathbf s^\prime)}{q_\theta(\mathbf s)} \frac{\mathrm T(\mathbf s^\prime \rightarrow \mathbf s)}{\mathrm T(\mathbf s \rightarrow \mathbf s^\prime)} \right)
$$ {#eq:metropolis_acceptance}
leads to the correct stationary distribution given that the system is ergodic.

<!-- TODO: While this algorithm has some disadvantages (dependence of samples, yada yada) it often is the only sampling strategy that works at all in very high-dimensional spaces or complex distributions (is also well parallelizable)... -->

![Visualization of the log-likelihoods $\ln\mathrm P(\mathbf x|\mathbf s)$ and the acceptance rates $\frac{\text{accepted}}{\text{accepted} + \text{rejected}}$ during individual Monte-Carlo runs. Each MCMC run used a different value for $\theta$, randomly chosen from the uniform distribution between 0 and 1. The $x$-axis shows the sample numbers 1 to 100. Between one sample and the next one we skip 1000 accept-reject steps for which no statistics were collected.](figures/monte_carlo_sims.svg){#fig:monte_carlo_sims}

![Comparison of the covariance matrices obtained a) by computing the empirical covariance of 1000 approximately uncorrelated samples taken from the MCMC procedure (for $\theta = 1$) and b) by analytically computing the covariance matrix of the normal distribution $\mathrm P(\mathbf s|\mathbf x)$. The proposal distribution is a multivariate normal distribution with covariance $\Sigma=\sigma^{-2} \mathbb I$, with a value of $\sigma=0.01$.](figures/mcmc_covariance_comparison.svg){#fig:mcmc_covariance}

For the Gaussian system we choose the proposal distribution $\mathrm T(\mathbf s \rightarrow \mathbf s^\prime)$ to be a multivariate normal distribution centered around $\mathbf s$ and with uniform covariance $\Sigma=\sigma^{-2} \mathbb I$. In @fig:mcmc_covariance we show that using the Metropolis-Hastings algorithm we can generate samples with an appropriate distribution that matches the analytical expectation. In @fig:thermodynamic_int_results we show the averaged potentials for 216 MCMC runs for different values of $\theta$. From these potentials we can the compute the marginal density using the estimator from @eq:lambda_mc. The results are very promising since the estimated value differs by merely 0.012 % from the analytically correct value of $\mathrm P(\mathbf x)$.

![Samples of the averaged potential for different values of $\theta$. There are 216 samples for values of $\theta$ chosen uniformly distributed in the interval $[0, 1]$. Every point is an individual MCMC simulation with 1000 approximately independent draws. The bars on the right show a histogram of the log-likelihoods. The TI estimate is the integral from $\theta=0$ to $1$ of the curve that the individual samples approximate. Using the estimate from @eq:lambda_mc, the estimated value differs by merely 0.012 % from the analytically correct value of $\mathrm P(\mathbf x)$. This shows that given enough samples, TI is able to provide very accurate results for the marginal density.](figures/mcmc_theromdynamic_integration.svg){#fig:thermodynamic_int_results}


## Estimating the Density of States

<!-- The goal of the Wang and Landau algorithm is to compute the density of states $\rho(E)$ for a system using an adaptive Metropolis scheme where the sampling distribution is changing throughout one simulation. We want to show that this algorithm can be used to get a better estimate of the marginal probability density for random trajectories. -->

In the context of statistical physics we often look at configurations of a system that can be described by a configuration $\mathbf{n}\in\Omega$ where $\Omega$ is the state space of the system. We can typically assign a probability (density) to each configuration. For example, let's consider the canonical ensemble for a given inverse temperature $\beta$ and Hamiltonian $\mathcal H$
$$
\mathrm{P}(\mathbf{n}) = \frac{1}{Z(\beta)} e^{-\beta \mathcal H(\mathbf{n})}
$$ {#eq:canonical_probability}
with the _partition function_ $Z(\beta)=\int \mathrm{d}\mathbf{n}\ e^{-\beta H(\mathbf{n})}$. The Hamiltonian assigns an energy to every state, i.e. for every state $\mathbf{n}$ we have an associated energy $\mathcal H(\mathbf{n})$. To learn more about the distribution of energies in our system we can now define the _density of states_ $g(E)$ which is proportional to the number of states with energy $E$. More precisely, let $\mathcal{N}$ be a random variable uniformly distributed in the state space, then
$$
g(E) = \mathrm{P}\left(\mathcal H(\mathcal N) = E\right)\,.
$$

While the density of states (_DOS_) thus describes the number of states at individual energy levels, the Boltzmann factor $e^{-\beta E}$ specifies the relative weight of states with energy $E$ in the canonical ensemble. That is, we can compute the ensemble average $\langle f(\mathcal H(\mathbf{n}))\rangle$ of any function $f$ that depends only on the energy of a given state as
$$
\langle f(\mathcal H(\mathbf{n}))\rangle = 
\frac{\int_\Omega\mathrm d\mathbf n\ f(\mathcal H(\mathbf n)) e^{-\beta \mathcal H(\mathbf{n})}}{\int_\Omega\mathrm d\mathbf n\ e^{-\beta \mathcal H(\mathbf{n})}} = 
\frac{
\int\mathrm dE\ f(E)\ g(E) e^{-\beta E}
 }{
   \int\mathrm dE\ g(E) e^{-\beta E}
 }
$$ {#eq:def_dos}
where $\int_\Omega\mathrm d\mathbf{n}$ denotes an integral over phase space. @Eq:def_dos motivates the common way of specifying the DOS using the Dirac delta function
$$
g(E) = \int\limits_\Omega \mathrm d\mathbf n\ \delta(\mathcal H(\mathbf n) - E)
$$ {#eq:dirac_dos}
which matches the intuition of plotting an energy histogram for uniformly chosen states. I.e. for discrete energies $E_1 < \cdots < E_n$ (the histogram bins) and random states $\mathbf n_1,\ldots,\mathbf n_N$ we can approximate the DOS as
$$
g_\text{discrete}(E_i) = \frac1N \sum\limits^N_{j=1} \delta_{\mathcal H(\mathbf n_j), E_i}
$$ {#eq:dos_histogram}
where $\delta_{\epsilon, E_i}$ is $1$ if the energy $\epsilon$ falls inside the $i$-th histogram bin and $0$ otherwise. As the number of random states and the number of histogram bins grow towards infinity, $g_\text{discrete}$ converges to @eq:dirac_dos.

For us, the DOS is of relevance because can be used to compute the partition function
$$
Z(\beta) = \int\limits_\Omega \mathrm{d}\mathbf{n}\ e^{-\beta H(\mathbf{n})} = \int \mathrm dE\ g(E) e^{-\beta E}
$$ {#eq:partition_fn_from_dos}
and thus the free energy. In the following we will discuss the Wang and Landau algorithm to estimate the DOS and evaluate its usefulness for the computation of the marginal density of a trajectory.

<!-- There are also good algorithms available to estimate the DOS of a system such that it appears a viable option to compute the marginal density of trajectories using these estimates. -->

### Wang and Landau Algorithm

Since the state spaces $\Omega$ are usually very large, one typically resorts to Monte-Carlo methods to estimate the density of states. There one generates a sequence of states $\mathbf{n}_i$ that are approximately independent and distributed according to $\mathrm{P}(\mathcal{N})$, e.g. by using the Metropolis-Hastings algorithm. For every sampled state $\mathbf{n}_i$ we can compute the Energy $\mathcal H(\mathbf{n}_i)$ and then approximate the density of states by a histogram of the energy values as in @eq:dos_histogram. To get an accurate estimate of the density of states for energy values $E$ where $g(E)$ is very small we need a lot of iterations since we will on average pick very few samples with low probability. 

The approach of Wang and Landau @2001.Wangg8b is instead to not generate samples that are distributed according to the equilibrium distribution $\mathrm{P}(\mathbf{n})$ but to adaptively vary the sampling distribution throughout the simulation. This is done such that for every energy value $E$ approximately the same number of samples are acquired such that the DOS can be accurately estimated even in regions of low density $g(E)$.

The main idea is to perform a random walk in energy space such that on average all energy levels in a predefined interval are visited equally often. If we switch from energy space to state space this implies that the probability density for this random walk to visit a state $\mathbf n$ is proportional to the reciprocal density of states $1/g[\mathcal H(\mathbf n)]$. Of course we can't sample directly using the reciprocal DOS since it is unknown. Instead, at each step of the algorithm we slightly alter our sampling distribution until the histogram of energy values becomes _flat_. Once the histogram is flat we conclude that the adaptively altered sampling distribution represents precisely the reciprocal DOS.

To perform Wang-Landau sampling we have to define energy bins $E_1,\ldots,E_n$ that represent the range of energies that we want to compute the DOS for. We start by setting $g(E_i)=1$ for all $i=1,\ldots,n$. Then—similarly to Metropolis-Hastings sampling—we iteratively propose and selectively accept new configurations such that the transition probability from energy level $E_i\rightarrow E_j$ is
$$
p(E_i\rightarrow E_j) = \min\left( 1, \frac{g(E_i)}{g(E_j)} \right)\,.
$$
After each proposal we update a histogram of visited energies $H(E_j)\leftarrow H(E_j) + 1$ and modify the density of states at $E_j$ by a constant factor $g(E_j)\leftarrow f\ g(E_j)$. This updating of the sampling distribution during the simulation is precisely what makes the random walk non-Markovian and promises fast convergence towards the correct DOS. We start the procedure with the factor $f=\exp(1)=e$. Once the histogram $H(E)$ is sufficiently flat (we use 95% flatness) we update $f\leftarrow \sqrt{f}$ and reset the histogram to continue sampling. We continue the simulation, iteratively reducing $f$ until it reaches a small predefined threshold value which allows us to adjust the tradeoff between accuracy and simulation speed.

One important consideration is the choice of energy bins. Since we are interested in the computation of the partition function using @eq:partition_fn_from_dos, the relevant energies are those where the product $g(E)e^{-\beta E}$ is not vanishingly small. Additionally we have to take into account that the estimate for the DOS is not normalized. To be able to correctly normalize the DOS we have to ensure that the range of energy bins includes all energies where the DOS is non-vanishing.

### Applying Wang-Landau to the Computation of the Marginal Density

#### Modified DOS

When working with statistical models such as the Ising model there is often a clear concept of the _state space_ $\Omega$ and an associated volume of regions in state space such that integrals of the form $\int_\Omega\mathrm d\mathbf n\ f(\mathbf n)$ are well-defined for any $f: \Omega\rightarrow\mathbb{R}$. However in the case where the individual states $\mathbf n$ represent stochastic trajectories it is not obvious what the meaning of such an integral should be. Therefore we use a modified DOS defined for a given response trajectory $\mathbf x$ by
$$
\rho(U) = \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\,\delta(U(\mathbf s, \mathbf x) - U)
$$ {#eq:modified_dos}
with $U(\mathbf s, \mathbf x) = -\ln\mathrm P(\mathbf x|\mathbf s)$. In other words we are assigning a measure $\mu$ to our state space which is defined by the inherent probability density of the signals, such that $\int\mu(\mathrm d\mathbf s) \equiv \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)$. Thus we can express the marginal density of $\mathbf x$ analogously to @eq:partition_fn_from_dos by
$$
\mathrm P(\mathbf x) = \int\mathrm dU\ \rho(U)\,e^{-U}\,.
$$ {#eq:modified_int}

#### Modified Wang-Landau algorithm

We have to slightly adapt the Wang-Landau procedure described above so that it produces an estimate of the modified DOS. To account for the density $\mathrm P(\mathbf s)$ in @eq:modified_dos we need ensure we propose configurations, asymptotically distributed according to $\mathrm P(\mathbf s)$, which we then—in a second step—accept or reject using the inverse DOS. However we can combine both of these steps into a single one by combining a Metropolis acceptance step with the usual Wang-Landau procedure. Our algorithm therefore consists of the following steps:

1. Set all entries of the modified DOS to 1, $\rho(U_i)=1, i=1,\ldots,n$.
1. Set all entries of the histogram to 0, $H(U_i)=0, i=1,\ldots,n$.
2. Loop until $f<\epsilon$.
    1. Propose a new configuration $\mathbf s^\prime\sim T(\mathbf s\rightarrow \mathbf s^\prime)$.
    2. Let $U_i$ be the potential of state $\mathbf s$ and $U_j$ the potential of state $\mathbf s^\prime$.
    3. Accept the new configuration with probability 
    $$A(\mathbf s^\prime, \mathbf s) = \min\left[1, \frac{\mathrm P(\mathbf s^\prime)}{\mathrm P(\mathbf s)} \frac{\rho(U_i)}{\rho(U_j)} \frac{T(\mathbf s^\prime\rightarrow \mathbf s)}{T(\mathbf s\rightarrow \mathbf s^\prime)} \right]\,.
    $$ {#eq:modified_acceptance_probability}
    5. Let $U^\star$ be either $U_j$ if $\mathbf s^\prime$ was accepted or $U_i$ otherwise.
    4. Update the histogram $H(U^\star)\leftarrow H(U^\star)+1$ and the DOS $\rho(U^\star)\leftarrow f\,\rho(U^\star)$.
    5. If the histogram $H$ is flat, set $f\leftarrow\sqrt f$ and reset $H(U_i)=0, i=1,\ldots,n$.

The proposal distribution $T$ can in principle be arbitrary. The definition of the acceptance probability in @eq:modified_acceptance_probability illustrates that—once the simulation is converged—$\rho$ describes how we have to modify the state space density $\mathrm P(\mathbf s)$ such that we sample uniformly over all potentials.

#### Comparison to usual Wang-Landau algorithm

Usually the Wang-Landau algorithm is used to estimate the regular DOS $g(E)=\int\mathrm d\mathbf s\ \delta(E(\mathbf s, \mathbf x) - E)$ where $E(\mathbf s, \mathbf x)=-\ln\mathrm P(\mathbf s, \mathbf x)$. We can easily verify that algebraically this definition of $E(\mathbf s, \mathbf x)$ satisfies equation @eq:partition_fn_from_dos such that
$$
Z = \mathrm P(\mathbf x) = \int\mathrm dE\ g(E)e^{-E}\,.
$$
The Wang-Landau algorithm to compute the DOS $g(E)$ would be the same as the ones described above with the sole difference being the acceptance rate which would instead be
$$
\tilde A(\mathbf s^\prime, \mathbf s) =  \min\left[1,\frac{\rho(U_i)}{\rho(U_j)} \frac{T(\mathbf s^\prime\rightarrow \mathbf s)}{T(\mathbf s\rightarrow \mathbf s^\prime)} \right]\,.
$$ {#eq:acceptance_probability}
The difference between @eq:modified_acceptance_probability and @eq:acceptance_probability precisely reflects the fact that $\rho(E)$ and $g(E)$ are defined using different integration measures ($\rho(E)$ includes $\mathrm P(\mathbf s)$ in its integral in @eq:modified_dos while $g(E)$ does not).
However the regular version of the Wang-Landau algorithm can't be used in practice to compute the marginal density in the multivariate normal system that we are using since the regular DOS grows without bound as $E\rightarrow\infty$: In our case the joint distribution $\mathrm P(\mathbf s, \mathbf x)$ is a multivariate normal distribution in $2d$ dimensions such that for any number $\epsilon>0$ there exists a $2d$-dimensional ellipsoid that contains all points such that $E(\mathbf s, \mathbf x)\leq\epsilon$. The DOS $g(\epsilon)$ is precisely the volume of the $2d-1$ dimensional surface of this ellipsoid. As we increase $\epsilon$ the size of the ellipsoid also keeps increasing such that $\lim_{\epsilon\rightarrow\infty} g(\epsilon)=\infty$. This behavior of the regular DOS makes it impossible to use the Wang-Landau algorithm to compute $g(E)$ since by its design the algorithm merely estimates the unnormalized DOS. For this reason we have to modify the DOS as is done in @eq:modified_dos which ensures that we can normalize the result of the Wang-Landau algorithm.

#### Connection to Standard Monte-Carlo Sampling

When we perform a standard Monte-Carlo estimate of $\mathrm{P}(\mathbf{x})$ we generate independent samples $\mathbf{s}_1,\ldots,\mathbf{s}_M$, all identically distributed according to $\mathrm P(\mathbf{s})$ and then compute
$$
\hat{\mathrm{P}}(\mathbf{x}) = \frac{1}{M} \sum\limits^M_{i=1} \mathrm{P}(\mathbf x|\mathbf s_i) \,.
$$
This estimate is essentially the same as performing the integral from @eq:modified_int where the density of states $\rho(U)$ is just approximated as the histogram of the potentials $U(\mathbf{s}_1),\ldots,U(\mathbf{s}_M)$. Specifically in the limit of the width of histogram bins approaching 0, the approximate modified density of states becomes $\hat{\rho}(U)=1/M\ \sum^M_{i=1} \delta(U-U(\mathbf{s}_i))$ and therefore
$$
\int\mathrm{d}U\ \hat{\rho}(U) e^{-U} = \frac{1}{M}\int\mathrm dU \left[\sum^M_{i=1} \delta(U-U(\mathbf{s}_i)) e^{-U}\right] = \frac{1}{M} \sum\limits^M_{i=1} e^{-U(\mathbf{s}_i)} = \hat{\mathrm{P}}(\mathbf{x})\,.
$$

For the purpose of comparing estimates we can therefore associate the standard MC approach with computing the empirical histogram of potential values when sampling signals according to their marginal distribution. In the next section we will compare this empirical histogram to the DOS as computed using the Wang-Landau algorithm.

### Example Results for a Wang-Landau Simulation

![Illustrating the benefit of using the Wang-Landau algorithm for the multivariate normal system at different dimensionalities: The blue lines show the modified density of states from @eq:modified_dos, computed using a conventional MC simulation. The orange line represents the Boltzmann factor $e^{-U}$ and the green line shows the integrand of @eq:modified_int which is the product of the other two quantities. For visualization purposes, all functions where rescaled such that their integrals over the displayed interval equal 1. We see that especially at higher dimensionality there is very little overlap between the green and the blue line which leads to high inaccuracy in the computation of the marginal density. For all simulations we chose the covariance matrix using $\Delta t = 64$.](figures/normalized_densities.svg){#fig:normalized_densities}

@Fig:normalized_densities makes it clear why we expect Wang-Landau sampling to lead to a better estimate of the marginal density than the brute-force Monte-Carlo computation, especially in high-dimensional state spaces. For $d=50$ and $d=200$ most of the weight of the integral is in regions where $\rho(E)\approx 0$. In these low density regions we usually get a very inaccurate estimation of the (modified) DOS by normal MC simulations since we only very occasionally sample a relevant state. The Wang-Landau algorithm ensures that for every energy there is a consistent sampling density and we get a good estimate of the DOS even in low-density regimes.

From @fig:normalized_densities we can also estimate for which range of potentials we must compute the DOS. Since the Boltzmann weight $e^{-U}$ strongly favors low-potential configurations it is important to compute the DOS for very low potentials even if it nearly vanishes there (i.e. in regions where the blue line vanishes but the green line has relevant weight).

![Plots of estimates of the modified DOS from @eq:modified_dos compared on both linear and log scales. The blue line shows the Wang-Landau estimate while the orange line is a histogram estimate using unbiased sampling according to $\mathrm P(\mathbf s)$. We see that especially in the low-potential regime the Wang-Landau estimate is much more accurate.](figures/wl_dos.svg){#fig:wl_dos}

In @fig:wl_dos we display the estimated DOS using the Wang-Landau algorithm compared with a histogram estimate of the DOS using unbiased sampling. We see that in the highly relevant regime of low potential the Wang-Landau procedure allows us to get an accurate estimate of the DOS even though its density is as low as $e^{-45}\approx 10^{-20}$. Using @eq:modified_int we compute the marginal density to be $-664.01$ whereas the "correct" value computed analytically is $-664.24$. We thus find a relative error of $0.03\%$ in this estimate.

We have found two methods, Thermodynamic Integration and Wang-Landau sampling to show very promising results for the estimation of the marginal density. So far however, we have only tested these methods using the Gaussian approximation where it is especially easy to generate trial moves to be used in the MCMC schemes. Since in the Gaussian approximation a trajectory is described by a $d$-dimensional vector $\mathbf s = (s^1,\ldots,s^d)^T\in\mathbb{R}^d$ we can simply generate a small _displacement vector_ $\xi\in\mathbb{R}^d$ from an appropriate multivariate normal with sufficiently small covariance such that we have a proposal $\mathbf s^\prime = (s^1 + \xi^1,\ldots,s^d + \xi^d)$. Given a symmetric distribution of displacements, the proposal is naturally also symmetric. With all that said, using the Gaussian approximation we can only describe a fairly limited range of stochastic processes. Since the signal dynamics are fully characterized by their first- and second-order statistics it fails to describe non-Markovian or discrete signals. Going beyond the Gaussian approximation, we are interested in signals that are described by a general stochastic differential equation (SDE). We have already discussed how to generate trajectories from SDEs in TODO. In the next section we will therefore explore some ideas for the generation of appropriate trial moves for such general stochastic trajectories.

## Generating Proposal Trajectories from General Stochastic Dynamics

- Our goal is, given an initial trajectory $\mathbf s$, to generate a correlated trajectory $\mathbf s^\prime$, distributed according to $e^{-E(\mathbf s)}$ where $E(\mathbf s) = -\ln \mathrm P(\mathbf s) + U(\mathbf s)$
- Generate trials from P(s), then accept/reject



- Other interesting approaches: Guided proposals

## Marginalizing Out Individual Components of the Biochemical Network

So far, the computation of the mutual information is limited to the case where the signal directly interacts with the response...

## Discussion

We have shown that methods originally developed for computing ensemble averages in statistical physics can be very successfully applied to the computation of marginal densities from a known joint density. Computing the marginal densities forms the basis for many Bayesian computations (including the computation of information theoretic quantities like the mutual information) and is therefore also relevant to a wide variety of researchers in fields outside of statistical physics. Before we can make use the algorithms from statistical physics we have to map quantities like _energy_ and the _partition function_ to the corresponding probability densities. Once these terms are properly defined we suddenly have access to a wide variety of literature on statistical physics to aid us with efficient computation of the marginal entropy. Specifically we evaluated the usefulness of thermodynamic integration and the Wang-Landau algorithm for multivariate normal distributions. 

Using Markov Chain Monte Carlo sampling together with thermodynamic integration we were able to achieve very good accuracy in estimates of $\mathrm P(\mathbf x)$ with a reasonable amount of computation time. We do not find any inherent bias in the estimates as we did for the brute-force Monte-Carlo estimate. The drawback of TI is that we need to produce samples using MCMC sampling. In practice to efficiently generate samples we need to find a good proposal distribution $\mathrm T(\mathbf s \rightarrow \mathbf s^\prime)$ that yields proposals $\mathbf s^\prime$ that are _nearby_ the previous signal $\mathbf s$ such that there is on average a good chance that the proposal will be accepted according to @eq:metropolis_acceptance. On the other hand it is necessary for fast convergence that $\mathbf s^\prime$ is not _too close_ to $\mathbf s$ such that the sample chain explores a sufficiently large portion of the state space. While for our toy model it is relatively easy to come up with reasonable proposals, for stochastic trajectories this is less clear. Therefore we have proposed some ideas for possible trial moves in the space of trajectories. Even though TI seems to be a promising method for the computation of marginal densities, we also tested another technique that was originally developed in physics of condensed matter, the _Wang-Landau algorithm_.

While we can achieve a very precise estimate for the marginal density $\mathrm P(\mathbf x)$ using the estimated DOS from the Wang-Landau algorithm there remain some practical difficulties. For maximum efficiency and accuracy the different parameters affecting the procedure such as the required histogram flatness, the updating scheme of the $f$ parameter, and the choice of potential bins have to be tuned for a given problem. After some tuning of these parameters for the Gaussian system for a fixed set of covariance matrices we still find the estimation to be at least one order of magnitude slower in CPU time compared to thermodynamic integration. Therefore, at least for the Gaussian system thermodynamic integration seems to be better suited to compute the marginal density.

With that said, we do expect the Wang-Landau procedure to perform especially well in cases were there are many local minima of the potential. Here using TI we might get _"stuck"_ in a specific minimum and thus not sample all relevant states. Therefore we suggest that while TI should be the method of choice for the computation of the marginal density for high dimensional systems, in specific cases it may make sense to try other approaches that are well-established in statistical physics, such as Wang-Landau. 

Even if we can compute the mutual information more efficiently by using TI without estimating the DOS for individual responses it is worth noting that the intuitive picture behind the DOS shows very clearly why the brute-force Monte-Carlo estimates of the marginal density can converge very slowly. Particularly [@fig:normalized_densities;@fig:wl_dos] illustrate why more advanced simulation methods are unavoidable for the computation of marginal densities in high-dimensional state spaces. Looking at plots of the DOS can help to understand structural properties about the system at hand and is much easier to visualize than the high-dimensional state space.

# Conclusion

We are developing a novel approach to estimate the mutual information between a environmental signal and a cellular response, fully taking into account the time-dependent stochastic dynamics. In principle, this approach is applicable to any biochemical signaling network that can be described by a master equation. One main issue that we found is the computation of accurate estimates of the marginal density which is defined by the integral $\mathrm P(\mathbf x) = \int \mathrm d\mathbf s\,\mathrm P(\mathbf s)\,\mathrm P(\mathbf x|\mathbf s)$. By understanding the close relationship between that integral and the estimation of free energy differences, we were able to employ powerful computational methods originating in statistical physics. In that regard we have produced some promising results, yet it is clear that the work is not completed until we can demonstrate an accurate estimate for the mutual information for a simple biochemical network.

## Summary of Main Results

To quantify the fidelity of biochemical networks, we motivated the computation of the _information rate_ between a time-varying environmental signal and a cellular response. The information rate describes the (asymptotic) increase of mutual information (MI) with the duration of the signal. Using the very mild assumptions that a) the signal can be described by a stochastic differential equation, and b) the biochemical information processing network can be modeled by a master equation, we derived a general Monte Carlo procedure to compute the MI between signals and responses of arbitrary duration. It is based on our ability to generate independent sample signals and—for any given signal—the ability to generate an appropriate response. Crucially, we have shown, how for a given signal $\mathbf s$ and response $\mathbf x$ we can compute the so-called log-likelihood $\ln\mathrm P(\mathbf x|\mathbf s)$ using only terms from the master equation. One part of the Monte Carlo procedure consists in averaging these log-likelihoods for independently generated signals and responses. The second part of the computation is much harder since it consists in computing the average over the logarithm of the integrated likelihood $\ln\int \mathrm d\mathbf s\,\mathrm P(\mathbf s)\,\mathrm P(\mathbf x|\mathbf s)$ for different responses. 

An important contribution of this thesis is the careful analysis of the issues that can arise while numerically computing integrals of that form. If we merely use a standard Monte Carlo computation to evaluate the integrated likelihood we find that we consistently over-estimate the overall information rate. By approximating the signals and responses as Gaussian processes we were able to understand where this problem comes from. When randomly generating signals, as the duration grows it becomes increasingly unlikely to occasionally stumble upon a signal that contributes to the integral $\int \mathrm d\mathbf s\,\mathrm P(\mathbf s)\,\mathrm P(\mathbf x|\mathbf s)$. It is easy to imagine that a given, time-varying response $\mathbf x$ could only have arisen from a very narrow set of similar signals with the only difference between them being some small, insignificant fluctuations. In principle however, we have a huge variety of possible signal realizations that could have happened. We can understand that combing through the vast set of possible signals in search of that narrow subset of signals which contribute to the integral is very inefficient if done by brute force. It is not unlike searching for a needle in the haystack.

By having understood the problem we were able to infer what changes would lead to an improved estimate. The solution lies in the modification of the sampling strategy for the signals in such a way that we are more likely to find the strongly contributing ones. Pictorially speaking, this is perhaps comparable to the use of a metal proximity detector to help with the search of the needle. In the literature, there have been suggested multiple different but related ideas to modify sampling strategies in order to get more accurate results. While reviewing a few of the common methods we realized that many of these were inspired by techniques that originated in statistical physics for the estimation of free energy differences.

Thus, by changing the mathematical formulation of our problem we are able to see the computation of the marginal entropy akin to the computation of a free energy. Using this insight we have shown that using either of two well-known methods, _thermodynamic integration_ and _Wang-Landau sampling_ we are able to estimate the marginal likelihood much more efficiently for trajectories in the Gaussian approximation. While these ideas have not yet been tested on real stochastic trajectories, the results so far are very promising.

We propose that this thesis is an important step towards a general algorithmic framework to compute mutual information for arbitrary stochastic biological systems.

## Outlook

- areas still to explore
  - test on "real systems"
  - efficient generation of proposal trajectories
  - look beyond biochemical networks, described by a master equation?

# References {.unnumbered}
