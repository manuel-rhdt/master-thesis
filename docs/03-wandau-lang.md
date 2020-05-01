---
author:
  - Manuel Reinhardt
  - Pieter Rein ten Wolde
title: Wand and Landau Algorithm
institute: AMOLF
bibliography: ["library.bib"]
link-citations: true
linkReferences: true
autoEqnLabels: true
cref: true
---

# Estimates Using the Density of States

The goal of the Wand and Landau algorithm is to compute the density of states $\rho(E)$ for a system using an adaptive Metropolis scheme where the sampling distribution is changing throughout one simulation. We want to show that this algorithm can be used to get a better estimate of the marginal probability density for random trajectories.

In the context of statistical physics we often look at configurations of a system that can be described by a parameter vector $\mathbf{n}\in\Omega$ where $\Omega$ is the state space of the system. We can typically assign a probability (density) to each configuration. For example, by considering the canonical ensemble for a given temperature and Hamiltonian $H$ we get
$$
\mathrm{P}(\mathbf{n}) = \frac{1}{Z(\beta)} e^{-\beta H(\mathbf{n})}
$$ {#eq:canonical_probability}
with the inverse temperature $\beta$ and the partition function $Z(\beta)=\int \mathrm{d}\mathbf{n}\ e^{-\beta H(\mathbf{n})}$. The Hamiltonian assigns an energy to every state, i.e. for every state $\mathbf{n}$ we have an associated energy $H(\mathbf{n})$. To learn more about the distribution of energies in our system we can now define the _density of states_ $\rho(E)$ at a given energy $E$ as the probability density of a random^[By random state we mean a state randomly sampled from the probability density function in @eq:canonical_probability and _not_ a state with a completely random parameter vector.] state $\hat{\mathbf{n}}$ to have energy $H(\hat{\mathbf{n}}) = E$. More precisely, let $\mathcal{N}$ be a random variable distributed according to @eq:canonical_probability, then
$$
\rho(E) = \mathrm{P}\left(H(\mathcal N) = E\right)\,.
$$

Since the state spaces $\Omega$ are usually very large, one typically resorts to Monte-Carlo methods to estimate the density of states. There we generate a sequence of states $\mathbf{n}_i$ that are approximately independent and distributed according to $\mathrm{P}(\mathcal{N})$, e.g. by using the Metropolis-Hastings algorithm. For every sampled state $\mathbf{n}_i$ we can compute the Energy $H(\mathbf{n}_i)$ and then approximate the density of states by a histogram of the energy values. To get an accurate estimate of the density of states for energy values $E$ where $\rho(E)$ is very small we need a lot of simulation time since we will on average pick very few samples with low probability. 

The idea of the Wand and Landau algorithm is instead to not generate samples that are distributed according to the equilibrium distribution $\mathrm{P}(\mathbf{n})$ but to adaptively vary the sampling distribution throuhgout the simulation. (TODO: I think I understand how the algorithm works but I have not yet found the time to write it down accurately.) 

We can use the Wand and Landau algorithm to compute the marginal probability density $\mathrm{P}(\mathbf{x})$ for a given response trajectory $\mathbf{x}$. To estimate this we make use of the marginalization
$\mathrm{P}(\mathbf{x}) = \int \mathrm{d}\mathbf{s}\ \mathrm{P}(\mathbf{s})\ \mathrm P(\mathbf{x}|\mathbf{s})$
over the signals $\mathbf s$. To make the connection to the statistical physics context for the Wand and Landau algorithm notationally clear, we formally introduce the _energy_ of a signal trajectory $\mathbf s$ with respect to a given response $\mathbf x$ and define it as $E_\mathbf{x}(\mathbf s) = -\ln\mathrm{P}(\mathbf{x}|\mathbf{s})$. Since there is little potential for confusion we will just drop the index $\mathbf x$ from now on. Then we can define the density of states by analogy to be
$$
\rho(E) = \mathrm P\left(E(\mathcal S) = E\right)
$$
which allows us to express the marginal probability as
$$
\mathrm{P}(\mathbf{x}) = \int \mathrm{d}\mathbf{s}\ \mathrm{P}(\mathbf{s})\ \mathrm P(\mathbf{x}|\mathbf{s}) = \int \mathrm{d}\mathbf{s}\ \mathrm{P}(\mathbf{s})\ e^{-E(\mathbf s)} = \int\mathrm{d}E\ \rho(E) e^{-E}\,.
$$ {#eq:integrated_dos}
Therefore a viable approach to estimate the marginal entropy might be to compute the approximate density of states using the Wand and Landau algorithm and then to directly perform the integral in @eq:integrated_dos.

When we perform a standard Monte-Carlo estimate of $\mathrm{P}(\mathbf{x})$ we generate indipendent samples $\mathbf{s}_1,\ldots,\mathbf{s}_M$, all identically distributed according to $\mathrm P(\mathcal{S})$ and then compute
$$
\hat{\mathrm{P}}(\mathbf{x}) = \frac{1}{M} \sum\limits^M_{i=1} \mathrm{P}(\mathbf x|\mathbf s_i) \,.
$$
This estimate is essentially the same as performing the integral from @eq:integrated_dos where the density of states $\rho(E)$ is just approximated as the histogram of the energies $E(\mathbf{s}_1),\ldots,E(\mathbf{s}_M)$. Specifically in the limit of the width of histogram bins approaching 0, the approximate density of states becomes $\hat{\rho}(E)=1/M\ \sum^M_{i=1} \delta(E-E(\mathbf{s}_i))$ and therefore
$$
\int\mathrm{d}E\ \hat{\rho}(E) e^{-E} = \frac{1}{M}\int\mathrm dE \left[\sum^M_{i=1} \delta(E-E(\mathbf{s}_i)) e^{-E}\right] = \frac{1}{M} \sum\limits^M_{i=1} e^{-E(\mathbf{s}_i)} = \hat{\mathrm{P}}(\mathbf{x})\,.
$$

![Histograms of the negative log-likelihoods per unit time for various signals and one given response. The y-axis shows probability densities which are akin to the density of states in statistical physics. The data are from Gillespie simulations using a single random trajectory $\mathbf x$ of duration $T$ and 10000 signal trajectories $\mathbf{s}_i$, each sampled from $\mathrm{P}(\mathbf s)$.](figures/density-of-states.pdf){#fig:density-of-states}

In @fig:density-of-states we show several histograms of the density of states (however with finite bin width) for a Gillespie simulation of $M=10000$ signals. In these simulations we typically under-estimate the marginal probability density (and therefore over-estimate the marginal entropy). We propose that the bias arises because we do not have good enough sampling in the low-density regions of the density of states. From @eq:integrated_dos we see that due to the factor $e^{-E}$ especially trajectories with low energy contribute strongly to the overall estimate. This suggests that we should try to bias our sampling distribution towards trajectories with low energies. One way to do this might be the Wand and Landau algorithm.
