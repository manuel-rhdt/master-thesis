---
author:
  - Manuel Reinhardt
  - Pieter Rein ten Wolde
title: Mutual Information for Trajectories in a Gaussian Framework
bibliography: ["library.bib"]
link-citations: yes
---

# Mutual Information for Trajectories in a Gaussian Framework

## Introduction

After taking measurements of an uncertain quantity $\mathcal S$ we hope that the produced observations $\mathcal X$ can lead to a reduction of our uncertainty. We use entropies to quantify uncertainty and by rephrasing Bayes' theorem as a relation of uncertainties
$$
\mathrm H(\mathcal S|\mathcal X) = \mathrm H(\mathcal X|\mathcal S) - \mathrm H(\mathcal X) + \mathrm H(\mathcal S) \equiv \mathrm H(\mathcal S) - \mathrm I(\mathcal S, \mathcal X)
$$
we can read off that the observations $\mathcal X$ decrease our uncertainty in $\mathcal S$ on average by $\mathrm I(\mathcal S, \mathcal X) = \mathrm H(\mathcal X) - \mathrm H(\mathcal X|\mathcal S)$. We call $\mathrm I(\mathcal S, \mathcal X)$ the _mutual information between $\mathcal S$ and $\mathcal X$_, $\mathrm H(\mathcal X)$ is the _marginal entropy of $\mathcal X$_ and $\mathrm H(\mathcal X|\mathcal S)$ the _conditional entropy of $\mathcal X$ given $\mathcal S$_. Therefore to understand e.g. the quality of our measurements (or of the measurement device) it can be insightful to be able to estimate the mutual information.

When we consider the point of view of a signal-processing device (e.g. a cell) we might be interested in cases where both the quantity $\mathcal S$ and the data $\mathcal X$ change over time. We then consider the values of the random variables $\mathcal S, \mathcal X$ to be trajectories or sequences of states over time. Trajectories are usually represented as high-dimensional vectors (e.g. as a sequence of states and transition times). Our motivation is to compute the mutual information between such trajectories. To do this we intend to generate trajectories using a fully stochastic model of a biochemical network based on its master equation to compute the likelihoods of individual trajectories. In these notes however we only consider a very simple multi-dimensional Gaussian system because it allows us to test and understand the pitfalls of information estimation for high-dimensional systems with relatively low amounts of computing power. Additionally for a multivariate Gaussian system there exists a simple analytical expressions for the mutual information that we use to verify our estimates.

Hence we consider the case where the joint probability distribution $\mathrm P(\mathbf{s}, \mathbf{x})$ is given by a multivariate normal distribution
$$
\mathrm P(\mathbf{s}, \mathbf{x}) = \frac{1}{\sqrt{\left( 2\pi  \right)^{2d} \det Z}} \;\exp\left[-\frac12\ (\mathbf s^T\; \mathbf x^T)\ Z^{-1}\ \binom{\mathbf s}{\mathbf x}\right]
$$
where $\mathbf s, \mathbf x \in \mathbb R^d$ are the signal and response vectors respectively and the symmetric positive-definite covariance matrix $Z\in\mathbb R^{2d\times 2d}$ has the block form
$$
Z =  \begin{pmatrix}
C_{ss} & C_{xs} \\
C_{sx} & C_{xx}
\end{pmatrix}
$$ {#eq:corr_z}
with $C_{\alpha\beta}\in\mathbb R^{d\times d}$. For this distribution there exists a simple analytical expression to compute the mutual information [@Tostevin2010;@Shannon1948]
$$
\mathrm I(\mathcal S, \mathcal X) = \frac 12 \ln\left( \frac{\det C_{ss} \det C_{xx}}{\det Z} \right)
$$
which will be our benchmark to compare the proposed Monte-Carlo estimation procedure against. In a similar way we can also acquire analytical equations for both the marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X | \mathcal S)$.

We want to estimate the mutual information by separately computing the marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X | \mathcal S)$ from simulated data. In the present case we have the full information about our system which allows us to generate correctly distributed random observations as test data.

As it turns out the main difficulty of our estimation procedure is to get an unbiased estimate of $\mathrm H(\mathcal X)$ whereas it is much more straightforward to get a good estimate for $\mathrm H(\mathcal X | \mathcal S)$. Therefore in these notes we focus on the computation of the marginal entropy and consider the conditional entropy at the end.

## Monte-Carlo Estimate for the Marginal Entropy

We compute the marginal entropy $\mathrm H(\mathcal X)$ using Monte-Carlo sampling to evaluate the necessary integrals. First we generate a number of samples $(\mathbf x_i)_{i=1,\ldots,N_x}$ that are distributed according to the distribution of $\mathcal X$. We use these to estimate the entropy
$$
\mathrm H(\mathcal X) = -\int\mathrm d\mathbf x\ \mathrm P(\mathbf x)\ln \mathrm P(\mathbf x) \approx \frac{\sum\limits_{i=1}^{N_x} \ln\mathrm P(\mathbf x_i)}{N_x} \,.
$$ {#eq:mc_entropy}

It is important to realize that we do not actually need to know the distribution of $\mathcal X$ to do create appropriate Monte-Carlo samples. Since a stochastic model for trajectories provides us with the distributions $\mathrm P(\mathbf s)$ and $\mathrm P(\mathbf x|\mathbf s)$ we can generate samples from $\mathrm P(\mathbf x)$ by first generating a sample $\mathbf s_j$ from $\mathrm P(\mathbf s)$ and then use $P(\mathbf x|\mathbf s_j)$ to generate a sample $\mathbf x_i$.

Nonetheless we see from +@eq:mc_entropy that we _do_ have to evaluate $\mathrm P(\mathbf x_i)$ for every generated sample. However, when we want to extend the method to stochastic trajectories then the distribution of the responses $\mathrm P(\mathbf x)$ is not known a priori anymore. Therefore we choose to evaluate $\mathrm P(\mathbf x_i)$ by doing a Monte-Carlo integration using signal samples $(\mathbf s_j)_{j=1,\ldots,N_s}$ that are distributed according to $\mathrm P(\mathcal S)$:
$$
\mathrm P(\mathbf x_i) = \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\ \mathrm P(\mathbf x_i|\mathbf s) \approx \frac{\sum\limits_{j=1}^{N_s} \mathrm P(\mathbf x_i | \mathbf s_j)}{N_s} \,.
$$ {#eq:mc_marginal}
While for a low-dimensional signal space it might be feasible to instead compute the marginalization integral using direct evaluation (see [@CepedaHumerez2019]) we choose to use MC evaluation to also be able to handle high-dimensional signal spaces. This is crucial since eventually we are interested in computing the mutual information between trajectories where the dimensionality of the distributions increases with trajectory length.

We can summarize the estimation procedure for the marginal entropy using the equation
$$
\mathrm H(\mathcal X) = -\left\langle \ln \left\langle \mathrm P(\mathbf x | \mathbf s) \right\rangle_{\mathrm P(\mathbf s)} \right\rangle_{\mathrm P(\mathbf x)}
$$ {#eq:mc_entropy_notation}
where we use the notation $\langle f(x)\rangle_{g(x)}$ for the expected value of $f(x)$ when $x$ is distributed according to the probability density given by $g(x)$. Thus when thinking in mathematical terms we have the shorthand $\langle f(x)\rangle_{g(x)} \equiv\int \mathrm dx\ g(x) f(x)$. We can also easily translate this notation into a Monte-Carlo estimate, i.e. $\langle f(x)\rangle_{g(x)} = \lim\limits_{N\rightarrow\infty}\frac{\sum_{i=1}^N f(x_i)}{N}$ where $x_1, x_2,\ldots$ are independent samples of the probability distribution given by $g(x)$.

The estimates for the entropies given by +@eq:mc_entropy and +@eq:mc_marginal together allow us to compute the mutual information of $\mathcal X$ and $\mathcal S$ _in principle_. In the following we discuss for which conditions and sample sizes you can expect a good entropy estimate from this method.

### Choice of Covariance Matrices

We want to carefully choose the covariance matrices such that we can expect any sampling issues that arise in the Gaussian framework to also be present when dealing with stochastic trajectories. Therefore we chose to model a very simple gene expression model described by the reaction equations
\begin{align*}
\emptyset &\xrightarrow{\kappa} S \xrightarrow{\lambda} \emptyset \\
S &\xrightarrow{\rho} S + X \\
X &\xrightarrow{\mu} \emptyset
\end{align*}
where $X$ are particles representing the cell response and $S$ are particles that will be interpreted as the signal. We describe the signal and response trajectories as a vector of values at discrete sample times, e.g. $\mathbf s = \left(s(t_1),\ldots,s(t_d)\right)^T$. For this model we can analytically compute the correlation functions. For simplicity we assume that the system is in steady state such that the correlation functions do only depend on time differences, i.e. $C_{\alpha\beta}(t, t^\prime) = C_{\alpha\beta}(t^\prime-t)$. The correlation functions then give us the elements of the covariance matrices
$$
C_{\alpha\beta}^{ij} = C_{\alpha\beta}(t_j - t_i) = \langle\alpha(t_i)\beta(t_j)\rangle\,.
$$

![Matrix plots of the full correlation matrix $Z$ from +@eq:corr_z for different values of dimensionality $d$ and $\Delta t$. Brighter colors indicate higher matrix element values. We can clearly observe the block structure of $Z$ in every matrix plot. For every matrix plot, the element with coordinates $(i,j)$ in the top left quadrant shows the correlations $\langle s(i\Delta t) s(j \Delta t)\rangle$. In the top right quadrant we see the correlations $\langle s(i\Delta t) x(j \Delta t)\rangle$ and in the lower quadrants we see $\langle x(i\Delta t) s(j \Delta t)\rangle$ and $\langle x(i\Delta t) x(j \Delta t)\rangle$ on the left and right side respectively. The product $d\Delta t$ is the duration of the signal and response trajectories. The quantity $d\Delta t$ also serves as a rough measure of the sparsity of the correlation matrices (i.e. the fraction of matrix elements lower than some cutoff). In the plot grid we see correlation matrices with equal sparsity diagonally adjacent along lines from top right to bottom left. As we move along such a diagonal of equal sparsity and increasing dimensionality, we see correlation matrices that display the same features in a gradually more refined and smooth way.](matrix_plots.png){#fig:corr}

Using this system we have two parameters left to tune. We can freely choose the number $d$ and offsets $\Delta t$ of our time samples. The duration of the trajectories $\mathbf s$ and $\mathbf x$ is given by the product $T=d\Delta t$. In figure @fig:corr we show matrix plots of the joint covariance matrix $Z$ for different values of $d$ and $\Delta t$. We can also observe that $d$ determines the dimensionality of the problem while the product $d \Delta t$ serves as a measure for the sparsity of the correlation matrices. Note that the choice of $\Delta t$ affects how well the discretized trajectories approximate physical continuous-time trajectories. However here we are not interested in comparing our results to physical systems (as is done in [@Tostevin2010]) and therefore we can simply regard $\Delta t$ as a parameter describing the sparsity of the covariance matrices.

| $\kappa$ | $\lambda$ | $\rho$ | $\mu$ |
|:--------:|:---------:|:------:|:-----:|
|   0.25   | 0.005     | 0.01   | 0.01  |

Table: Values of the reaction constants used for all computations. {#tbl:k}

### Results

We want to use this fully Gaussian model to understand how the sample sizes of the different Monte Carlo steps affect the estimate and whether there exists a bias in the approximation. We calculate the marginal entropy as a Monte Carlo average over the logarithms of the marginal distribution densities of $N_x$ sampled responses as shown in +@eq:mc_entropy. The evaluation of the marginal density itself requires a Monte Carlo average over $N_s$ sampled signals (+@eq:mc_marginal). Hence to evaluate the marginal density we need to perform nested averaging as shown in +@eq:mc_entropy_notation. We performed this procedure for various values of $N_s$ and $N_x$ and compared the estimate with reference results using the analytical expression for the entropy of a multivariate Gaussian distribution.

Both, increase of $N_x$ and increase of $N_s$ should lead to an improved estimate of $\mathrm H(\mathcal X)$. To understand the accuracy of an estimate with a given $N_s$ and $N_x$ we repeat the estimation procedure multiple times and compute the mean and the standard deviation of the individual estimation results.

![Top: relative error $\frac{\mathrm H_\text{estimate}}{\mathrm H_\text{analytical}} - 1$ for the marginal entropy as a function of $1/N_x$. Bottom: empirical variance of ensembles of 144 estimates. The solid lines show a linear extrapolation of the data points for $N_x \rightarrow\infty$. All estimates were performed using a constant number of signal samples $N_s = 400$ and $d = 200$ dimensional covariance matrices. The linear extrapolation in the bottom plot indicates that we do predict the variance of the results to vanish in the limit of infinite sampling. This behavior is generally expected for Monte Carlo estimates. Strikingly however, we find that there is a consistent offset of the average estimate from the correct result, even in the limit $N_x \rightarrow\infty$. We see that the bias scales with the sparsity of the covariance matrices.](relative_error_responses.svg){#fig:rel_err_responses}

In +@fig:rel_err_responses we see how the relative error of our estimate varies with the number of simulated responses $N_x$. Here use the same number of signals per response $N_s$ for all estimates. While---as expected---the variance of the estimate decreases when we increase $N_x$ we find that especially for very sparse covariance matrices we consistently over-estimate the marginal entropy. Indeed, we find that the systematic bias in our results seems to be independent of $N_x$.

We found that the bias is stronger for correlation matrices with higher sparsities $d\Delta t$. Since the sparsity grows with trajectory duration we can expect an increasingly strong over-estimation for longer trajectories. The sparsity can be increased either by decreasing the time-resolution or by increasing the dimensions of the covariance matrix. To understand how these parameters relate to each other we tested how the estimation error changes when we increase the dimensionality of the correlation matrices while the sparsity remains constant.

*@fig:sparsity shows how large the estimation error for the marginal is on average for different levels of sparsity. We see that in all cases an increase of sparsity corresponds to an increase of estimation error. Additionally we find that increasing the dimensions of the covariance matrices while keeping the sparsity constant tends to slightly worsen the bias of the estimates as well. As we keep increasing the dimensionality at constant sparsity the matrices gradually become a more faithful representation of the continuous correlation functions of the system (see +@fig:corr). Extrapolating the lines in +@fig:sparsity we project that for very large covariance matrices, the sparsity is the only determining factor of the estimation bias.

![Absolute error of marginal entropy estimates for different values of the sparsity $d\Delta t$ of the correlation matrices. We see that for high dimensionality the lines of constant sparsity become increasingly flat. This indicates that for high-dimensional systems the sparsity of the covariance matrix is a good measure for the difficulty of correct estimation. We therefore claim that the bias of the entropy estimate for the Gaussian system primarily depends on the sparsity of the covariance matrix. Note that for lower numbers of dimensions the covariance matrices of along the diagonals of equal sparsity look more blocky (see +@fig:corr). That may be an indicator why the estimation error is not constant for a given sparsity at lower dimensions.](sparsity.svg){#fig:sparsity}

![Relative error $\frac{\mathrm H_\text{estimate}}{\mathrm H_\text{analytical}} - 1$ as a function of $\frac{1}{N_S}$. We can see that the relative error in the marginal entropy estimate increases with the sparsity (i.e. with trajectory duration). The linear extrapolating lines emphasize that there is a noticable but very slight decrease in error as $N_s\to\infty$. This seems puzzling since for infinite sampling we should expect the error to vanish. Apparently for high sparsity covariance matrices we need extraordinarly many signal samples to achieve unbiased estimates.](error_grid.svg){#fig:error_regression}

As a next step we investigated how changes in the sampling for the marginal density $\mathrm P(\mathbf x_i)$ affect the estimation bias. Thus in +@fig:error_regression we show how the increase of simulated signals per response $N_s$ improves the estimate of the marginal entropy. Here we again see that for high number of dimensions we over-estimate the marginal entropy. An increase of $N_s$ does lead to slightly less over-estimation but the linear extrapolation indicates that even if we choose enormously high values for $N_s$ we can not expect to reduce the bias substantially.

For a given number of samples, the fraction of the trajectory space probed by the Monte Carlo scheme decreases with the duration of the trajectories. Therefore, for a given number of Monte Carlo samples we expect the estimate to become worse for longer trajectories . This is confirmed by our results. Furthermore and more surprisingly we find that we consistently over-estimate the marginal entropy and while increasing $N_s$ _does_ reduce the bias slightly it appears to require an astronomically high sampling in signal trajectory space to reach arbitrary low errors. An increase in $N_x$ however reduces the variance of the results but does not influence the bias at all.

Thus we are lead to believe that the main difficulty in estimating the marginal entropy is the Monte-Carlo marginalization of the probability density function. To estimate $\mathrm P(\mathbf x)$ we sample signals from the marginal distribution $\mathrm P(\mathbf s)$ and average over the likelihoods $\mathrm P(\mathbf x | \mathbf s)$. However as the space of signals becomes increasingly vast for longer trajectories, it becomes more and more unlikely to sample a signal where $\mathbf x$ has a non-vanishing likelihood of occurring. Indeed, when we allow ourselves to evaluate $\mathrm P(\mathbf x)$ directly from the Gaussian PDF we find that the bias of the marginal entropy estimate does not appear.

To get a better estimate of $\mathrm P(\mathbf x)$ we decided to use _importance sampling_, i.e. to bias our sampling strategy towards signals that we expect to have a high likelihood for the given response. Since Monte-Carlo estimates depend on the distribution of the chosen samples we must correct our estimate by re-weighing the samples accordingly.

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
where in the last step we applied Bayes' rule. Since the expression above is completely independent of the chosen signal samples the result is deterministic and thus has zero variance. We also see from +@eq:opt_sampling that in practice we can't directly use $w_\text{opt}(\mathbf s) = \mathrm P(\mathbf s | \mathbf x_i)$ as our sampling distribution since the evaluation of $\mathrm P(\mathbf s | \mathbf x_i) = \frac{\mathrm P(\mathbf x_i|\mathbf s) \mathrm P(\mathbf s)}{\mathrm P(\mathbf x_i)}$ itself depends on $\mathrm P(\mathbf x_i)$ which is precisely the quantity we are interested in estimating.

![Relative error as a function of the dimensionality $d$. The solid lines show the results using non-optimized sampling while the dashed lines show the results when using a sampling distribution close to the optimal distribution $\mathrm P(\mathbf s|\mathbf x)$. We see that with optimized sampling there is no consistent over-estimation anymore. All estimated were done using $d = 200$ dimensional covariance matrices.](sampling2.svg){#fig:rel_err_opt}

Instead, we can try to obtain a sampling distribution that is as close as possible to $w_\text{opt}(\mathbf s)$. A known approach involves using random samples from $\mathrm P(\mathbf s | \mathbf x_i)$ to pick the most optimal sampling distribution from a family of candidate distributions [@Chan2012]. Generating so-called _posterior samples_ from $\mathrm P(\mathbf s | \mathbf x_i) \sim \mathrm P(\mathbf x_i|\mathbf s) \mathrm P(\mathbf s)$ is generally possible without knowledge of the normalization factor $\mathrm P(\mathbf x_i)$ e.g. by using Metropolis-Sampling [@Mueller1991;@Tierney1994]. To test within the Gaussian framework whether such an approach to importance sampling could work in principle, we generate 400 posterior samples by directly sampling from the analytically known posterior distribution $\mathrm P(\mathbf s | \mathbf x_i)$. We compute the empircial mean $\bar{\mathbf s}$ and the empirical covariance $\bar C_{\mathbf s|\mathbf x_i}$ of these samples as parameter estimates for a multivariate Gaussian $\mathcal N(\bar{\mathbf s}, \bar C_{\mathbf s|\mathbf x_i})$ and use it as an optimized sampling distribution.

In +@fig:rel_err_opt we show that using optimized sampling we can strongly reduce the systematic bias in marginal entropy estimation. As expected, importance sampling is especially useful when the sparsity is very high, i.e. the trajectories are long. It is clear that for longer trajectories we expect $\mathrm P(\mathbf s | \mathbf x_i)$ to be a much more narrow sampling distribution than $\mathrm P(\mathbf s)$ whenever the $\mathcal S$ and $\mathcal X$ are not completely independent. Consequently, it becomes more and more unlikely to obtain a sample $s^\prime$ from $\mathrm P(\mathbf s)$ such that $\mathrm P(s^\prime | \mathbf x_i) > \epsilon$ for any $\epsilon > 0$ and therefore more difficult to accurately estimate $\mathrm P(\mathbf x_i)$ using an unbiased sampling distribution.

<!-- To choose a sensible sampling distribution we use the fact that we generate the responses $\mathbf x_i$ by first picking a signal $\mathbf s_i^\star$ and subsequently picking from the likelihood distribution $\mathrm P(\mathbf x|\mathbf s_i^\star)$. As a matter of fact, $\mathbf s_i^\star$ can be regarded as a single random sample from $\mathrm P(\mathbf s | \mathbf x_i)$, i.e. the “optimal” sampling distribution. It is not unreasonable to assume in this context that a multivariate Gaussian distribution centered at $\mathbf s_i^\star$ could turn out to be a good weighing function. -->

## Estimating the Conditional Entropy

We can also estimate the _conditional entropy_ and thus the mutual information within the Gaussian Framework. We express the conditional entropy using the notation introduced above
$$
\mathrm H(\mathcal X|\mathcal S) = -\int \mathrm d\mathbf s\mathrm d\mathbf x\ \mathrm P(\mathbf s)\mathrm P(\mathbf x | \mathbf s) \ln\mathrm P(\mathbf x|\mathbf s) = -\left\langle\langle\ln\mathrm P(\mathbf x | \mathbf s)\rangle_{\mathrm P(\mathbf x | \mathbf s)} \right\rangle_{\mathrm P(\mathbf s)}
$$
to show that we require nested Monte Carlo integrations to evaluate the integral. We first generate signal samples $\mathbf s_1, \ldots, \mathbf s_{N_s}$ from the density $\mathrm P(\mathbf s)$. Let $\mathbf x_i^1,\ldots,\mathbf x_i^{N_x}$ be response samples generated from $\mathrm P(\mathbf x | \mathbf s_i)$. The Monte Carlo estimate for the conditional entropy then reads
$$
\mathrm H(\mathcal X|\mathcal S) \approx - \frac1{N_s N_x} \sum\limits_{i=1}^{N_s} \sum\limits_{j=1}^{N_x} \ln\mathrm P(\mathbf x_i^j | \mathbf s_i)\,.
$$

![Comparison of the relative error of conditional entropy estimates versus marginal error estimates. The relative errors are shown on a logarithmic scale as a function of the sparsity. We can see that the relative error for the estimate of the conditional entropy is a few orders of magnitude smaller than the estimates of the marginal entropy. All estimates were performed with $N_x=25600$ and $N_s=1000$.](conditional.svg){#fig:conditional}

For both, marginal entropy and conditional entropy we have to evaluate the likelihood $\mathrm P(\mathbf x| \mathbf s)$ a total of $N_s N_x$ times. To compare the accuracy of we performed estimates of the marginal entropy with and without optimized sampling together with estimates of the conditional entropy for $N_s = 1000$ and $N_x = 25600$. In +@fig:conditional we show the relative error of both, marginal and conditional entropy estimates as a function of the sparsity. We find that the estimate of the conditional entropy is very accurate regardless of sampling size. Even with optimized sampling the marginal entropy estimate is roughly two orders of magnitude worse than a comparable conditional entropy estimate.

We can thus conclude that the main difficulty of obtaining a good estimate for the mutual information between trajectories lies in the efficient and accurate computation of the marginal entropy. A viable approach for this seems to be to use importance sampling in the signal space in the computation of the marginal probability density. Our results indicate that such an approach could also work for trajectories generated using a fully stochastic model of a biochemical network.

## References