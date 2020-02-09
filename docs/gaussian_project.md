---
author: Manuel Reinhardt
title: Notes on Estimating Mutual Information from Experimental Data without Binning
---

# The model-independent approach to entropy estimation

In the common cases where people want to estimate entropies the only data available is a list of observations for the quantity of interest. The hope is that these observations are <span class='check'>approximately independent</span> samples from an underlying probability distribution that governs the observations. We assume there are $K$ different possible observations (bins) and $N$ observations. Then the naïve maximum likelihood estimator for the entropy is
$$
S^{\mathrm{ML}} = -\sum\limits_{i=1}^K \frac{n_i}{N} \ln \frac{n_i}{N}
$$
where $n_i \in \{0,\ldots, N\}$ is the number of observations in bin $i$. While the maximum likelihood estimate of the entropy is not unbiased there have been shown to exist ways to reduce or eliminate the bias even in the deeply undersampled regime [@Paninski2003;@Nemenman2004].

There are two problems with this approach. One is that these estimates typically do not take into account any problem-specific prior knowledge that may be available to improve the estimate of the entropy. However in practice this is usually not a big issue since even very weakly informative priors can lead to sensible estimates for the entropy [@Nemenman2004]. A more fundamental issue arises if the dimensionality of the observed data is very large. This happens naturally if we want to e.g. analyze the entropy for a system where every observation is an entire trajectory $(t_i, x_i)_{i=1,\ldots,L}$. As $L$ gets larger we need absurd amounts of samples to be able to define the bins in a meaningful way. We have to conclude that in this regime it is not feasible anymore to estimate entropy in a model-independent way. Therefore we will show how we can combine knowledge of the concrete system of interest with Monte-Carlo sampling to estimate the entropy and the mutual information of high-dimensional data. There is a comparison of a related technique and other approaches to estimate the mutual information for time-varying signals in [@CepedaHumerez2019].

# Entropy of the Output of a Noisy Channel

One of the main advances of Shannon's information theory is the realization that the efficiency of an information transmission channel crucially depends on the statistical properties of the input signal. Therefore it makes little sense to study a lossy channel without any model of the signal.

From this point of view we model the system of interest as a random variable $\mathcal{S}$ representing the “signal” and a random variable $\mathcal{X}$ representing the “response” or output of the channel. We assume that we have a model (i.e. a probability distribution) for the signal but that we are only able to observe the response. Additionally we need to have a model that describes the noise introduced by the channel. Such a model defines the likelihood $\mathrm P(\mathcal X | \mathcal{S}=s)$ of the responses for any given signal $s$.

In a typical system we expect that by knowing the outcome of $\mathcal{X}$ we gain _additional_ knowledge of the outcome of $\mathcal{S}$. Note that we already have some information about $\mathcal{S}$, namely we know the marginal probability distribution $\mathrm P(\mathcal{S})$ _a priori_. All the additional information that we gain by having observed the outcome $x$ of $\mathcal{X}$ must therefore be embedded in $\mathrm P(\mathcal S | \mathcal{X} = x)$. It follows that if we wish to gain any knowledge on the signal, $\mathcal{X}$ and $\mathcal{S}$ must not be independent. If they _were_ independent then $\mathrm P(\mathcal S) \equiv \mathrm P(\mathcal S | \mathcal{X} = x)$, i.e. the marginal distribution $\mathrm P(\mathcal S)$ contains all obtainable information on the signal regardless of the outcome of $\mathcal{X}$. The _amount_ of information about the signal that we gain by observing the response is quantified by the _mutual information_ $\mathrm{I}(\mathcal{X}, \mathcal{S}) = \mathrm H(\mathcal X) - \mathrm H(\mathcal X | \mathcal S)$.

We can calculate the marginal distribution of the responses by marginalizing out the signals
$$
\mathrm P(x) = \int \mathrm ds\; \mathrm P(x, s) = \int \mathrm ds\; \mathrm P(s)\,\mathrm P(x|s)
$$
and use that result to compute the (differential) entropy of $\mathcal X$
\begin{align}
\mathrm H(\mathcal X) &= - \int \mathrm dx\; \mathrm P(x)\,\ln\mathrm P(x)\\
&= - \int \mathrm dx\int \mathrm ds\; \mathrm P(s)\,\mathrm P(x|s)\,\ln\left[\int \mathrm ds^\prime\; \mathrm P(s^\prime)\,\mathrm P(x|s^\prime)\right]
\label{eq:marginal}
\end{align}
using only the distributions we assume to know a priori. In the same fashion we can also compute the conditional entropy
$$
\mathrm H(\mathcal X | \mathcal S) = - \int \mathrm ds\; \mathrm P(s)\int \mathrm dx\,\mathrm P(x|s)\,\ln\mathrm P(x|s)
$$
using only the quantities contained in our model.

## Relation to Bayesian Updating

We can also reframe this same problem as finding the entropy of measured data while having a prior belief for the statistics of the measured quantity. In typical bayesian reasoning we have a prior belief for the distribution of a quantity of interest $\mathrm P_\text{prior}(s)$. We then collect data $x$ to improve our estimate for the distribution of $\mathcal S$. The optimal (or _rational_) way of updating our expectations from the observed data is given by Bayes' rule
$$
\mathrm P_\text{posterior}(s|x) = \frac{\mathrm P (x | s)}{\mathrm P (x)} \mathrm P_\text{prior}(s)\,.
$$ {#eq:bayes_update}
Note that while the probability of the data $\mathrm P(x)$ is not explicitly conditioned on any variable, it does depend on how we choose our prior belief for $\mathcal S$ (since that determines what probability we assign to the measured data)! This becomes more clear if we rewrite +@eq:bayes_update as
$$
\mathrm P_\text{posterior}(s|x) = \frac{\mathrm P (x | s)}{\int\mathrm ds^\prime\ \mathrm P_\text{prior}(s^\prime)\mathrm P(x|s^\prime)} \mathrm P_\text{prior}(s)\,.
$$
To understand how much information we gain from the observed data, we can also rephrase +@eq:bayes_update in terms of marginal and conditional entropies
$$
\mathrm H(\mathcal S | \mathcal X) = \mathrm H(\mathcal X | \mathcal S) - \mathrm H(\mathcal X) + \mathrm H(\mathcal S) = \mathrm H(\mathcal S) - \mathrm I(\mathcal X, \mathcal S) \,.
$${#eq:bayes_entropies}

The interpretation of +@eq:bayes_entropies is that _on average_ our uncertainty of $\mathcal S$ is reduced by $\mathrm I(\mathcal X, \mathcal S)$ if we optimally update our beliefs given the observed data. Importantly the equation does not imply that every possible observation reduces our uncertainty of $\mathcal S$ by the same amount! <span class='check'>Indeed, the more surprised we are by the observed data (given our prior belief) the more drastic should be the update to our belief.</span> However before having obtained any data we will _expect_ the experimental results to reduce our uncertainty by roughly $\mathrm I(\mathcal X, \mathcal S)$.

---

Now imagine a situation where you are presented with a number of experimental observations of some quantity and are asked to estimate the entropy of the acquired samples. In other words you are tasked with estimating $\mathrm H(\mathcal X)$ where $\mathcal X$ is the random variable from which the experimental observations are drawn. The discussion above makes it clear that you need _some_ prior belief about the statistics of the underlying quantity $\mathcal S$ that the experimentalist tried to measure. Perhaps surprisingly, you find that your estimate of $\mathrm H(\mathcal X)$ _depends_ on what your prior beliefs are. This is a manifestation of the fact that entropy is a measure of uncertainty which necessarily depends on what your preexisting knowledge is.

From +@eq:bayes_entropies we can read off
$$
\mathrm H(\mathcal X) = \mathrm H(\mathcal X | \mathcal S) + \mathrm I(\mathcal X, \mathcal S)\,.
$${#eq:mutual_information}

So the entropy of the data is always _at least_ the amount of information that we gain by observing it. In this context the conditional entropy $\mathrm H(\mathcal X | \mathcal S)$ describes the quality of the measurement apparatus: The lower the conditional entropy is, the more accurately the experimental observations model the physical quantity $\mathcal S$. An experimenter that is very confident in the quality of their measurements would therefore be able to approximate the information gain of their data by estimating its entropy. However this does not free them from considering what their prior beliefs are.

The imaginary situation described above can be compared to the reality of a living cell. Since cells are interested in making correct predictions about the current (and future) state of their environment we would expect that their signalling networks evolved in order to perform bayesian inference. Hence by analyzing the cell's biochemical signalling network (obtaining the corresponding $\mathrm P(x | s)$) and the typical environmental conditions ($\mathrm P(s)$) we should in principle be able to estimate the amount of information $\mathrm I(\mathcal X, \mathcal S)$ that the cell gains by observing its environment.

# Estimating the Mutual Information from Sampled Data

The previous discussions serve to understand which quantities are required to calculate the mutual information. We need

- a prior belief for the statistics of the signal, i.e. a prior $\mathrm P(s)$ and
- a model for the measurement process (a noisy channel), i.e. the likelihood $\mathrm P(x|s)$.

To show that we can in principle correctly estimate the mutual information we first analyze a simple toy problem. We consider the case where the joint probability distribution $\mathrm P(\mathbf{s}, \mathbf{x})$ is given by a multivariate normal distribution
$$
\mathrm P(\mathbf{s}, \mathbf{x}) = \frac{1}{\sqrt{\left( 2\pi  \right)^{2d} \det Z}} \;\exp\left[-\frac12\ (\mathbf s^T\; \mathbf x^T)\ Z^{-1}\ \binom{\mathbf s}{\mathbf x}\right]
$$
where $\mathbf s, \mathbf x \in \mathbb R^d$ are the signal and response vectors respectively and the symmetric covariance matrix $Z\in\mathbb R^{2d\times 2d}$ has the block form
$$
Z =  \begin{pmatrix}
C_{ss} & C_{xs} \\
C_{sx} & C_{xx}
\end{pmatrix}
$$
with $C_{ij}\in\mathbb R^{d\times d}$. For this distribution there exists a simple analytical expression to compute the mutual information [@Tostevin2010]
$$
\mathrm I(\mathcal S, \mathcal X) = \frac 12 \ln\left( \frac{\det C_{ss} \det C_{xx}}{\det Z} \right)
$$
which will be our benchmark to compare our statistical analysis against.

We want to estimate the mutual information using +@eq:mutual_information, i.e. by separately computing the marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X | \mathcal S)$ from observed data. In the present case we have the full information about our system which allows us to generate correctly distributed random observations as test data.

## Marginal Entropy

We compute the marginal entropy using Monte-Carlo sampling to evaluate the integrals in +@eq:marginal. To do that we first generate a number of samples $(\mathbf x_i)_{i=1,\ldots,N_x}$ that are distributed according to the distribution of $\mathcal X$ and then make the estimate
$$
\mathrm H(\mathcal X) = -\int\mathrm d\mathbf x\ \mathrm P(\mathbf x)\ln \mathrm P(\mathbf x) \approx \frac{\sum\limits_{i=1}^{N_x} \ln\mathrm P(\mathbf x_i)}{N_x} \,.
$${#eq:mc_entropy}

It is important to realize that we do not actually need to know the distribution of $\mathcal X$ to do create appropriate Monte-Carlo samples. Since our model provides us with the distributions $\mathrm P(\mathbf s)$ and $\mathrm P(\mathbf x|\mathbf s)$ we can generate samples from $\mathrm P(\mathbf x)$ by first generating a sample $\mathbf s_j$ from $\mathrm P(\mathbf s)$ and then use $P(\mathbf x|\mathbf s_j)$ to generate a sample $\mathbf x_i$.

Nonetheless we see from +@eq:mc_entropy that we _do_ have to evaluate $\mathrm P(\mathbf x_i)$ for every generated sample. We are only allowed to use probabilities that are part of our prior beliefs. Therefore we can only evaluate $\mathrm P(\mathbf x_i)$ by doing a Monte-Carlo using signal samples $(\mathbf s_j)_{j=1,\ldots,N_s}$ that are distributed according to $\mathrm P(\mathcal S)$:
$$
\mathrm P(\mathbf x_i) = \int\mathrm d\mathbf s\ \mathrm P(\mathbf s)\ \mathrm P(\mathbf x|\mathbf s) \approx \frac{\sum\limits_{j=1}^{N_s} \mathrm P(\mathbf x_i | \mathbf s_j)}{N_s} \,.
$${#eq:mc_marginal}

Eq. @eq:mc_entropy and +@eq:mc_marginal together allow us to find a good estimate for the marginal entropy of $\mathcal X$ _in principle_. In the following we show for which conditions and sample sizes you can expect a good entropy estimate from this method.

## Conditional Entropy

_TODO_

<!-- blabla



$$
\mathrm P(\mathbf{x}) = \frac{1}{\sqrt{\left( 2\pi  \right)^d \det C_{xx}}} \;\exp\left[-\frac12\ \mathbf x^T\ C_{xx}^{-1}\ \mathbf x\right]
\,.
$$
Note that in this case it is n -->

# References

<div id="refs"></div>
