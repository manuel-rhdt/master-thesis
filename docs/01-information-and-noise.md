---
bibliography: ["library.bib"]
autoEqnLabels: true
link-citations: true
linkReferences: true
cref: true
---

# Acknowledgements {.unnumbered}

I would like to thank...

# Introduction

Information processing...

## Goal of the Thesis

## Related Work

## Structure

# Information and Noise in Biological Systems

_Noise_ is inherent across diverse biological systems and remains relevant at all biological scales. From _stochastic gene expression_ and random _action potential spikes_ in neuronal networks at the cellular scale to the _development of multicellular organisms_ all the way to the _variations in population level_ of competing species in whole ecosystems, we find examples of processes which can only be described precisely by taking into account noise as an intrinsic feature [@2002:Elowitz;@2008:Faisal;@1990:Parsons;@2011:Hallatschek;@2014:Tsimring]. In this thesis we focus on the smaller end of this scale, namely on the stochastic description of _biochemical networks_. These comprise among others _gene expression_, _gene regulatory networks_, and _cell signaling_ networks, all of which exhibit noise due to small copy numbers of participating components. Additionally, since biological processes often happen out of equilibrium even macroscopic quantities can exhibit large fluctuations. The main source of noise at the cellular level may be fluctuations in _gene expression_ which propagate to higher levels of biological organization @2014:Tsimring. The abundance of noise in all these systems invites the questions of how cells can reliably make correct decisions, even in complex and changing environments and how cells are able to _encode_ the information about their environment using biochemical networks.

To successfully function, cells must generally correctly react to environmental changes. This requires processing the environmental cues they receive through their receptors and thereby filtering out the useful information from the noisy signal. The processing of information through biochemical networks can be quite elaborate with individual network components being able to perform analogous functions as silicon computational devices such as transistors @2020:Leifer. It is thus tempting to think that optimization of information processing drives the evolution of cellular signaling networks. To understand information processing from the cell's point of view we employ the very general framework of _information theory_ which was originally developed to address questions of the reliability of telecommunications @1948:Shannon. Information theory has been successfully used in biology to study cellular communication, embryonic development and other biological systems [@2009:Tkačik;@2020:Uda] [citation needed]. Notably, there are many parallels between biological signal processing and _noisy channels_ which are used for example to describe information transmission across a telephone line.

A crucial feature of _information theory_ is that its results are broadly applicable, irrespective of the nature of the communication channel or the medium used to transmit a signal. The communication channel merely describes any kind of abstract device that processes an _input_ in a probabilistic way to yield a corresponding _output_. It turns out that the study of information transmission through such a channel is closely related to the study of noise since the amount of noise introduced by a communication channel sets an upper bound to the amount of data that can be transmitted through it, known as the _channel capacity_ @2006:Cover. Consequently, the output can be described as a deterministic, lossless transformation of the input _plus_ some random noise from the channel which leads to a loss of information. Note that in biological systems the signal itself is typically a fluctuating quantity such that the noise in the output is a combination of the channel noise and the signal noise. Since in cell signaling both, input and output are time-varying quantities we require a description of our system that allows for deterministic _and_ stochastic time evolution.

_Differential equations_ are generally extremely useful to describe any kind of system that evolves deterministically. Therefore it is natural to try to extend the framework of _ordinary differential equations_ (ODEs) to also include the ability to describe the effects of noise. Historically, this approach to modeling stochastic dynamics has first been formulated heuristically by Langevin to describe _Brownian motion_ @1908:Langevin. Later the theory of _stochastic differential equations_ (SDEs) was put on solid mathematical footing by Itô and Stratonovich through the development of _stochastic calculus_ which has been successfully used for applications in physics, biology, economics and others [@2010:Kunita;@1997:Bunkin] [citation needed]. The solutions to SDEs are not ordinary functions like for ODEs but _stochastic processes_ that describe the probabilities for the system to be in any state for every instant in time. Consequently, a stochastic process contains the probabilities for any possible individual sequence of states in time, i.e. the probabilities for individual _trajectories_. Since SDEs contain a complete account of noise in the system, information theoretic concepts like the _entropy_ and the _mutual information_---which we are going to use to understand information transmission in cell signaling---can be applied to stochastic processes. While SDEs can be formulated to describe the evolution of biochemical networks in a discrete state space it is generally more useful to use a less general but simpler _chemical master equation_ for these kinds of problems @2009:Gardiner. 

The chemical master equation is a description of a subset of stochastic processes by deriving the _time-evolution_ of the probability distribution over the discrete state space. I.e. instead of describing the stochastic change to an individual state at a given time it focuses on the _deterministic_ evolution of the whole probability distribution over all states. Conveniently, for a given set of chemical reactions that form a reaction network, we can easily find the corresponding chemical master equation that describes the stochastic dynamics of this network given some assumptions of homogeneity. The stochastic process that emerges of this formulation describes the probabilities for the individual counts of all species and how these probabilities change with time. The ease with which the chemical master equation allows the construction of a stochastic process for any kind of biochemical netork makes the idea very attractive to try to use _master equations_ as the basis for information theoretic computations. If we can _solve_ the master equation we in principle have access to all stochastic (and therefore information theoretic) properties of the corresponding biochemical network. E.g. in @2010:Tostevin it is shown how by analytically solving some very simple biochemical networks (using some approximations) it is possible to compute many quantities related to information transmission for time-varying signals and corresponding responses of these networks.

In most cases however, chemical master equations cannot be solved analytically and instead require computing averages for ensembles of _stochastic trajectories_. For instance, in Shannon's information theory the amount of communicated information is not a function of individual signal-response pairs but an averaged quantity that depends on the probability of seeing _any_ random signal-response-pair. Hence computationally efficient generation of stochastic realizations for a given master equation is a central requirement for the exact computation of information processing in chemical networks. A very well known algorithm for the _exact_ generation of trajectories from a given initial condition is the _stochastic simulation algorithm_ (SSA) also known by the name of its inventor as the _Gillespie algorithm_ @1976:Gillespie. The most widely used variant is the _direct Gillespie method_ which works by alternatingly a) computing the time during which no reactions happen and b) choosing which of the available reactions to perform next. As a result we generate a list of times where some reaction happens and a corresponding list of reactions that specifies the exact trajectory that was generated. This algorithm works quite well in practice and is also used for the work presented in this thesis. It is still worth mentioning that for systems that evolve at many different time scales simultaneously, the direct Gillespie method can be computationally inefficient since by its design it always operates at the smallest time scale. Therefore there have been developed further trajectory-generation algorithms that can generate _approximately_ correct trajectories by accumulating various reactions into a single time step such as the $\tau$-leap method @2001:Gillespie.

The paragraphs above provide an introductory overview of techniques and issues related to information processing in biological contexts. In this thesis we aim to present computational advances that allow the efficient estimation of the amount of information that a cell can use from its environment. In the following sections of this chapter we explain in more detail how we can model biochemical signaling networks in cells and how to understand information transmisssion in that context. The following chapters...

## Information Theory in the context of cellular signaling networks

In this section we motivate the use of a quantities known as the _mutual information_ and the _information rate_ as relevant properties to understand optimality criteria for biochemical signaling networks.

### Mutual Information as an Efficiency Measure in Cell signaling

![Abstracting cell signaling as an information channel. The channel's input is an environmental signal that the cell needs to respond to. The signal processing happens through a biochemical network which "computes" a response which is the output of the information channel. The mutual information between signals and responses quantifies the cell's ability to discern between different signals and choose appropriate responses.](figures/information_cartoon.svg){#fig:information_cartoon}

In a general sense, cells sense chemical _signals_ from their environment e.g. through receptors on their membrane. These signals provide the cell with important information e.g. about current environmental conditions, their position inside a structure or the location of food. To translate the signal into a useful response (such as expressing a certain gene or changing movement direction) cells have evolved biochemical signaling networks that recognize and process the signals. We use a general description of the cell that is depicted in @fig:information_cartoon where the signal acts as the input of the information channel. The processing of the signal that yields a response is assumed to be a known biochemical network and the response is one species of the biochemical network that acts as the "result" of the computation and represents the reaction of the cell to the signal.

For any given signal there are many stochastically possible responses. Conversely, for any given response there is a range of signals that could have produced it. Both of these statements of uncertainty can be quantified using a single concept: the _mutual information_. In information theoretic terms we quantify the _uncertainty_ of a random variable $\mathcal S$ by the _entropy_ $$
\mathrm H(\mathcal S)=-\int\limits_{\sigma(\mathcal S)} 
\mathrm d\mathbf s\ \mathrm P(\mathbf s)\,\ln\mathrm P(\mathbf s)
$$ {#eq:signal_entropy}
where $\sigma(\mathcal S)$ is the set of possible realizations of $\mathcal S$ and we use $\mathrm P(\mathbf s)$ to denote the probability (density) of $\mathbf s$ with respect to the distribution of $\mathcal S$.  If $\mathcal S$ describes the signal then a large entropy signifies that there is a large range of possible signals that could be expected by the cell. Now given the response $\mathbf x$ to a signal $\mathbf s$, we expect $\mathbf x$ to contain information about $\mathbf s$ such that the uncertainty about the signal is reduced. The conditional entropy $\mathrm H(\mathcal S|\mathcal X)$ denotes the _average_ remaining uncertainty of a signal after observing the response, hence it reads
$$
\mathrm H(\mathcal S|\mathcal X)=-
\int\limits_{\sigma(\mathcal X)} 
\mathrm d\mathbf x\ \mathrm P(\mathbf x)
\int\limits_{\sigma(\mathcal S)} 
\mathrm d\mathbf s\ \mathrm P(\mathbf s|\mathbf x)\ln\mathrm P(\mathbf s|\mathbf x) \,.
$$ {#eq:conditional_entropy}
$\mathrm P(\mathbf s|\mathbf x)$ is the conditional distribution of the signals for a given response $\mathbf x$ which encodes the transmission characteristics of the communication channel. So we have @eq:signal_entropy which describes the distribution of signals and @eq:conditional_entropy which describes how information is processed. Combining [@eq:signal_entropy;@eq:conditional_entropy] we can express the _average_ amount of information gained on the signal by observing the response
$$
\mathrm I(\mathcal S,\mathcal X) = \mathrm H(\mathcal S) - \mathrm H(\mathcal S|\mathcal X) = 
\int\limits_{\sigma(\mathcal X)} 
\mathrm d\mathbf x\ \mathrm P(\mathbf x)
\int\limits_{\sigma(\mathcal S)} 
\mathrm d\mathbf s\ \mathrm P(\mathbf s|\mathbf x)
\ln \frac{\mathrm P(\mathbf s|\mathbf x)}{\mathrm P(\mathbf s)}
$$ {#eq:mi_form1}
which is precisely the _mutual information between $\mathcal S$ and $\mathcal X$_. It depends on both, characteristics of the communication channel and the statistics of the input signal. Notably the $\mathrm I(\mathcal S,\mathcal X)$ is symmetric under exchange of $\mathcal S$ and $\mathcal X$ such that we can express it as
$$
\mathrm I(\mathcal S,\mathcal X) = \mathrm H(\mathcal X) - \mathrm H(\mathcal X|\mathcal S) = 
\int\limits_{\sigma(\mathcal S)} 
\mathrm d\mathbf s\ \mathrm P(\mathbf s)
\int\limits_{\sigma(\mathcal X)} 
\mathrm d\mathbf x\ \mathrm P(\mathbf x|\mathbf s)
\ln \frac{\mathrm P(\mathbf x|\mathbf s)}{\mathrm P(\mathbf x)}\,,
$$ {#eq:mi_form2}
resulting in a more useful formula for the Monte-Carlo estimation of the MI. Therefore, on average, a response reduces the uncertainty about the signal by exactly the same amount that a signal reduces the uncertainty about the possible responses.

<!-- The mutual information has been used previously to understand TODO elaborate -->

### Information Transmission for Time-Varying Signals

Biochemical networks may store information about the signal in the time-dependency of the response, for example in cellular Ca^2+^ signaling information seems to be encoded in the timing and duration of Calcium bursts [@2010:Tostevin;@2008:Boulware;@2020:Richards]. For the case where the signal can be regarded as slowly changing with respect to the response, Cepeda-Humerez, et. al. propose a Monte-Carlo technique for the estimation of the MI that includes information that is stored in the full temporal dynamics of the response @2019:Cepeda-Humerez. As described in their work, that method is limited to situations that are well-described by a finite number of discrete signals.

Often however, biochemical networks not only respond to instantaneous signal levels but also to changes in the signal over time. [reference needed]. Therefore, we build on the technique in @2019:Cepeda-Humerez by extending it to allow for time-varying, stochastic signals as well. In this way we aim to find a novel way to compute the MI for time-varying signals _and_ responses for general biochemical networks.

The study of time-varying quantities motivates the use of the _information rate_ which is the asymptotic rate at which the MI between signal and response increases @2010:Tostevin
$$
\mathrm I_R = \lim\limits_{T\rightarrow\infty} \frac{\mathrm I(\mathcal S_T,\mathcal X_T)}{T}
$$ {#eq:information_rate}
where $\mathcal S_T$ and $\mathcal X_T$ are random variables over _trajectories_ of length $T$. By _trajectories_ we denote an entire time-trace of the signal or response, instead of singular values at specific times. Since @eq:information_rate describes the information gained by the cell in a unit time interval it may be an important quantity for the cell to optimize for.

Hence we find that the fundamental building block for the computation of the information rate at the cellular level is to estimate the mutual information between trajectories of finite length. In the remainder of this thesis, we describe methods to compute the mutual information between trajectories $\mathcal S_T,\mathcal X_T$ based on a stochastic model of the biochemical signaling network.

<!-- ### Mutual Information between Stochastic Trajectories

A trajectory $X$ with $N$ steps is defined by a set of pairs $X=\{(t_i, \mathbf{x}_i)\; |\; i=0\ldots N-1 \}$ where $\mathbf{x}_i$ defines the trajectory value at time $t_i$. We can also have random variables over trajectories and therefore probability distributions over the space of all trajectories.

As a next step we can make sense of the entropy of a trajectory. Let $\mathcal{X}_N$ be a random variable over trajectories of length $N$. We call

$$
\mathrm H(\mathcal{X}_N) = - \int\limits_{X\in \sigma(\mathcal{X}_N)} dX\; \mathrm{P}(\mathcal{X}_N = X)\; \ln \mathrm{P} (\mathcal{X}_N = X)
$$ {#eq:entropy_integral}

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

These averages are defined as integrals over the very high-dimensional space of trajectories and thus very hard to evaluate analytically or numerically in the general case. Our goal is use *Monte-Carlo sampling* in the trajectory space to evaluate the above averages. To do this we have to sample trajectories from their probably distribution and we need to evaluate the likelihood for a response given a signal. -->

## Stochastic Modeling of Biochemical Networks

- There exist model-free methods to quantify the entropies associated with signals and responses
- These can only yield an upper bound for the entropy
- In practice, biochemical networks may operate far from that bound(?)

### Markov Processes

- memoryless

### Chemical Master Equation

As a model for the biochemical processing that takes place inside a cell we suppose that all interactions can be described by a chemical networks composed of different molecular species and reactions between them. Such networks can be described by a _chemical master equation_ which makes it possible to compute all the probabilities associated with the time-evolution of such a system.

For illustrative purposes, let's consider a highly simplified model of gene expression consisting of two components and four reactions 
$$ \begin{gathered}
\emptyset \xrightarrow{\kappa} S \xrightarrow{\lambda} \emptyset\\
S \xrightarrow{\rho} S + X\\
X \xrightarrow{\mu}\emptyset\,.
\end{gathered} $$ {#eq:reaction_network1}
Here $S$ and $X$ are arbitrary chemical species whose particle counts we want to describe stochastically. The constants $\kappa, \lambda, \rho, \mu$ determine the rates at which the individual reactions occur. Hence, for every signal particle there is a constant rate $\rho$ to be sensed by the cell which triggers the creation of an $X$. Additionally, $X$ particles decays by themselves over time with a per-particle rate of $\mu$. Assuming a well stirred system in thermal equilibrium, it can be shown that the probabilities for the individual reactions happening at least once in the time interval $[t, t+\mathrm\delta t]$ are
$$
\begin{aligned}
p^{(\kappa)}_{[t, t+\mathrm\delta t]} &= \kappa\delta t + \mathcal{O}(\delta t^2)\\
p^{(\lambda)}_{[t, t+\mathrm\delta t]}(s) &= s\lambda\delta t + \mathcal{O}(\delta t^2)\\
p^{(\rho)}_{[t, t+\mathrm\delta t]}(s) &= s\rho\delta t + \mathcal{O}(\delta t^2)\\
p^{(\mu)}_{[t, t+\mathrm\delta t]}(x) &= x\mu\delta t + \mathcal{O}(\delta t^2)
\end{aligned}
$$ {#eq:transition_probabilities}
where $s$ and $x$ denote the particle numbers of the respective species at time $t$. Consequently, the probability for _any_ of the reactions to occur at least once in the time interval $[t, t+\mathrm\delta t]$ is
$$p_{[t, t+\mathrm\delta t]}(s, x) = (\kappa + s\lambda + s\rho + x\mu)\ \mathrm \delta t + \mathcal{O}(\delta t^2)$$ {#eq:exit_probability}

Using these expressions we can write down the so-called _chemical master equation_ for this network. Let $\mathrm P_{s,x}(t)$ be the probability that the system is in state $(s, x)$ at time $t$. Assuming that at most one reaction happens in the small time interval $[t, t+\delta t]$ we can use the transition probabilities from [@eq:transition_probabilities;@eq:exit_probability] to write
$$
\begin{aligned}
\mathrm P_{s,x}(t + \delta t) =& \phantom{+}p^{(\kappa)}_{[t, t+\mathrm\delta t]}\ \mathrm P_{s-1,x}(t)\\ &+
p^{(\lambda)}_{[t, t+\mathrm\delta t]}(s + 1)\ \mathrm P_{s+1,x}(t)\\ &+
p^{(\rho)}_{[t, t+\mathrm\delta t]}(s)\ \mathrm P_{s,x-1}(t)\\ &+
p^{(\mu)}_{[t, t+\mathrm\delta t]}(x + 1)\ \mathrm P_{s, x + 1}(t)\\ &+ 
\left[1 - p_{[t, t+\mathrm\delta t]}(s, x)\right]\ \mathrm P_{s,x}(t)
\end{aligned}
$$
and by taking the limit $\delta t\rightarrow 0$ we arrive at the chemical master equation
$$
\begin{aligned}
\frac{\partial \mathrm P_{s,x}(t)}{\partial t} &= \lim\limits_{\delta t\rightarrow 0} \frac{\mathrm P_{s,x}(t + \delta t) -  \mathrm P_{s,x}(t)}{\delta t}\\
&= \kappa\ \mathrm P_{s-1,x}(t) +
(s+1)\lambda\ \mathrm P_{s+1,x}(t) +
s\rho\ \mathrm P_{s,x-1}(t) +
(x+1)\mu\ \mathrm P_{s, x + 1}(t)\\ &\phantom{=} - 
(\kappa + s\lambda + s\rho + x\mu)\ \mathrm P_{s,x}(t)
\end{aligned}
$$ {#eq:chemical_master_equation}
In an analogous way a chemical master equation can be derived for any biochemical network [@2009:Gardiner] and thus forms the basis for our further computations. Stochastic processes that can be described by a master equation are generally called _jump processes_ and provide a useful abstraction for the processes that we want to consider.

### Jump Processes

Since particle counts can't ever become negative, @eq:chemical_master_equation describes a Markov process in continuous time with the state space $\{(s, x) | s\in\mathbb{N}_0, x\in\mathbb{N}_0\}$. In general, every continuous-time Markov process with a discrete state space obeys a master equation. Such processes are also commonly called _jump processes_ since they generate discontinuous sample paths [@2017:Weber].

A jump process with state space $\mathcal{U}$ and an initial state $\mathbf{x}_0\in\mathcal{U}$ at time $t_0$ generates trajectories that can be described by a sequence of pairs $(\mathbf{x}_i, t_i)_{i=1,2,\ldots}$ where at every _transition time_ $t_i$ there occurs a jump in state space $\mathbf{x}_{i-1}\rightarrow \mathbf{x}_{i}$. As illustrated in @eq:chemical_master_equation the master equation encodes the transition rates for all possible state changes in the system. Similarly, for an arbitrary jump process we can formulate the master equation
$$
\frac{\partial \mathrm P(\mathbf x, t|\mathbf x_0, t_0)}{\partial t} = \sum\limits_{\mathbf x^\prime\in\mathcal U} \left[
w_t(\mathbf x, \mathbf x^\prime)\ \mathrm P(\mathbf x^\prime, t|\mathbf x_0, t_0)
- w_t(\mathbf x^\prime, \mathbf x)\ \mathrm P(\mathbf x, t|\mathbf x_0, t_0)
\right]
$$ {#eq:general_master_eq}
where $w_t(\mathbf x^\prime, \mathbf x)$ specifies the rate for the transition $\mathbf x \rightarrow \mathbf x^\prime$ at time $t$ and $w_t(\mathbf x, \mathbf x) = 0$. The first term of the sum in @eq:general_master_eq is known as the _gain_ since it describes how probability flows from other states into the current one while the second term is called _loss_ as it expresses how much probability flows away from the current state. This form of the master equation allows us to find a relatively simple expression for the probability of individual trajectories of jump processes.

### Probability Densities of Trajectories

From [@eq:mi_form1;@eq:mi_form2] we can see immediately that the computation of the mutual information between two random variables relies on the calculation of (conditional) probability densities for these variables. Our interest in this thesis lies in the computation of the mutual information for time-varying signals and responses. Hence we need to evaluate the mutual information between random variables whose individual samples are _entire trajectories_. To do this we require the notion of probability for a given trajectory, based on a given biochemical model for a cell. In the following derivation we show how the master equation as formulated in @eq:general_master_eq makes it possible to compute the trajectory probability exactly for the corresponding reaction network. In the next chapter we will then discuss how to use these results in practice to compute the conditional probability of a response $\mathrm P(\mathbf x|\mathbf s)$ for a biochemical network and a known signal.

By the probability of a trajectory with $N$ jumps we denote joint probability density $\mathrm P(\mathbf x_0, t_0;\ldots;\mathbf x_N,t_N)$. We include the initial state $\mathbf x_0, t_0$ in the probability since the initial condition itself is usually given by a probability distribution. Using conditional probabilities we can write
$$
\begin{aligned}
\mathrm P(\mathbf x_0, t_0;\ldots;\mathbf x_N,t_N) =\; &\mathrm P(\mathbf x_0,t_0)\,\mathrm P(\mathbf x_1,t_1|\mathbf x_0,t_0)\,\mathrm P(\mathbf x_2,t_2|\mathbf x_1,t_1;\mathbf x_0,t_0)\cdots\\
&\mathrm P(\mathbf x_N,t_N|\mathbf x_{N-1},t_{N-1};\ldots;\mathbf x_0,t_0)\,.
\end{aligned}
$$ {#eq:trajectory_probability}
We can make use of the fact that jump processes are Markov processes to simplify @eq:trajectory_probability such that every conditional probability only explicitly depends on the immediately preceding state
$$
\begin{aligned}
\mathrm P(\mathbf x_0, t_0;\ldots;\mathbf x_N,t_N) &= \mathrm P(\mathbf x_0,t_0)\,\mathrm P(\mathbf x_1,t_1|\mathbf x_0,t_0)\,\mathrm P(\mathbf x_2,t_2|\mathbf x_1,t_1)\cdots \mathrm P(\mathbf x_N,t_N|\mathbf x_{N-1},t_{N-1})\\
&= \mathrm P(\mathbf x_0,t_0)\,\prod\limits^{N}_{i=1} \mathrm P(\mathbf x_i,t_i|\mathbf x_{i-1},t_{i-1})
\,.
\end{aligned}
$$ {#eq:trajectory_probability_product}
Hence the expression for the probability of a trajectory contains two distinct kinds of probability distributions, a) the probability of the initial condition $\mathrm P(\mathbf x_0,t_0)$ and b) the transition probabilities for every step in the trajectory given by $\mathrm P(\mathbf x_i,t_i|\mathbf x_{i-1},t_{i-1})$. The initial condition $\mathrm P(\mathbf x_0,t_0)$ depends on the specific problem and will often be taken to be the steady-state distribution. For the transition probabilities however, we can find a simple expression in terms of the master equation.

The transition probability $\mathrm P(\mathbf x_i,t_i|\mathbf x_{i-1},t_{i-1})$ describes the probability for a small segment of the trajectory from $t_{i-1}$ to $t_i$ where first, the system is at state $\mathbf x_{i-1}$ for a duration $t_i - t_{i-1}$ and then makes a jump $\mathbf x_{i-1}\rightarrow \mathbf x_{i}$. We are now going to derive how to express the probabilities of the constant part as well as the probability of the jump using only terms from the master equation given in @eq:general_master_eq. At the start of the segment, the probability for there to be no jump until at least $t > t_{i-1}$ is called the _survival probability_ $S(\mathbf x_{i-1}, t_{i-1}, t)$. By noticing that the change in the survival probability is given by the _loss_ term of the master equation
$$
S(\mathbf x_{i-1}, t_{i-1}, t+\delta t) = S(\mathbf x_{i-1}, t_{i-1}, t) \left[ 1-\delta t
\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t_{i-1}}(\mathbf x^\prime, \mathbf x_{i-1}) + \mathcal O(\delta t^2) \right]
$$
we can motivate the differential equation
$$
\frac{\partial S(\mathbf x_{i-1}, t_{i-1}, t)}{\partial t} = - S(\mathbf x_{i-1}, t_{i-1}, t) \sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t_{i-1}}(\mathbf x^\prime, \mathbf x_{i-1})
$$
with the solution
$$
S(\mathbf x_{i-1}, t_{i-1}, t) = \exp\left( -\int\limits^{t}_{t_{i-1}}\mathrm dt^\prime 
\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t^\prime}(\mathbf x^\prime, \mathbf x_{i-1})
\right)\,.
$$ {#eq:survival_probability}
The expression $S(\mathbf x_{i-1}, t_{i-1}, t)$ is the probability for no reaction happening from time $t_{i-1}$ until _at least_ $t$. The survival probability therefore is the cumulative probability distribution for the _waiting time_ $\tau_{i-1} = t - t_{i-1}$, i.e. using standard notation we can write $\mathrm P(\tau_{i-1} \geq \delta) = S(\mathbf x_{i-1}, t_{i-1}, t_{i-1}+\delta)$. Since the probability for the waiting time to be _exactly_ $\tau_{i-1}=t_i-t_{i-1}$ is _zero_ we instead compute the probability density
$$\mathrm P(\tau_{i-1} = t_i-t_{i-1}) =
\left. 
\frac{\partial\mathrm P(\tau_{i-1} < \delta)}{\partial\delta} 
\right|_{\delta = t_i-t_{i-1}} = 
\left. 
\frac{\partial [1 - S(\mathbf x_{i-1}, t_{i-1}, t_{i-1} + \delta)]}{\partial\delta}
\right|_{\delta = t_i-t_{i-1}}\,.
$$ {#eq:survival_probability_density}

The jump $\mathbf x_{i-1}\rightarrow \mathbf x_{i}$ at time $t_i$ itself happens with probability
$$\mathrm P(\mathbf x_{i-1}\rightarrow \mathbf x_{i}, t_i) = \frac{w_t(\mathbf x_{i}, \mathbf x_{i-1})}{\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t}(\mathbf x^\prime, \mathbf x_{i-1})}\,.$$ {#eq:jump_probability}
Combining the jump probability with the survival probability we can thus express the transition probability density by the multiplication
$$
\mathrm P(\mathbf x_i,t_i|\mathbf x_{i-1},t_{i-1}) = 
\mathrm P(\tau_{i-1} = t_i-t_{i-1})\;
\mathrm P(\mathbf x_{i-1}\rightarrow \mathbf x_{i}, t_i)
$$
and by inserting the results from [@eq:survival_probability;@eq:survival_probability_density;@eq:jump_probability] we arrive at the expression
$$
\mathrm P(\mathbf x_i,t_i|\mathbf x_{i-1},t_{i-1}) = w_t(\mathbf x_{i}, \mathbf x_{i-1})
\exp\left( -\int\limits^{t_i}_{t_{i-1}}\mathrm dt^\prime 
\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t^\prime}(\mathbf x^\prime, \mathbf x_{i-1})
\right)
\,.
$$ {#eq:transition_probability}

Therefore @eq:transition_probability plugged into @eq:trajectory_probability_product allows the exact computation of probability densities for arbitrary stochastic trajectories. The results were derived in enough generality such that time-dependent transition rates are explicitly allowed. We will see that the analytical expression for the probability of a stochastic trajectory is central for the estimation of the mutual information throughout this thesis.  [@Eq:trajectory_probability_product;@eq:transition_probability] will be used to compute the conditional probability $\mathrm P(\mathbf x|\mathbf s)$ for stochastically generated signal and response trajectories $\mathbf s$ and $\mathbf x$. Furthermore, the survival probability and the jump probability introduced in this section form the basis of the _Stochastic Simulation Algorithm_ that is introduced in the following subsection.

### Generating Stochastic Trajectories for Jump Processes

One main ingredient for a computational recipe to calculate the mutual information is the formula for the probability of a trajectory, as derived in the previous section. Of equal importance is the correct generation of stochastic realizations for the biochemical model we intend to study. So far we showed how to model any biochemical network as a _jump process_ through the formulation of the chemical master equation. In the previous section we were able to derive an equation for the probability of a trajectory using only terms from the master equation. Similarly, we can formulate the _stochastic simulation algorithm_ which efficiently generates trajectories for any jump process with a known master equation. In the following we explain the variant named _direct Gillespie alogithm_ after its inventor @1976:Gillespie.

The direct Gillespie algorithm is an efficient and exact procedure to generate stochastic trajectories for jump processes. Starting from an initial state $\mathbf x_0$ at time $t_0$ it generates a trajectory $\mathbf x = (\mathbf x_0, t_0;\ldots;\mathbf x_N, t_N)$ that is a realization of the corresponding stochastic process. Since the trajectories of a jump process are constant segments, separated by discontinuous jumps, at every point in time the system either remains unchanged or it performs an instantaneous jump. Therefore the algorithm works in two phases that are repeated over and over until a trajectory of the desired length is generated.

In the first phase we compute the stochastic waiting time $\tau$ until the next _event_ occurs. Given that we are currently at time $t_{i-1}$ and state $\mathbf x_{i-1}$, the probability for the waiting time being $\tau\geq t-t_{i-1}$ is given by the survival probability from @eq:survival_probability. By the inverse function method we can generate a correctly distributed waiting time $\tau$ using a random variate $u$ that is uniformly distributed between 0 and 1 and then solving the equation
$$
u = 1-\exp\left( -\int\limits^{t_{i-1} + \tau}_{t_{i-1}}\mathrm dt^\prime 
\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t^\prime}(\mathbf x^\prime, \mathbf x_{i-1}) \right)
$$ {#eq:inverse_function_method}
for $\tau$. We then set $t_i = t_{i-1} + \tau$. For general transition rates $w_{t^\prime}$ @eq:inverse_function_method will have no analytical solution and it may be necessary to use approximate numerical integration techniques to solve for $\tau$. If however the transition rates are time-independent for a given state we find the simple solution
$$
\tau = -\frac{\ln(1-u)}{\sum\limits_{\mathbf x^\prime\in\mathcal U} w_{t_{i-1}}(\mathbf x^\prime, \mathbf x_{i-1})}\,.
$$

After the first phase we already know the time of the next event. Therefore all that remains to do is to decide which reaction should happen. The probability for the individual jumps is given by @eq:jump_probability and thus we can make a weighted choice between all reactions according to their probabilities. Using this choice we have found the next state $\mathbf x_i$ and we have finished one step in the stochastic simulation.

Using both of 

<!-- ### Simulating a Biochemical Network Driven by an External Signal

TODO
If the signal itself is described by a _jump process_ we can even generate _exact_ response trajectories
TODO -->

## A Simple Model of Gene Expression

In the course of this thesis we explore computational methods for the estimation of the mutual information in cellular reaction networks. To test and analyse our algorithms we perform computations for a very simple biochemical network for gene expression, consisting in merely four reactions and two species. It represents one of the simplest possible problems that a Monte-Carlo approach should be able to solve. Consequently, an estimation procedure that fails even on this simple model will almost certainly not provide satisfactory results for more realistic and complex networks. So on one hand we use a very simple biochemical network as a minimal system that an algorithm must be able to solve. On the other hand this model for gene expression has already been studied from an information theoretical point of view. Its simple structure allowed Tostevin, et. al. @2010:Tostevin to derive precise analytical approximations for the mutual information that we can use to test the results of our computational estimates against. In the following we describe the biochemical network and derive some useful analytical results.

We begin by specifying the 4 reactions that make up a simple model for gene expression
$$
\begin{gathered}
\emptyset \xrightarrow{\kappa} S \xrightarrow{\lambda} \emptyset\\
S \xrightarrow{\rho} S + X\\
X \xrightarrow{\mu}\emptyset
\end{gathered}
$$ {#eq:simple_reaction_network}
and note that it is the same network that was already used in @eq:reaction_network1. It includes two species, $S$ and $X$ that represent the signal and the response, respectively. The signal dynamics are fully described by the first two reactions, i.e. a birth-death process where signals are "born" stochastically with a constant rate $\kappa$ and every signal particle has a constant decay rate of $\lambda$. The response dynamics are given by the other two reactions in @eq:simple_reaction_network. Each signal particle can stochastically produce a response with rate $\rho$ and every response itself decays with rate $\mu$. The stochastic process that corresponds to the full reaction network under well-mixed conditions is described by a chemical master equation that we anticipatingly already derived for this network in @eq:chemical_master_equation. While the stochastic formulation describes the full reaction dynamics (assuming the system is well-mixed) we can understand many of its properties from a deterministic description.

The average number of signal particles $s(t)$ is described by the deterministic rate equation
$$
\partial_t s(t) = \kappa - \lambda\ s(t)\,.
$$ {#eq:det_s}
This ODE has a fixed point that represents the steady-state average of the signal $\bar s = \kappa / \lambda$. In a similar way we can deterministically describe the average number of response particles $x(t)$ using the ODE
$$
\partial_t x(t) = \rho\ s(t) - \mu\ x(t)\,.
$$ {#eq:det_x}
Analogous to the signal we have a steady state average for the response $\bar x = \rho\bar s / \mu$. In this case, the deterministic equations are not only useful to derive steady-state averages but can also be used to understand _correlations_ between the components of the system. If $S_t$ is the stochastic count of signal particles at time $t$ we call $C_{ss}(t, t^\prime) = \langle S_t S_{t^\prime}\rangle$ the _autocorrelation function_ of $\mathcal S$. Similarly we can also define correlation functions between the signal and the response, such as $C_{sx}(t, t^\prime) = \langle S_t X_{t^\prime}\rangle$, where by $X_t$ we denote the stochastic number of X molecules at time $t$. Thus, in total we have four different correlation functions $C_{ss}, C_{sx}, C_{xs},$ and $C_{xx}$. If we assume the system is in steady state, the correlation functions do only depend on the time-difference $t^\prime - t$ such that we can write $C_{\alpha\beta}(t, t^\prime) = C_{\alpha\beta}(t^\prime - t)$ @2009:Gardiner. For simple biochemical networks like @eq:simple_reaction_network it is relatively straightforward to analytically derive the correlation functions for the steady state.

Since the deterministic [@eq:det_s;@eq:det_x] are both _first-order_ ODEs we say that the biochemical system has _linear equations of motions_ such that for some matrix $A(t)$ and vector $\xi(t)$ we can write
$$
\frac{\partial \mathbf z(t)}{\partial t} = A(t)\ \mathbf z(t) + \xi(t)
$$
with $\mathbf z(t) = (s(t),x(t))^T$. We say that the corresponding Markov process describes a _linear system_. In the following we make use of the fact that for linear systems we can apply the _regression  theorem_ @2009:Gardiner $$\partial_tC_{ij}(t) = -\sum_k A_{ik}(t)C_{kj}(t)$$ where $C_{ij}$ are the correlation functions. Specifically, for our biochemical network we find the following set of coupled differential equations for the correlation functions
$$
\begin{aligned}
\partial_t C_{ss}(t) &= -\lambda\ C_{ss}(t) \\
\partial_t C_{sx}(t) &= \rho\ C_{ss}(t)-\mu\ C_{sx}(t) \\
\partial_t C_{xs}(t) &= -\lambda\ C_{xs}(t) \\
\partial_t C_{xx}(t) &= \rho\ C_{xs}(t)-\mu\ C_{xx}(t)
\end{aligned}
$$
with the solutions for $t\geq0$ (assuming steady-state initial conditions)
$$
\begin{aligned}
C_{ss}(t) &= \frac\kappa\lambda e^{-\lambda t} \\
C_{sx}(t) &= \frac{\rho\kappa}{\lambda(\lambda - \mu)} 
\left[ \left( 1 + \frac{\lambda - \mu}{\lambda + \mu} \right) e^{-\mu t} - e^{-\lambda t} \right] \\
C_{xs}(t) &= \frac{\rho\kappa}{\lambda(\lambda + \mu)} e^{-\lambda t}  \\
C_{xx}(t) &= \frac{\rho^2 \kappa}{\lambda(\lambda^2 - \mu^2)} \left( e^{-\mu t} - e^{-\lambda t} \right) + \left( 1 + \frac{\rho}{\lambda + \mu} \right) \frac{\rho \kappa}{\lambda \mu} e^{-\mu t}\,.
\end{aligned}
$$ {#eq:correlation_functions}
Since by definition $C_{\alpha\beta}(-t) = C_{\beta\alpha}(t)$, from @eq:correlation_functions we know the complete correlation functions for the system in steady state.

The steady-state averages and the correlation functions represent the first and second moments of the trajectories of $\mathcal S$ and $\mathcal X$. By discarding all higher-order moments and discretizing time we build an approximate stochastic model for the chemical reaction network that allows further analytic results. Since we will approximate the trajectories using samples form multivariate normal distributions, this method is often called the _Gaussian approximation_. For a given discretization $t_1<t_2<\cdots<t_d$ we describe the signal and response trajectories as a vector of values at discrete sample times, e.g. $\mathbf s = \left(s(t_1),\ldots,s(t_d)\right)^T$. In the following we define a multivariate probability density $\mathrm P(\mathbf{s}, \mathbf{x})$ for the signal and response trajectories such that the resulting approximate system has identical correlation functions to the full biochemical network.

Hence we consider the case where the joint probability distribution $\mathrm P(\mathbf{s}, \mathbf{x})$ is given by a multivariate normal distribution
$$
\mathrm P(\mathbf{s}, \mathbf{x}) = \frac{1}{\sqrt{\left( 2\pi  \right)^{2d} \det Z}} \;\exp\left[-\frac12\ (\mathbf s^T\; \mathbf x^T)\ Z^{-1}\ \binom{\mathbf s}{\mathbf x}\right]
$$ {#eq:joint_multivariate}
where $\mathbf s, \mathbf x \in \mathbb R^d$ are the signal and response vectors respectively and the symmetric positive-definite covariance matrix $Z\in\mathbb R^{2d\times 2d}$ has the block form
$$
Z =  \begin{pmatrix}
C_{ss} & C_{xs} \\
C_{sx} & C_{xx}
\end{pmatrix}
$$ {#eq:corr_z}
with matrices $C_{\alpha\beta}\in\mathbb R^{d\times d}$. The correlation functions in @eq:correlation_functions then give us the elements of the matrix blocks
$$
C_{\alpha\beta}^{ij} = C_{\alpha\beta}(t_j - t_i)\,.
$$
For the joint distribution in @eq:joint_multivariate there exists a simple analytical expression to compute the mutual information between $\mathcal S$ and $\mathcal X$ [@2010:Tostevin;@1948:Shannon]
$$
\mathrm I(\mathcal S, \mathcal X) = \frac 12 \ln\left( \frac{\det C_{ss} \det C_{xx}}{\det Z} \right)
$$ {#eq:analytical_mi}
which will be our benchmark to compare the proposed Monte-Carlo estimation procedure against. In a similar way we can also acquire analytical equations for both the marginal entropy $\mathrm H(\mathcal X)$ and the conditional entropy $\mathrm H(\mathcal X | \mathcal S)$.

While we can use @eq:analytical_mi to directly evaluate the mutual information for trajectories using discretized time, Tostevin, et. al. @2010:Tostevin derive—specifically for this reaction network—a formula to analytically compute the mutual information rate (in nats per unit time) in the limit of infinitely fine discretization
$$
I_R = \frac{\lambda}{2}\left( \sqrt{1 + \frac{\rho}{\lambda}} -1 \right)\,.
$$ {#eq:analytical_rate}

In this subsection we analyzed a simple, yet non-trivial example of a biochemical network for which we intend to compute the mutual information using Monte-Carlo techniques. Since we were able to derive the correlation functions for this network, we can make use of the Gaussian approximation to estimate the mutual information using @eq:analytical_mi. That approximation will turn out to be very useful to understand the accuracy of Monte-Carlo estimates. To make the different results in this thesis comparable we used consistent parameter values for the reaction constants that are shown in @tbl:k.

| $\kappa$ | $\lambda$ | $\rho$ | $\mu$ | $\bar s$ | $\bar x$ | $\tau_s$ | $\tau_x$ | $I_R$ |
|:--------:|:---------:|:------:|:-----:|:--------:|:--------:|:--------:|:--------:|:-----:|
|   0.25   | 0.005     | 0.01   | 0.01  | 50       | 50       | 200      | 180      | 0.00183 |

Table: Values of the reaction constants along with the steady-state averages, correlation times and approximate information rate from @eq:analytical_rate. The values for the reaction constants were used for all computations unless stated otherwise. {#tbl:k}

## References


