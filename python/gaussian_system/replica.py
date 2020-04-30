import numpy
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from .correlation_funcs import System, time_matrix


def log_marginal_entropy_power(sys: System, t, num_x: int, num_s: int, power: int):
    marg_x = multivariate_normal(cov=sys.corr_xx(t))
    marg_s = multivariate_normal(cov=sys.corr_ss(t))

    x_samples = marg_x.rvs(num_x).reshape((num_x, 1, 1, -1))

    s_shape = (num_x, num_s, power)
    s_samples = marg_s.rvs(s_shape).reshape(s_shape + (-1,))

    log_l = sys.log_likelihood(x_samples, s_samples, t)
    return logsumexp(numpy.sum(log_l, axis=-1), axis=1) - numpy.log(numpy.double(num_s))


def replica_estimate_sim(
    sys: System, dim: int, delta_t: float, num_x: int, num_s: int, max_n: int
):
    t = time_matrix(dim, delta_t)
    c_ss = sys.corr_ss(t)
    c_sx = sys.corr_sx(t)
    c_xs = sys.corr_xs(t)
    c_xx = sys.corr_xx(t)

    n = numpy.arange(1, max_n + 1)

    data = []
    for n in n:
        # this performs the monte carlo estimate of ln P(x)^n
        z_n = log_marginal_entropy_power(sys, t, num_x, num_s, n)
        data.append(
            pd.DataFrame(
                {
                    "n": n,
                    "log_marginal_power": z_n,
                    "num_responses": 1,
                    "num_signals": num_s * n,
                    "dim": dim,
                    "delta_t": delta_t,
                }
            )
        )

    return pd.concat(data, ignore_index=True)


def estimate_entropy_replica(data, allow_intercept=True):
    def logmeanexp(data):
        import scipy

        return scipy.special.logsumexp(data, b=1 / len(data))

    def replica_marginal_entropy(data):
        import statsmodels.formula.api as smf

        model = "log_marginal_power ~ n"
        if not allow_intercept:
            model += " - 1"
        res = smf.ols(model, data=data).fit()
        return pd.Series(
            {
                "num_signals": data.num_signals.sum(),
                "num_responses": data.num_responses.sum(),
                "intercept": res.params[0] if allow_intercept else 0.0,
                "marginal_entropy": -res.params[-1],
                "stderr": res.bse[-1],
            },
            dtype=object,
        )

    data = (
        data.groupby(["dim", "delta_t", "n"])
        .agg(
            log_marginal_power=("log_marginal_power", logmeanexp),
            num_responses=("num_responses", "sum"),
            num_signals=("num_signals", "sum"),
        )
        .reset_index(level=2)
    )
    return (
        data.groupby(["dim", "delta_t"]).apply(
            replica_marginal_entropy).reset_index(), data.reset_index()
    )
