import numpy as np
from matplotlib import pyplot
from matplotlib.colors import Normalize
from scipy import stats
from scipy.special import logsumexp
from IPython.display import Markdown
import seaborn as sns

lamda = 0.005
kappa = 0.25
rho = 0.01
mu = 0.01

def sigma_ss(rho, mu, lamda, kappa):
    return kappa / lamda

def sigma_xs(rho, mu, lamda, kappa):
    return sigma_ss(rho, mu, lamda, kappa) / (lamda + mu)

def sigma_xx(rho, mu, lamda, kappa):
    return kappa / lamda * rho / mu * (1 + rho / (lamda + mu))

def corr_ss(t, rho, mu, lamda, kappa):
    return kappa / lamda * np.exp(-np.abs(t) * lamda)

def corr_xs_pos(t, rho, mu, lamda, kappa):
    return rho * kappa / lamda / (lamda + mu) * np.exp(-lamda * t)

def corr_sx_pos(t, rho, mu, lamda, kappa):
    a = rho * kappa / lamda / (lamda - mu)
    b1 = (1 + (lamda - mu)/(lamda + mu))*np.exp(-mu * t)
    b2 = - np.exp(-lamda * np.abs(t))
    return a * (b1 + b2)

def corr_xs(t, rho, mu, lamda, kappa):
    return np.where(t >= 0, corr_xs_pos(t, rho, mu, lamda, kappa), corr_sx_pos(-t, rho, mu, lamda, kappa))

def corr_sx(t, rho, mu, lamda, kappa):
    return np.where(t >= 0, corr_sx_pos(t, rho, mu, lamda, kappa), corr_xs_pos(-t, rho, mu, lamda, kappa))

def corr_xx(t, rho, mu, lamda, kappa):
    c1 = np.exp(-mu * np.abs(t)) - np.exp(-lamda * np.abs(t))
    c2 = np.exp(-mu * np.abs(t))
    d1 = rho**2 / (lamda**2 - mu**2) * kappa / lamda
    d2 = (1+rho/(lamda+mu)) * kappa / lamda * rho / mu
    return (d1*c1 + d2*c2)

def corr_z(t, rho, mu, lamda, kappa):
    c_ss = corr_ss(t, rho, mu, lamda, kappa)
    c_sx = corr_sx(t, rho, mu, lamda, kappa)
    c_xs = corr_xs(t, rho, mu, lamda, kappa)
    c_xx = corr_xx(t, rho, mu, lamda, kappa)
    return np.block([[c_ss, c_xs], [c_sx, c_xx]])

def time_matrix(N, delta_t):
    time_stamps = np.expand_dims(np.linspace(0, (N-1)*delta_t, N), 0)
    return time_stamps - time_stamps.T

data = []
for n_dim in [10, 80, 640]:
    for dt in [1, 8, 64]:        
        t = time_matrix(n_dim, dt)
        c_z = corr_z(t, rho, mu, lamda, kappa)
        data.append({'n_dim': n_dim, 'dt': dt, 'correlation': c_z})

def draw_heatmap(*args, **kw_args):
    return sns.heatmap(args[0].values[0], **kw_args)
        
data = pd.DataFrame(data)
g = sns.FacetGrid(data, col='dt', row='n_dim', sharex=False, sharey=False, aspect=1, height=5/2.54, gridspec_kws={"wspace": 0.1, "hspace": 0.2})
g = g.map(draw_heatmap, 'correlation', cbar=False, xticklabels=False, yticklabels=False).set_titles("$N$ = {row_name} | $\Delta t$ = {col_name}").set(xlabel=None)
pyplot.savefig('matrix_plots.png', dpi=300)