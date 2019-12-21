import numpy as np


# from https://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-
# matrix.html
def logdet_symmetric(matrix):
    lower_triangular = np.linalg.cholesky(matrix)
    return 2 * np.sum(np.log(np.diag(lower_triangular)))


def mutual_information_from_matrix(c_ss, c_sx, c_xx):
    z = np.block([[c_ss, c_sx.T], [c_sx, c_xx]])

    det_c_ss = logdet_symmetric(c_ss)
    det_c_xx = logdet_symmetric(c_xx)
    det_z = logdet_symmetric(z)

    return 0.5 * (det_c_ss + det_c_xx - det_z)
