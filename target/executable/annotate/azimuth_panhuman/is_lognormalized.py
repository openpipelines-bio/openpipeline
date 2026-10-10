import numpy as np


def is_lognormalized(count_matrix, target_sum=None):
    if not target_sum:
        # If no target sum unknown, use first observation's sum as reference
        target_sum = np.asarray(np.expm1(count_matrix).sum(axis=1)).ravel()[0]
    exp_row_sums = np.expm1(count_matrix).sum(axis=1)
    # Check if exponentiated row sums match target (log-normalized data)
    # and don't overflow (exponentiated row sums become infinity for raw data)
    if np.isfinite(exp_row_sums).all() and np.allclose(exp_row_sums, target_sum):
        return True

    return False
