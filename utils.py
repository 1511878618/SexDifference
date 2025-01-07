from scipy import stats
import numpy as np

#### cal genome control lambda ####
def calculate_lambda_from_z(z_scores):
    chisq = z_scores**2
    median_chisq = np.median(chisq)
    lambda_gc = median_chisq / 0.4549364
    return lambda_gc


def calculate_lambda_from_p(p_values):
    z_scores = np.abs(stats.norm.ppf(p_values / 2))
    return calculate_lambda_from_z(z_scores)


#### cal p-value from z-score ####
from scipy.stats import norm, chi2

def p_z_norm(est, se):
    """Convert estimate and se to Z-score and P-value."""
    try:
        Z = est / se
    except (FloatingPointError, ZeroDivisionError):
        Z = float("inf")

    P = chi2.sf(Z**2, 1, loc=0, scale=1)  # 0 if Z=inf
    return P, Z