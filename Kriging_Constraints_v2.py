import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import erf

# -----------------------------
# Keep ORIGINAL function names:
#   - neg_log_likelihood
#   - predictor
# And keep ORIGINAL "Option" usage:
#   - Option == 'Pred'         -> return mean prediction (scalar)
#   - Option == 'NegLogExpImp' -> return negative (EI * PoF) for minimization
# -----------------------------

def neg_log_likelihood(log_theta, X, y):
    theta = 10 ** log_theta
    p = 1.99
    n = X.shape[0]

    Psi = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            r = np.sum(theta * np.abs(X[i] - X[j]) ** p)
            Psi[i, j] = Psi[j, i] = np.exp(-r)

    # nugget
    Psi += np.eye(n) * 1e-8

    try:
        U = cholesky(Psi, lower=False)
    except np.linalg.LinAlgError:
        # keep your old behavior: return big penalty
        return 1e4, None

    one = np.ones((n, 1))
    Uy = solve_triangular(U, y, lower=False)
    U1 = solve_triangular(U, one, lower=False)

    mu = (one.T @ Uy) / (one.T @ U1)
    residual = y - mu * one
    Ur = solve_triangular(U, residual, lower=False)

    sigma_sq = (Ur.T @ Ur) / n
    neg_ln_like = -(-n / 2 * np.log(sigma_sq) - 0.5 * np.sum(np.log(np.diag(U) ** 2))).item()

    return neg_ln_like, U


# ---- internal helper: GP/Kriging mean & std at x ----
def _kriging_mean_std(x, model_info):
    """
    Returns (f_mean, f_std) at x for a Kriging model defined by model_info:
      model_info['X'], model_info['y'], model_info['Theta'], model_info['U']
    """
    X = model_info['X']
    y = model_info['y']
    theta = 10 ** np.array(model_info['Theta'])
    p = 1.99
    U = model_info['U']

    x = np.atleast_2d(x)  # ensure (1,k)
    n = X.shape[0]

    one = np.ones((n, 1))

    # GLS mean mu
    mu = (one.T @ solve_triangular(U, solve_triangular(U.T, y, lower=True))) / \
         (one.T @ solve_triangular(U, solve_triangular(U.T, one, lower=True)))

    # Process variance sigma_sqr
    r = y - one * mu
    sigma_sqr = (r.T @ solve_triangular(U, solve_triangular(U.T, r, lower=True))) / n

    # Correlation vector psi
    psi = np.array([np.exp(-np.sum(theta * np.abs(X[i] - x) ** p)) for i in range(n)]).reshape(-1, 1)

    # Predictive mean
    f = mu + psi.T @ solve_triangular(U, solve_triangular(U.T, r, lower=True))
    f_mean = float(f)

    # Predictive variance
    ssqr = sigma_sqr * (1 - psi.T @ solve_triangular(U, solve_triangular(U.T, psi, lower=True)))
    f_std = float(np.sqrt(np.abs(ssqr)) + 1e-12)

    return f_mean, f_std


# Predictor function including EI and constrained EI (EI * PoF)
def predictor(x, model_info):
    """
    Keeps your ORIGINAL interface:
      - if model_info['Option'] == 'Pred' -> returns mean prediction (scalar)
      - else computes NEGATIVE constrained EI (EI * PoF) for minimization
        if 'Constraint' is in model_info.
    Notes:
      - x is expected in the SAME space as model_info['X'] (unit or physical).
      - This version uses GP-based PoF using each constraint model's predicted std
        (instead of constant 1e-3).
    """
    # --- mean/std for objective model ---
    f_mean, f_std = _kriging_mean_std(x, model_info)

    if model_info.get('Option', 'Pred') == 'Pred':
        return f_mean

    # --- Probability of feasibility (PoF) ---
    p_feasible = 1.0
    if 'Constraint' in model_info:
        for model_con in model_info['Constraint']:
            g_mu, g_sigma = _kriging_mean_std(x, model_con)
            # P(g <= 0) = Phi((0 - mu)/sigma) with Phi via erf
            zc = (0.0 - g_mu) / (g_sigma * np.sqrt(2))
            p_single = 0.5 + 0.5 * erf(zc)
            p_feasible *= float(np.clip(p_single, 0.0, 1.0))

    # --- Expected Improvement (EI) for minimization ---
    y_train = model_info['y']
    y_best = float(np.min(y_train))

    s = f_std
    z = (y_best - f_mean) / (s * np.sqrt(2))

    Phi = 0.5 + 0.5 * erf(z)
    phi = (1.0 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((y_best - f_mean) / s) ** 2)

    EI = (y_best - f_mean) * Phi + s * phi
    EI = float(max(EI, 0.0))  # safety

    # Return negative for DE minimization
    return -(EI * p_feasible)
