import numpy as np
import ast
from fluent_run_con_final import run_fluent_simulation
from nastran_run_con import run_nastran_simulation
from scipy.optimize import differential_evolution
from Kriging_Constraints_v2 import neg_log_likelihood, predictor

# -----------------------------
# 1) Read existing samples log
# -----------------------------
X_list, Y_list, G_list = [], [], []

with open("samples_log.txt", "r") as f:
    for line in f:
        # x
        x_str = line.split("x = ")[1].split("],")[0] + "]"
        x = ast.literal_eval(x_str)
        X_list.append(x)

        # objective y
        y_val = float(line.split("y = ")[1].split()[0])
        Y_list.append(y_val)

        # constraints (assumes g4 exists)
        g1 = float(line.split("g1 = ")[1].split(",")[0])
        g2 = float(line.split("g2 = ")[1].split(",")[0])
        g3 = float(line.split("g3 = ")[1].split(",")[0])
        g4 = float(line.split("g4 = ")[1].split()[0])
        G_list.append([g1, g2, g3, g4])

X = np.array(X_list, dtype=float)   # (n,12) physical
Y = np.array(Y_list, dtype=float)   # (n,)
G = np.array(G_list, dtype=float)   # (n,4)

# Sign is wrong in the samples! care
Y = -Y
# column vectors for your Kriging code
y = Y.reshape(-1, 1)                # (n,1)
constraint_Y = [
    G[:, 0].reshape(-1, 1),         # g1
    G[:, 1].reshape(-1, 1),         # g2
    G[:, 2].reshape(-1, 1),         # g3
    G[:, 3].reshape(-1, 1),         # g4
]

# -----------------------------
# 2) Bounds (physical)
# -----------------------------
lb = np.array([6.5, 0.2, 0 , -4, 0.15, 0.6 , 0.1, 0.1, 0.000202, 0.000202, 0.000202, 0.000202], dtype=float)
ub = np.array([15 , 1  , 15,  4, 0.3 , 0.75, 0.3, 0.3, 0.000606, 0.000606, 0.000606, 0.000606], dtype=float)

k = X.shape[1]  # should be 12

# -----------------------------
# 3) Convert X to unit space for modeling/EI
# -----------------------------
den = (ub - lb)
if np.any(den == 0):
    raise ValueError("Some bounds have ub == lb; cannot scale to unit space.")

X_unit = (X - lb) / den

# Optional sanity check: are all existing samples within bounds?
if np.any(X_unit < -1e-8) or np.any(X_unit > 1 + 1e-8):
    bad_idx = np.where((X_unit < -1e-8) | (X_unit > 1 + 1e-8))
    raise ValueError(f"Some loaded X are outside [lb,ub]. Check samples_log.txt. Example indices: {bad_idx}")

# -----------------------------
# 4) Infill settings
# -----------------------------
n_infill = 24
bounds_log_theta = [(-3, 2)] * k

# -----------------------------
# 5) Infill loop (Kriging in unit space, CFD in physical)
# -----------------------------
for it in range(n_infill):

    # ---- Fit objective model on (X_unit, y) ----
    def nll_obj(log_theta):
        nll, _ = neg_log_likelihood(log_theta, X_unit, y)
        return nll

    res_obj = differential_evolution(nll_obj, bounds_log_theta)
    theta = res_obj.x
    _, U = neg_log_likelihood(theta, X_unit, y)

    # ---- Fit constraint models on (X_unit, g) ----
    constraint_models = []
    for g in constraint_Y:
        def nll_con(log_theta, g_local=g):
            return neg_log_likelihood(log_theta, X_unit, g_local)[0]

        res_c = differential_evolution(nll_con, bounds_log_theta)
        theta_c = res_c.x
        _, U_c = neg_log_likelihood(theta_c, X_unit, g)

        constraint_models.append({
            'X': X_unit,
            'y': g,
            'Theta': theta_c,
            'U': U_c,
            'Option': 'Pred'
        })

    # ---- Constrained EI in unit space ----
    model_info = {
        'X': X_unit,
        'y': y,
        'Theta': theta,
        'U': U,
        'Option': 'NegLogExpImp',
        'Constraint': constraint_models
    }

    # differential_evolution gives (k,) -> predictor should see (1,k)
    ei_func = lambda x: predictor(np.atleast_2d(x), model_info)

    # Optimize in unit hypercube
    res_ei = differential_evolution(ei_func, bounds=[(0.0, 1.0)] * k, maxiter=1000)
    x_new_unit = res_ei.x
    x_new = lb + x_new_unit * den

    # ---- Run expensive coupled solvers in physical space ----
    y1, g1_new, half_span, a_stall = run_fluent_simulation(x_new[0], x_new[1], x_new[2], x_new[3], it)
    y2, g2_new, g3_new, g4_new = run_nastran_simulation(x_new[4], x_new[5], x_new[6], x_new[7], x_new[8], x_new[9], x_new[10], x_new[11],half_span, a_stall)

    # Combined objective (must match what you stored in samples_log.txt)
    y_new = -(y1 * (106.2517 / (12.65 + (2 * y2))))

    # ---- Append to datasets ----
    X_unit = np.vstack([X_unit, x_new_unit])
    y = np.vstack([y, [[y_new]]])

    constraint_Y[0] = np.vstack([constraint_Y[0], [[g1_new]]])
    constraint_Y[1] = np.vstack([constraint_Y[1], [[g2_new]]])
    constraint_Y[2] = np.vstack([constraint_Y[2], [[g3_new]]])
    constraint_Y[3] = np.vstack([constraint_Y[3], [[g4_new]]])

    # (Optional) also keep physical X if you want it later
    X = np.vstack([X, x_new])

    # ---- Log ----
    with open("infill_log.txt", "a") as f:
        f.write(
            f"Infill {it + 1}, x = {x_new.tolist()}, "
            f"y1 = {y1:.6f}, y2 = {y2:.6f}, "
            f"y = {y_new:.6f}, "
            f"g1 = {g1_new:.6f}, g2 = {g2_new:.6f}, "
            f"g3 = {g3_new:.6f}, g4 = {g4_new:.6f}\n"
        )