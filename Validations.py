import numpy as np
import ast
from scipy.optimize import differential_evolution
from Kriging_Constraints_v2 import neg_log_likelihood, predictor

# ============================================================
# 1) Read ORIGINAL samples from samples_log.txt
#    These are used to train the surrogate
# ============================================================
X_list, Y_list = [], []

with open("samples_log.txt", "r") as f:
    for line in f:
        # Extract x
        x_str = line.split("x = ")[1].split("],")[0] + "]"
        x = ast.literal_eval(x_str)
        X_list.append(x)

        # Extract objective y
        y_val = float(line.split("y = ")[1].split()[0])
        Y_list.append(y_val)

X = np.array(X_list, dtype=float)   # physical design variables
Y = np.array(Y_list, dtype=float)   # objective values

# Sign correction, same as in your infill code
Y = -Y

# Column vector for Kriging
y = Y.reshape(-1, 1)

# ============================================================
# 2) Bounds (physical space)
# ============================================================
lb = np.array([6.5, 0.2, 0, -4, 0.15, 0.6, 0.1, 0.1,
               0.000202, 0.000202, 0.000202, 0.000202], dtype=float)

ub = np.array([15, 1, 15, 4, 0.3, 0.75, 0.3, 0.3,
               0.000606, 0.000606, 0.000606, 0.000606], dtype=float)

den = ub - lb
if np.any(den == 0):
    raise ValueError("Some bounds have ub == lb; cannot scale.")

# ============================================================
# 3) Convert training samples to unit space
# ============================================================
X_unit = (X - lb) / den

if np.any(X_unit < -1e-8) or np.any(X_unit > 1 + 1e-8):
    raise ValueError("Some training samples are outside bounds.")

# ============================================================
# 4) Fit objective Kriging model on samples_log.txt
# ============================================================
k = X.shape[1]
bounds_log_theta = [(-3, 2)] * k

def nll_obj(log_theta):
    nll, _ = neg_log_likelihood(log_theta, X_unit, y)
    return nll

print("Fitting Kriging model on training samples...")
res_obj = differential_evolution(nll_obj, bounds_log_theta)
theta = res_obj.x
_, U = neg_log_likelihood(theta, X_unit, y)

if U is None:
    raise RuntimeError("Cholesky factorization failed for the objective model.")

model_info = {
    'X': X_unit,
    'y': y,
    'Theta': theta,
    'U': U,
    'Option': 'Pred'
}

# ============================================================
# 5) Read INFILL points from infill_log.txt
#    These are used as test/validation points for RMSE
# ============================================================
X_infill_list, Y_infill_list = [], []

with open("infill_log.txt", "r") as f:
    for line in f:
        # Extract x
        x_str = line.split("x = ")[1].split("],")[0] + "]"
        x = ast.literal_eval(x_str)
        X_infill_list.append(x)

        # Extract objective y
        # line format contains "... y = value, g1 = ..."
        y_val = float(line.split("y = ")[1].split(",")[0])
        Y_infill_list.append(y_val)

X_infill = np.array(X_infill_list, dtype=float)            # physical infill points
y_infill = np.array(Y_infill_list, dtype=float).reshape(-1, 1)  # true objective values

# If needed, uncomment this if sign must also be corrected here
# y_infill = -y_infill

# ============================================================
# 6) Convert infill points to unit space
# ============================================================
X_infill_unit = (X_infill - lb) / den

if np.any(X_infill_unit < -1e-8) or np.any(X_infill_unit > 1 + 1e-8):
    raise ValueError("Some infill points are outside bounds.")

# ============================================================
# 7) Predict surrogate values at infill points
# ============================================================
y_pred = np.array([
    predictor(np.atleast_2d(xi), model_info)
    for xi in X_infill_unit
], dtype=float).reshape(-1, 1)

# ============================================================
# 8) Compute RMSE
# ============================================================
rmse = np.sqrt(np.mean((y_infill - y_pred) ** 2))

print("\n================ RESULTS ================\n")
print(f"Number of training samples : {X.shape[0]}")
print(f"Number of infill points    : {X_infill.shape[0]}")
print(f"RMSE on infill points      : {rmse:.6e}")

# ============================================================
# 9) Optional: print true vs predicted values
# ============================================================
print("\nTrue vs Predicted on infill points:\n")
for i in range(len(y_infill)):
    err = y_infill[i, 0] - y_pred[i, 0]
    print(
        f"Infill {i+1:2d}: "
        f"true = {y_infill[i,0]: .6e}, "
        f"pred = {y_pred[i,0]: .6e}, "
        f"error = {err: .6e}"
    )