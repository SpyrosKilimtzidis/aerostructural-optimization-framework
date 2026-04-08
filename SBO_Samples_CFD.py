from fluent_run_con_test import run_fluent_simulation
from pyDOE import lhs
import numpy as np
from scipy.optimize import differential_evolution
from Kriging_Constraints import neg_log_likelihood,predictor

# Main SBO procedure with constraints
k = 12
n_init = 120
X = lhs(k, samples=n_init)
# Bounds
lb = np.array([6.5, 0.2, 0 , -4, 0.15, 0.6 , 0.1, 0.1, 0.000202, 0.000202, 0.000202, 0.000202])
ub = np.array([15 , 1  , 15,  4, 0.3 , 0.75, 0.3, 0.3, 0.000606, 0.000606, 0.000606, 0.000606])

# Scale to actual variable space
X_real = lb + X * (ub - lb)

y_values, g1_values, g2_values, g3_values = [], [], [], []

for i, x in enumerate(X_real):
    y1, g1, g2, g3, half_span, a_stall = run_fluent_simulation(x[0], x[1], x[2], x[3], i)
    y_values.append(y1)
    g1_values.append(g1)
    g2_values.append(g2)
    g3_values.append(g3)
    with open("samples_log.txt", "a") as f:
        f.write(f"Sample {i + 1}, x = {x.tolist()}, y1 = {y1:.6f}, g1 = {g1:.6f}, g2 = {g2:.6f}, g3 = {g3:.6f}\n")

# Ensure y is a 2D column vector
y = np.array(y_values).reshape(-1, 1)
g1 = np.array(g1_values).reshape(-1, 1)
g2 = np.array(g2_values).reshape(-1, 1)
g3 = np.array(g3_values).reshape(-1, 1)
constraint_Y = [g1, g2, g3]



