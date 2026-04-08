import numpy as np
def createmat8(name, E1, E2, nu12, G12, G23, G13, rho):
    return (
        f'material.create( "Analysis code ID", 1, "Analysis type ID", 1, "{name}", 0, '
        f'"Date: 22-Jul-23           Time: 19:01:56", "2d Orthotropic", 5, '
        f'"Directionality", 4, "Linearity", 1, "Homogeneous", 0, "Linear Elastic", 1, '
        f'"Model Options & IDs", ["", "", "", "", ""], [0, 0, 0, 0, 0], "Active Flag", 1, '
        f'"Create", 10, "External Flag", FALSE, "Property IDs", '
        f'["Elastic Modulus 11", "Elastic Modulus 22", "Poisson Ratio 12", '
        f'"Shear Modulus 12", "Shear Modulus 23", "Shear Modulus 13", "Density"], '
        f'[2, 3, 5, 8, 9, 10, 16, 0], "Property Values", '
        f'["{E1}", "{E2}", "{nu12}", "{G12}", "{G23}", "{G13}", "{rho}", ""] )'
    )

def createmat8_fail(name, Xt, Xc, Yt, Yc, S, SB, iterm):
    return (
        f'material.create( "Analysis code ID", 1, "Analysis type ID", 1, "{name}", 0, '
        f'"Date: 22-Jul-23           Time: 19:01:56", "2d Orthotropic", 5, '
        f'"Directionality", 4, "Linearity", 0, "Homogeneous", 0, "Failure", 4, '
        f'"Model Options & IDs", ["Stress", "Tsai-Wu", "", "", ""], [6, 4, 0, 0, 0], '
        f'"Active Flag", 1, "Create", 10, "External Flag", FALSE, '
        f'"Property IDs", ["Tension Stress Limit 11", "Tension Stress Limit 22", '
        f'"Compress Stress Limit 11", "Compress Stress Limit 22", '
        f'"Shear Stress Limit", "Interaction Term", "Bonding Shear Stress Limit"], '
        f'[99, 102, 100, 103, 101, 133, 132, 0], '
        f'"Property Values", ["{Xt}", "{Yt}", "{Xc}", "{Yc}", "{S}", "{iterm}", "{SB}", ""] )'
    )

def eq_moduliv1(layup, matid, hply, E1, E2, G12, nu12):
    thetadt = np.array(layup)

    if matid == 'nmat1':
        hply_value = hply[0]
    else:
        hply_value = hply[1]

    nu21 = (E2 / E1) * nu12
    Nplies = len(thetadt)
    thetadb = thetadt[::-1]  # flipped array
    h = Nplies * hply_value

    zbar = np.array([-(h + hply_value) / 2 + i * hply_value for i in range(1, Nplies + 1)])

    denom = 1 - nu12 * nu21
    Q11 = E1 / denom
    Q12 = nu12 * E2 / denom
    Q22 = E2 / denom
    Q66 = G12
    Q = np.array([[Q11, Q12, 0],
                  [Q12, Q22, 0],
                  [0, 0, Q66]])

    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    for i in range(Nplies):
        theta = thetadb[i] * np.pi / 180
        m = np.cos(theta)
        n = np.sin(theta)
        T = np.array([
            [m ** 2, n ** 2, 2 * m * n],
            [n ** 2, m ** 2, -2 * m * n],
            [-m * n, m * n, m ** 2 - n ** 2]
        ])
        Qbar = np.linalg.inv(T) @ Q @ np.linalg.inv(T).T

        A += Qbar * hply_value
        B += Qbar * hply_value * zbar[i]
        D += Qbar * (hply_value * zbar[i] ** 2 + hply_value ** 3 / 12)

    ABD = np.block([[A, B], [B, D]])
    invABD = np.linalg.inv(ABD)
    E1meq = 1 / h * ((A[0, 0] * A[1, 1] - A[0, 1] ** 2) / A[1, 1])
    E1beq = 12 / h ** 3 * ((D[0, 0] * D[1, 1] - D[0, 1] ** 2) / D[1, 1])
    G12meq = 1 / (h * invABD[2, 2])
    G12beq = 12 / (h ** 3 * invABD[5, 5])

    E1eq = max(E1meq, E1beq)
    G12eq = max(G12meq, G12beq)

    return E1eq, G12eq

def createmat1(name, E, G, rho, Xt, Yt, S):
    # Elastic properties
    key1 = (
        f'material.create( "Analysis code ID", 1, "Analysis type ID", 1, "{name}", 0, "", '
        f'"Isotropic", 1, "Directionality", 1, "Linearity", 1, "Homogeneous", 0, '
        f'"Linear Elastic", 1, "Model Options & IDs", ["", "", "", "", ""], '
        f'[0, 0, 0, 0, 0], "Active Flag", 1, "Create", 10, "External Flag", FALSE, '
        f'"Property IDs", ["Elastic Modulus", "Shear Modulus", "Density"], '
        f'[2, 8, 16, 0], "Property Values", ["{E:.6f}", "{G:.6f}", "{rho:.6f}", ""] )'
    )

    # Failure properties
    key2 = (
        f'material.create( "Analysis code ID", 1, "Analysis type ID", 1, "{name}", 0, "", '
        f'"Isotropic", 1, "Directionality", 1, "Linearity", 0, "Homogeneous", 0, '
        f'"Failure", 4, "Model Options & IDs", ["Tsai-Wu", "", "", "", ""], '
        f'[4, 0, 0, 0, 0], "Active Flag", 1, "Create", 10, "External Flag", FALSE, '
        f'"Property IDs", ["Tension Stress Limit", "Compression Stress Limit", "Shear Stress Limit"], '
        f'[99, 100, 101, 0], "Property Values", ["{Xt:.6f}", "{Yt:.6f}", "{S:.6f}", ""] )'
    )

    return key1, key2


