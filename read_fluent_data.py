import numpy as np

def read_aero_data(filename):
    """
    Reads CFD data and extracts aerodynamic cell center coordinates (with Y-negated)
    and normalized face area vectors as normal vectors.

    Parameters:
        filename (str): Path to the CFD results file.

    Returns:
        aero_points (np.ndarray): Nx3 array of cell center coordinates with Y flipped.
        aero_normals (np.ndarray): Nx3 array of normalized normal vectors.
    """
    # Load the full matrix
    data = np.loadtxt(filename, skiprows = 1)

    # Extract x, y, z cell center coordinates and flip y
    x = data[:, 1]
    y = -data[:, 2]  # Y flipped for symmetry
    z = data[:, 3]
    aero_points = np.column_stack((x, y, z))

    # Extract face area components
    fx = data[:, 10]
    fy = -data[:, 11]
    fz = data[:, 12]
    mag = data[:, 9]
    # Avoid divide by zero
    mag[mag == 0] = 1e-12
    nx = fx / mag
    ny = fy / mag
    nz = fz / mag
    aero_normals = np.column_stack((nx, ny, nz))
    # Compute static forces
    static_forces = data[:, 4]*mag
    # Compute viscous forces
    viscous_x = data[:, 6] * data[:, 9]
    viscous_y = data[:, 7] * data[:, 9]
    viscous_z = data[:, 8] * data[:, 9]
    static_forces_x = static_forces*nx
    static_forces_y = static_forces*ny
    static_forces_z = static_forces*nz
    aero_forces_x  = 1.5*(static_forces_x + viscous_x)
    aero_forces_y  = 1.5*(static_forces_y + viscous_y)
    aero_forces_z  = 1.5*(static_forces_z + viscous_z)
    aero_forces = np.column_stack((aero_forces_x, aero_forces_y, aero_forces_z))

    return aero_points, aero_forces, aero_normals
