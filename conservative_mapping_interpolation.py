import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from findelementcenters import findelementcenters2


def conservative_mapping_projection(aero_points, aero_forces, aero_normals, struct_bdf):
    """
    Conservative CFD->STRUCT mapping (force-conserving) using kNN + inverse-distance weights.
    Keeps the SAME inputs/outputs as your original function.

    Inputs (unchanged):
        aero_points, aero_forces, aero_normals, struct_bdf

    Outputs (unchanged):
        node_forces, grd
    """
    # Load structural element centers and areas
    centerpoints, struct_areas, quads, _, grd = findelementcenters2(struct_bdf)

    Nnodes = len(grd)
    struct_points = centerpoints[:, :3]  # Structural element centers
    struct_areas = np.asarray(struct_areas).ravel()

    Na = aero_points.shape[0]
    Ns = struct_points.shape[0]

    # ---- Basic checks
    if aero_points.ndim != 2 or aero_points.shape[1] != 3:
        raise ValueError("aero_points must be (Na,3).")
    if aero_forces.shape != aero_points.shape:
        raise ValueError("aero_forces must have same shape as aero_points (Na,3).")
    if struct_points.ndim != 2 or struct_points.shape[1] != 3:
        raise ValueError("struct_points must be (Ns,3).")
    if struct_areas.shape[0] != Ns:
        raise ValueError("struct_areas length must match number of structural centroids.")

    # ---- Diagnostics (safe to remove)
    span = np.linalg.norm(struct_points.max(axis=0) - struct_points.min(axis=0))
    if span > 0:
        dmin = cdist(aero_points, struct_points).min()
        if dmin > 0.2 * span:
            print("WARNING: aero_points and struct_points seem far apart (frame mismatch likely).")

    # ---- Conservative kNN mapping (normals not needed; kept only to preserve signature)
    k_neigh = min(16, Ns)   # you can tune 8–32
    p = 2.0                 # inverse-distance power
    eps = 1e-12

    tree = cKDTree(struct_points)
    dist, idx = tree.query(aero_points, k=k_neigh)

    # Ensure (Na,k)
    if dist.ndim == 1:
        dist = dist[:, None]
        idx = idx[:, None]

    # Weights: area / (r^p)
    w = struct_areas[idx] / (np.power(dist, p) + eps)  # (Na,k)

    # Normalize weights per aero point => each row sums to 1 (force-conservative)
    w_sum = np.sum(w, axis=1, keepdims=True)
    bad = (w_sum[:, 0] <= 0.0) | (~np.isfinite(w_sum[:, 0]))
    if np.any(bad):
        w[bad, :] = 1.0
        w_sum[bad, :] = np.sum(w[bad, :], axis=1, keepdims=True)
    w = w / w_sum

    # Accumulate mapped element forces
    struct_forces = np.zeros((Ns, 3), dtype=float)
    for j in range(idx.shape[1]):
        np.add.at(struct_forces, idx[:, j], aero_forces * w[:, j:j+1])

    # NOTE: No correction_factor needed; global force is already conserved.

    # ---- Element -> node distribution (QUAD4)
    Nfem = quads[:, 0:4].astype(int)

    # Create mapping from node ID to index
    node_map_dict = {nid: idxn for idxn, nid in enumerate(grd[:, 3].astype(int))}
    node_map = np.vectorize(lambda nid: node_map_dict.get(nid, -1))(Nfem)

    node_forces = np.zeros((Nnodes, 3), dtype=float)
    shape_functions = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)

    # IMPORTANT ASSUMPTION (same as your original code):
    # struct_forces[e] corresponds to the same element ordering as quads[e].
    # If your centroids are not in the same order as quads, you must remap by element IDs (EIDs).
    for e in range(Nfem.shape[0]):
        nodes = node_map[e, :]
        if np.any(nodes == -1):
            raise ValueError(f"Invalid node mapping in element {e}. Node not found in grd.")

        for j in range(4):
            node_forces[nodes[j], :] += struct_forces[e, :] * shape_functions[j]

    return node_forces, grd
