import numpy as np
def createmseed(curveid, numelms):
    # Ensure curveid is a scalar integer
    if isinstance(curveid, (list, np.ndarray)):
        curveid = int(np.ravel(curveid)[0])
    return f'mesh_seed_create( "Curve {curveid}", 1, {numelms}, 0., 0., 0. )'

def createmeshsurf(surfs, num):
    # Convert list of surface IDs to string
    curves_str = ' '.join(str(s) for s in surfs)

    key1 = 'INTEGER fem_create_mesh_surfa_num_nodes'
    key2 = 'INTEGER fem_create_mesh_surfa_num_elems'
    key3 = 'STRING fem_create_mesh_s_nodes_created[VIRTUAL]'
    key4 = 'STRING fem_create_mesh_s_elems_created[VIRTUAL]'
    key5 = f'fem_create_mesh_surf_4( "IsoMesh", 49152, "Surf {curves_str}", 1, ["{num:.6f}"], "Quad4", "#", "#", "Coord 0", "Coord 0", fem_create_mesh_surfa_num_nodes, fem_create_mesh_surfa_num_elems, fem_create_mesh_s_nodes_created, fem_create_mesh_s_elems_created )'

    return key1, key2, key3, key4, key5

def equivalencenodes(tol):
    key1 = 'REAL fem_equiv_all_x_equivtol_ab'
    key2 = 'INTEGER fem_equiv_all_x_segment'
    key3 = f'fem_equiv_all_group4( [" "], 0, "", 1, 1, {tol:.6f}, FALSE, fem_equiv_all_x_equivtol_ab, fem_equiv_all_x_segment )'
    return key1, key2, key3

def createmeshcurve(curves, num):
    key1 = 'INTEGER fem_create_mesh_curve_num_nodes'
    key2 = 'INTEGER fem_create_mesh_curve_num_elems'
    key3 = 'STRING fem_create_mesh_c_nodes_created[VIRTUAL]'
    key4 = 'STRING fem_create_mesh_c_elems_created[VIRTUAL]'

    # Convert list of curves to space-separated string
    if isinstance(curves, (list, np.ndarray)):
        curves_str = ' '.join(str(c) for c in np.ravel(curves))
    else:
        curves_str = str(curves)

    key5 = (
        f'fem_create_mesh_curv_1( "Curve {curves_str}", 16384, {num}, '
        f'"Bar2", "#", "#", "Coord 0", "Coord 0", '
        f'fem_create_mesh_curve_num_nodes, fem_create_mesh_curve_num_elems, '
        f'fem_create_mesh_c_nodes_created, fem_create_mesh_c_elems_created )'
    )

    return key1, key2, key3, key4, key5
