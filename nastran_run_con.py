import numpy as np
from thick2lammsymm import select_laminates_single
from geo_settol import settol
import pcl_geo_create_curves as pcl_curves
import pcl_geo_create_points as pcl_points
import pcl_geo_create_surfs as pcl_surfs
import pcl_geo_create_planes as pcl_planes
import pcl_meshing as pcl_mesh
import pcl_create_materials_failure as pcl_create_matfail
import pcl_properties as pcl_prop
import pcl_geo_create_coords as pcl_coords
import pcl_bc
import pcl_groups
from writebdfoutput import write2bdfoutput
import subprocess
import time
import psutil
import math
from collect_output import collect_output
import conservative_mapping_interpolation as cfd_interpolation
import read_fluent_data as read_fluent
import generate_loadcase_commands as gen_lc_comms
from delete_files import delete_files

def run_nastran_simulation(fs_pos, rs_pos,rib_space,stringer_space, t_skins, t_fs, t_rs, t_ribs, half_span, a_stall):

    mat = {
        'WBX': {
            'E1': 152e9,
            'E2': 6.9e9,
            'nu12': 0.4,
            'G12': 3.3e9,
            'G23': 1.9e9,
            'G13': 3.3e9,
            'rho': 1520,
            'Xt': 2255e6,
            'Xc': 667.9e6,
            'Yt': 17.3e6,
            'Yc': 83.1e6,
            'S': 42.7e6,
            'SB': 39.4e6,
            'iterm': 0,
            'nmat': 'mat1',
            'hply': 0.202e-3
        },
        'LE': {
            'E1': 64.3e9,
            'E2': 64.3e9,
            'nu12': 0.08,
            'G12': 4.7e9,
            'G23': 5.2e9,
            'G13': 5.2e9,
            'rho': 1510,
            'Xt': 727.3e6,
            'Xc': 342.2e6,
            'Yt': 727.3e6,
            'Yc': 342.2e6,
            'S': 47.4e6,
            'SB': 22.9e6,
            'iterm': 0,
            'nmat': 'mat2',
            'hply': 0.363e-3
        }
    }
    # Extract design variables from x
    d1 = fs_pos
    d2 = rs_pos
    strspacing0 = stringer_space
    w1 = 0.5e-3
    w2 = 0.5e-3
    w3 = 0.5e-3
    nr_ribs = round(half_span/rib_space)
    vtot = np.linspace(0, 1,nr_ribs )
    # Stringer spacing vector
    strspacingtot = np.arange(0, 1 + strspacing0, strspacing0)
    # Safeguards you can tune
    gmin = 0.05  # min distance to each spar in parametric space (0..1)
    hmin = 0.02  # optional: minimum allowed gap between partitions (prevents slivers)
    force_one = True  # keep at least one stringer (midpoint) if otherwise none

    # usable interval after clearance
    L_eff = 1.0 - 2.0 * gmin
    if L_eff <= 0.0:
        # no room at all
        strspacing = np.array([])
        strspacingtot = np.array([0.0, 1.0])
    else:
        # Choose N (number of interior stringers) so that uniform gap <= strspacing0
        # gap = L_eff / (N+1) <= strspacing0  -> N >= L_eff/strspacing0 - 1
        N_min = max(0, int(math.ceil(L_eff / strspacing0) - 1))

        # Also enforce gap >= hmin (optional but recommended)
        # gap = L_eff / (N+1) >= hmin -> N <= L_eff/hmin - 1
        N_max = max(0, int(math.floor(L_eff / hmin) - 1))

        if N_min > N_max:
            # No feasible N that satisfies both max-spacing and min-gap constraints.
            # In practice: keep fewer stringers (respect hmin) and accept slightly larger gap.
            N = N_max
        else:
            # Pick minimal N that satisfies gap <= strspacing0 (stable topology, fewer stringers)
            N = N_min

        if N <= 0:
            if force_one and (L_eff > 0.0):
                # one stringer at the center of the usable region
                strspacing = np.array([0.5])
            else:
                strspacing = np.array([])
        else:
            # symmetric interior points in [gmin, 1-gmin]
            strspacing = np.linspace(gmin, 1.0 - gmin, N + 2)[1:-1]
        # Full partition vector including spars
        strspacingtot = np.concatenate(([0.0], strspacing, [1.0]))
    # --- End robust stringer spacing generation ---
    num1 = 0.01
    num2 = 0.01
    tol = 5e-3
    dextracut = 0.9

    xn = np.array([t_skins, t_fs, t_rs, t_ribs, 0.606e-3, 0.606e-3, 0.606e-3])                  # Target thicknesses
    hply = (0.202e-3, 0.363e-3)                                   # WBX, LE ply thicknesses
    compid = np.array([1, 1, 1, 1, 1, 1, 1])                      # WBX or LE for each panel
    laminate, thickness = select_laminates_single(xn, hply, compid)
    # Initialize
    data1 = [[]]
    surfid = list(range(1, 7))
    direction = 2
    inst = 0
    nr = 0
    n_ribs = len(vtot)
    lepoint = np.zeros((n_ribs, 1), dtype=int)
    tepoint = np.zeros((n_ribs, 2), dtype=int)
    cte = np.zeros((n_ribs, 1), dtype=int)
    wgbxp = np.zeros((n_ribs, 2 * (len(strspacingtot) + 1)), dtype=int)
    cint = np.zeros((n_ribs, len(strspacingtot) + 1), dtype=int)
    n_spanwise_curves = 1 + 2 + 2 * (len(strspacingtot) + 1)
    n_spanwise_curves_upper_lower = len(strspacingtot) + 2
    csect = np.zeros((n_ribs - 1, n_spanwise_curves), dtype=int)
    csectu = np.zeros((n_ribs - 1, n_spanwise_curves_upper_lower+1), dtype=int)
    csectl = np.zeros((n_ribs - 1, n_spanwise_curves_upper_lower+1), dtype=int)
    cs1 = np.zeros((n_ribs, n_spanwise_curves_upper_lower), dtype=int)
    cs2 = np.zeros((n_ribs, n_spanwise_curves_upper_lower ), dtype=int)
    ribs = np.zeros((n_ribs, n_spanwise_curves_upper_lower), dtype=int)
    us = np.zeros((n_ribs-1, n_spanwise_curves_upper_lower), dtype=int)
    ls = np.zeros((n_ribs-1, n_spanwise_curves_upper_lower), dtype=int)
    fs = np.zeros((n_ribs-1, n_spanwise_curves_upper_lower), dtype=int)
    rs = np.zeros((n_ribs-1, n_spanwise_curves_upper_lower), dtype=int)
    ste = np.zeros((n_ribs-1, n_spanwise_curves_upper_lower), dtype=int)
    cursurf = [1, 2, 3]
    axis = 2
    curveid = [1]
    lecurve = 1
    pointid = [1, 2]
    identu = 1
    identl = 1
    coordid = []
    planeid = []
    # Main Loop
    while nr < len(vtot):
        data1[0].extend(settol(5e-4))
        if nr == 0:
            curpoint = pointid[0]
        elif nr == len(vtot)-1:
            curpoint = pointid[1]
        else:
            pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, lecurve, vtot[nr])
            data1[0].extend([cmd1, cmd2])
            curpoint = pointid[-1]
        if nr == 0:
            unq = 1
            inst = 1
            dist = 0
            direction = 2
            # Extract upper curve
            pointid, curveid, cmd1, cmd2 = pcl_curves.extractcurvefromsurf(dist, direction, pointid, curveid, cursurf[0], unq, inst)
            data1[0].extend([cmd1, cmd2])
            # Save upper surface curve ID
            c1 = curveid[-1]
            # Extract lower curve
            pointid, curveid, cmd1, cmd2 = pcl_curves.extractcurvefromsurf(dist, direction, pointid, curveid,cursurf[1], unq, inst)
            data1[0].extend([cmd1, cmd2])
            c2 = curveid[-1]
            lepoint[nr, 0] = curpoint
            tepoint[nr, 0] = pointid[-2]
            tepoint[nr, 1] = pointid[-1]
            # Create TE curve
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, tepoint[nr, 0],  tepoint[nr, 1])
            data1[0].extend([cmd1, cmd2])
            cte[nr, 0] = curveid[-1]
        else:
            coordid, cmd1, cmd2 = pcl_coords.createcoordeuler(coordid, curpoint)
            data1[0].extend([cmd1, cmd2])
            planeid, cmd1, cmd2 = pcl_planes.createplanepoints(planeid, coordid[-1], curpoint, axis)
            data1[0].extend([cmd1, cmd2])
            k1, k2, k3 = settol(4e-5)
            data1[0].extend([k1, k2, k3])
            pointid, curveid, cmd1, cmd2 = pcl_curves.createcurveintersectplanesurf(pointid, curveid, planeid[-1], cursurf[0], identu)
            data1[0].extend([cmd1, cmd2])
            c1 = curveid[-1]
            pointid, curveid, cmd1, cmd2 = pcl_curves.createcurveintersectplanesurf(pointid, curveid, planeid[-1], cursurf[1], identl)
            data1[0].extend([cmd1, cmd2])
            c2 = curveid[-1]
            lepoint[nr, 0] = curpoint
            tepoint[nr, 0] = pointid[-2]
            tepoint[nr, 1] = pointid[-1]
            # Create TE curve
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, tepoint[nr, 0], tepoint[nr, 1])
            data1[0].extend([cmd1, cmd2])
            cte[nr, 0] = curveid[-1]
        # Create midpoint in TE upper and lower
        k1, k2, k3 = settol(3e-5)
        data1[0].extend([k1, k2, k3])
        pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, curveid[-1], 0.5)
        data1[0].extend([cmd1, cmd2])
        k1, k2, k3 = settol(5e-4)
        data1[0].extend([k1, k2, k3])
        # Create local chord line (from LE to midpoint)
        curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, lepoint[nr, 0], pointid[-1])
        data1[0].extend([cmd1, cmd2])
        # Extract spar point 1 (at d1 along chord line)
        pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, curveid[-1], d1)
        data1[0].extend([cmd1, cmd2])
        # Extract spar point 2 (at d2 along chord line)
        pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, curveid[-1], d2)
        data1[0].extend([cmd1, cmd2])
        # Create spar curve (between the last two points)
        curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, pointid[-2], pointid[-1])
        data1[0].extend([cmd1, cmd2])
        curcur = curveid[-1]  # Store current spar curve ID
        curp = [pointid[-2], pointid[-1]]  # Start and end points of spar
        curptot1 = []  # Intermediate stringer points
        # Extract stringer points along the spar
        idx1 = 0
        for spacing in strspacing:
            pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, curveid[-1], spacing)
            data1[0].extend([cmd1, cmd2])
            curptot1.append(pointid[-1])
        # TE Break point at location `dextracut` from the second-to-last curve (spar curve)
        pointid, cmd1, cmd2 = pcl_points.extractpointfromcurve(pointid, curveid[-2], dextracut)
        data1[0].extend([cmd1, cmd2])
        # Combine all points: start/end of spar, stringer points, TE break point
        curptot = curp + curptot1 + [pointid[-1]]
        k1, k2, k3 = settol(8e-6)
        data1[0].extend([k1, k2, k3])
        for pt in curptot:
            planeid, cmd1, cmd2 = pcl_planes.createplanepoints(planeid, 0, pt, 1)
            data1[0].extend([cmd1, cmd2])
            pointid, cmd1, cmd2 = pcl_points.createpointintersectcurve2plane(pointid, c1, planeid[-1])
            data1[0].extend([cmd1, cmd2])
            wgbxp[nr, idx1] = pointid[-1]
            idx1 += 1
            pointid, cmd1, cmd2 = pcl_points.createpointintersectcurve2plane(pointid, c2, planeid[-1])
            data1[0].extend([cmd1, cmd2])
            wgbxp[nr, idx1] = pointid[-1]
            idx1 += 1
        wgbxp[nr, :] = np.concatenate([
            wgbxp[nr, 0:2],  # first two points
            wgbxp[nr, 4:-2],  # mid points excluding swapped ones
            wgbxp[nr, 2:4],  # swapped two mid points
            wgbxp[nr, -2:]  # last two points
        ])
        ccint = 0
        # Create chordwise curves
        for kk in range(0, wgbxp.shape[1], 2):  # kk = 0, 2, 4, ...
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, wgbxp[nr, kk], wgbxp[nr, kk + 1])
            data1[0].append(cmd1)
            data1[0].append(cmd2)
            cint[nr][ccint] = curveid[-1]
            ccint += 1
        if nr > 0:
            countccc = 0
            # Leading Edge spanwise connection
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, lepoint[nr - 1, 0], lepoint[nr, 0])
            data1[0].extend([cmd1, cmd2])
            csect[nr-1, countccc] = curveid[-1]
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, tepoint[nr - 1, 0], tepoint[nr , 0])
            data1[0].extend([cmd1, cmd2])
            csect[nr - 1, countccc + 1] = curveid[-1]
            curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, tepoint[nr - 1, 1], tepoint[nr, 1])
            data1[0].extend([cmd1, cmd2])
            csect[nr - 1, countccc + 2] = curveid[-1]
            countccc += 3
            # Odd indices (kk = 0, 2, 4, ...)
            for kk in range(0, wgbxp.shape[1], 2):
                curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, wgbxp[nr - 1, kk], wgbxp[nr, kk])
                data1[0].extend([cmd1, cmd2])
                csect[nr - 1, countccc] = curveid[-1]
                countccc += 1
            # Even indices (kk = 1, 3, 5, ...)
            for kk in range(1, wgbxp.shape[1], 2):
                curveid, cmd1, cmd2 = pcl_curves.createcurve2p(curveid, wgbxp[nr - 1, kk], wgbxp[nr, kk])
                data1[0].extend([cmd1, cmd2])
                csect[nr - 1, countccc] = curveid[-1]
                countccc += 1
            csect[nr - 1, :] = np.concatenate((
                [csect[nr - 1, 0]],
                csect[nr - 1, 3:],  # 4th to end
                csect[nr - 1, 1:3]  # 2nd and 3rd
            ))
            csectu[nr - 1, :] = np.concatenate((
                csect[nr - 1, :cs1.shape[1]],  # First `n` curves (cs1.shape[1] == # of upper ribs)
                [csect[nr - 1, -2]]  # Second-to-last for wraparound
            ))
            csectl[nr - 1, :] = np.concatenate((
                [csect[nr - 1, 0]],  # First LE curve
                csect[nr - 1, cs1.shape[1]:-2],  # Lower rib region
                [csect[nr - 1, -1]]  # Last TE curve
            ))
        # Break upper airfoil section curves
        [cmd1, cmd2] = pcl_curves.reversecurve(c2)
        data1[0].extend([cmd1, cmd2])
        countc = 0
        k1, k2, k3 = settol(5e-4)
        data1[0].extend([k1, k2, k3])
        for kk in range(0, wgbxp.shape[1], 2):
            if kk == 0:
                curveid, k1, k2, k3, k4, k5 = pcl_curves.breakcurves(curveid, wgbxp[nr, kk], c1)
            else:
                curveid, k1, k2, k3, k4, k5 = pcl_curves.breakcurves(curveid, wgbxp[nr, kk], curveid[-2])

            data1[0].extend([k1, k2, k3, k4, k5])
            cs1[nr, countc] = curveid[-1]
            countc += 1
            # Add last upper curve
            cs1[nr, countc] = curveid[-2]
        # Break lower airfoil section curves
        countc = 0
        for kk in range(1, wgbxp.shape[1], 2):  # kk = 1, 3, 5, ...
            if kk == 1:
                curveid, k1, k2, k3, k4, k5 = pcl_curves.breakcurves(curveid, wgbxp[nr, kk], c2)
            else:
                curveid, k1, k2, k3, k4, k5 = pcl_curves.breakcurves(curveid, wgbxp[nr, kk], curveid[-2])
            data1[0].extend([k1, k2, k3, k4, k5])
            cs2[nr, countc] = curveid[-1]
            countc += 1
        # Add last lower curve
        cs2[nr, countc] = curveid[-2]
        cint1 = np.hstack([cint, cte.reshape(-1, 1)])  # Ensure cte is column vector
        k1, k2, k3 = settol(5e-5)
        data1[0].extend([k1, k2, k3])
        for kk in range(cs1.shape[1]):
            if kk == 0:
                surfid, k1, k2, k3, k4, k5, k6 = pcl_surfs.createsurf3(
                    surfid, [cs1[nr, kk], cs2[nr, kk], cint1[nr, 0]]
                )
                data1[0].extend([k1, k2, k3, k4, k5, k6])
            else:
                surfid, k1, k2 = pcl_surfs.createsurf4(
                    surfid, [cs1[nr, kk], cint1[nr, kk], cs2[nr, kk], cint1[nr, kk - 1]]
                )
                data1[0].extend([k1, k2])

            ribs[nr, kk] = surfid[-1]
        k1, k2, k3 = settol(5e-4)
        data1[0].extend([k1, k2, k3])
        if nr > 0:
            for kk in range(cs1.shape[1]):
                # Upper Spanwise Surface (US)
                surfid, k1, k2 = pcl_surfs.createsurf4(surfid, [csectu[nr - 1, kk], csectu[nr - 1, kk + 1], cs1[nr, kk], cs1[nr - 1, kk]])
                data1[0].extend([k1, k2])
                us[nr - 1, kk] = surfid[-1]

                # Lower Spanwise Surface (LS)
                surfid, k1, k2 = pcl_surfs.createsurf4(
                    surfid, [csectl[nr - 1, kk], csectl[nr - 1, kk + 1], cs2[nr, kk], cs2[nr - 1, kk]]
                )
                data1[0].extend([k1, k2])
                ls[nr - 1, kk] = surfid[-1]
            # Forward Spar (FS)
            surfid, k1, k2 = pcl_surfs.createsurf4(
                surfid, [cint[nr - 1, 0], cint[nr, 0], csectu[nr - 1, 1], csectl[nr - 1, 1]]
            )
            data1[0].extend([k1, k2])
            fs[nr - 1, kk] = surfid[-1]

            # Rear Spar (RS)
            surfid, k1, k2 = pcl_surfs.createsurf4(
                surfid, [cint[nr - 1, -2], cint[nr, -2], csectu[nr - 1, -3], csectl[nr - 1, -3]]
            )
            data1[0].extend([k1, k2])
            rs[nr - 1, kk] = surfid[-1]

            # Trailing Edge (TE)
            k1, k2, k3 = settol(1e-5)
            data1[0].extend([k1, k2, k3])
            surfid, k1, k2 = pcl_surfs.createsurf4(
                surfid, [cte[nr - 1, 0], cte[nr, 0], csectu[nr - 1, -1], csectl[nr - 1, -1]]
            )
            data1[0].extend([k1, k2])
            k1, k2, k3 = settol(5e-4)
            data1[0].extend([k1, k2, k3])
            ste[nr - 1, kk] = surfid[-1]
        nr = nr+1
    # Meshing
    numelte = 1
    # Create mesh seeds for trailing edge curves
    for kk in range(len(cte)):
        cmd = pcl_mesh.createmseed(cte[kk], numelte)
        data1[0].append(cmd)
    # Flatten all ribs into a 1D array
    ribs_tot = ribs.flatten()
    # Get last column (trailing edge ribs)
    ribs_te = ribs[:, -1]
    # Subtract trailing edge ribs from all ribs
    ribs_red = np.setdiff1d(ribs_tot, ribs_te)
    # Mesh TE surfaces (ribs_red)
    cmd1, cmd2, cmd3, cmd4, cmd5 =  pcl_mesh.createmeshsurf(ribs_red, num2)
    data1[0].extend([cmd1, cmd2, cmd3, cmd4, cmd5])
    # Equivalence nodes
    cmd1, cmd2, cmd3 = pcl_mesh.equivalencenodes(tol)
    data1[0].extend([cmd1, cmd2, cmd3])
    # Mesh TE ribs (ribs_te)
    cmd1, cmd2, cmd3, cmd4, cmd5 = pcl_mesh.createmeshsurf(ribs_te, num2)
    data1[0].extend([cmd1, cmd2, cmd3, cmd4, cmd5])
    # Equivalence again
    cmd1, cmd2, cmd3 = pcl_mesh.equivalencenodes(tol)
    data1[0].extend([cmd1, cmd2, cmd3])
    # Mesh remaining surfaces
    Cexc = np.setdiff1d(surfid[6:], ribs_tot)  # Python uses 0-based indexing
    key1, key2, key3, key4, key5 = pcl_mesh.createmeshsurf(Cexc, num1)
    data1[0].extend([key1, key2, key3, key4, key5])
    # Equivalence nodes
    key1, key2, key3 = pcl_mesh.equivalencenodes(tol)
    data1[0].extend([key1, key2, key3])
    # Mesh beam curves
    beamcurves_u = csectu[:, 1:-2].reshape(-1)
    beamcurves_l = csectl[:, 1:-2].reshape(-1)
    beamcurves = np.concatenate((beamcurves_u, beamcurves_l))
    key1, key2, key3, key4, key5 = pcl_mesh.createmeshcurve(beamcurves, num1)
    data1[0].extend([key1, key2, key3, key4, key5])
    # Final node equivalence
    key1, key2, key3 = pcl_mesh.equivalencenodes(tol)
    data1[0].extend([key1, key2, key3])
    # MATERIALS
    # %% Materials & Properties
    # WBX
    data1[0].append(pcl_create_matfail.createmat8(
        mat["WBX"]["nmat"], mat["WBX"]["E1"], mat["WBX"]["E2"], mat["WBX"]["nu12"],
        mat["WBX"]["G12"], mat["WBX"]["G23"], mat["WBX"]["G13"], mat["WBX"]["rho"]))
    data1[0].append(pcl_create_matfail.createmat8_fail(
        mat["WBX"]["nmat"], mat["WBX"]["Xt"], mat["WBX"]["Xc"], mat["WBX"]["Yt"],
        mat["WBX"]["Yc"], mat["WBX"]["S"], mat["WBX"]["SB"], mat["WBX"]["iterm"]))
    # Fabric
    data1[0].append(pcl_create_matfail.createmat8(
        mat["LE"]["nmat"], mat["LE"]["E1"], mat["LE"]["E2"], mat["LE"]["nu12"],
        mat["LE"]["G12"], mat["LE"]["G23"], mat["LE"]["G13"], mat["LE"]["rho"]))
    data1[0].append(pcl_create_matfail.createmat8_fail(
        mat["LE"]["nmat"], mat["LE"]["Xt"], mat["LE"]["Xc"], mat["LE"]["Yt"],
        mat["LE"]["Yc"], mat["LE"]["S"], mat["LE"]["SB"], mat["LE"]["iterm"]))
    cmd1, cmd2 = pcl_coords.createcoord(coordid, lepoint[0, 0], tepoint[0, 0], lepoint[1, 0])
    data1[0].extend([cmd1, cmd2])
    # PROPERTIES
    data1[0].append(pcl_prop.createpcompprops(laminate[0], mat['WBX']['nmat'], 'laminateus', hply[0], 2))
    data1[0].append(pcl_prop.createpcompprops(laminate[0], mat['LE']['nmat'], 'laminateus_LE_TE', hply[1], 2))
    # Upper Skin
    skinsurfs_1_us = np.concatenate([us.flatten(), ste[:, -1]])
    le_te_surfs_upper = np.concatenate([us[:, 0], us[:, -1]])
    skinsurfs_2_us = np.setdiff1d(skinsurfs_1_us, le_te_surfs_upper)
    data1[0].append(pcl_prop.createpshellprops('us', 'laminateus', 1, skinsurfs_2_us))
    data1[0].append(pcl_prop.createpshellprops('us_LE_TE', 'laminateus_LE_TE', 1, le_te_surfs_upper))
    # Lower Skin
    data1[0].append(pcl_prop.createpcompprops(laminate[0], mat['WBX']['nmat'], 'laminatels', hply[0], 2))
    data1[0].append(pcl_prop.createpcompprops(laminate[0], mat['LE']['nmat'], 'laminatels_LE_TE', hply[1], 2))
    skinsurfs_1_ls = ls.flatten()
    le_te_surfs_lower = np.concatenate([ls[:, 0], ls[:, -1]])
    skinsurfs_2_ls = np.setdiff1d(skinsurfs_1_ls, le_te_surfs_lower)
    data1[0].append(pcl_prop.createpshellprops('ls', 'laminatels', 1, skinsurfs_2_ls))
    data1[0].append(pcl_prop.createpshellprops('ls_LE_TE', 'laminatels_LE_TE', 1, le_te_surfs_lower))
    # Front Spar
    data1[0].append(pcl_prop.createpcompprops(laminate[1], mat['WBX']['nmat'], 'laminatefs', hply[0], 2))
    data1[0].append(pcl_prop.createpshellprops('fsw', 'laminatefs', 1, fs[:, -1]))
    # Rear Spar
    data1[0].append(pcl_prop.createpcompprops(laminate[2], mat['WBX']['nmat'], 'laminaters', hply[0], 2))
    data1[0].append(pcl_prop.createpshellprops('rsw', 'laminaters', 1, rs[:, -1]))
    # Ribs
    data1[0].append(pcl_prop.createpcompprops(laminate[3], mat['WBX']['nmat'], 'laminaterw', hply[0], 2))
    data1[0].append(pcl_prop.createpshellprops('riw', 'laminaterw', 0, ribs.flatten()))
    # FSC
    orientsc = [0, 0, 1]

    Efsc, Gfsc = pcl_create_matfail.eq_moduliv1(laminate[4], mat['WBX']['nmat'], hply,
                             mat['WBX']['E1'], mat['WBX']['E2'],
                             mat['WBX']['G12'], mat['WBX']['nu12'])
    [cmd1, cmd2] = pcl_create_matfail.createmat1('matfsc', Efsc, Gfsc, mat['WBX']['rho'],
                               mat['WBX']['Xt'], mat['WBX']['Xc'], mat['WBX']['S'])
    data1[0].extend([cmd1, cmd2])
    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectu[:, 1], w1, thickness[4], 'matfsc',
                                     'fscu', 'sect1u', orientsc, thickness[0], 0)
    data1[0].extend([cmd1, cmd2])

    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectl[:, 1], w1, thickness[4], 'matfsc',
                                     'fscl', 'sect1l', orientsc, thickness[0], 1)
    data1[0].extend([cmd1, cmd2])
    # RSC
    Ersc, Grsc = pcl_create_matfail.eq_moduliv1(laminate[5], mat['WBX']['nmat'], hply,
                             mat['WBX']['E1'], mat['WBX']['E2'],
                             mat['WBX']['G12'], mat['WBX']['nu12'])
    [cmd1, cmd2] = pcl_create_matfail.createmat1('matrsc', Ersc, Grsc, mat['WBX']['rho'],
                               mat['WBX']['Xt'], mat['WBX']['Xc'], mat['WBX']['S'])
    data1[0].extend([cmd1, cmd2])
    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectu[:, -3], w2, thickness[5], 'matrsc',
                                     'rscu', 'sect2u', orientsc, thickness[0], 0)
    data1[0].extend([cmd1, cmd2])
    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectl[:, -3], w2, thickness[5], 'matrsc',
                                     'rscl', 'sect2l', orientsc, thickness[0], 1)
    data1[0].extend([cmd1, cmd2])
    # STR
    orientst = [1, 0, 0]

    Est, Gst =  pcl_create_matfail.eq_moduliv1(laminate[6], mat['WBX']['nmat'], hply,
                           mat['WBX']['E1'], mat['WBX']['E2'],
                           mat['WBX']['G12'], mat['WBX']['nu12'])
    [cmd1, cmd2] = pcl_create_matfail.createmat1('matst', Est, Gst, mat['WBX']['rho'],
                               mat['WBX']['Xt'], mat['WBX']['Xc'], mat['WBX']['S'])
    data1[0].extend([cmd1, cmd2])
    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectu[:, 2:-3].ravel(), w3, thickness[6], 'matst','stu', 'sect3u', orientst, thickness[0], 0)
    data1[0].extend([cmd1, cmd2])
    [cmd1, cmd2] = pcl_prop.createcbeamprops(csectl[:, 2:-3].ravel(), w3, thickness[6], 'matst','stl', 'sect3l', orientst, thickness[0], 1)
    data1[0].extend([cmd1, cmd2])
    data1[0].append(pcl_bc.createclamp(ribs[0, :]))
    # Create groups needed
    skins_total = np.concatenate([us.flatten(),ls.flatten()])
    cmd1, cmd2, cmd3 = pcl_groups.creategroup('wing_skins', skins_total)
    data1[0].extend([cmd1, cmd2, cmd3])
    ste_1 = ste[:, -1].flatten()
    cmd1, cmd2, cmd3 = pcl_groups.creategroup('wing_te', ste_1)
    data1[0].extend([cmd1, cmd2, cmd3])
    tiprib=ribs[-1,:].flatten()
    cmd1, cmd2, cmd3 = pcl_groups.creategroup('wing_tiprib', tiprib)
    data1[0].extend([cmd1, cmd2, cmd3])
    # Write groups to bdf
    data1 = write2bdfoutput('wing_skins', 'wing_skins', data1)
    data1 = write2bdfoutput('wing_te', 'wing_te', data1)
    data1 = write2bdfoutput('wing_tiprib', 'wing_tiprib', data1)
    # Read the base SES file
    with open("patran_ini.ses", "r") as f_ini:
        data_ini = f_ini.readlines()
    # Merge and write to new SES file
    with open("patran_upd_ini.ses", "w") as f_out:
        # Write existing lines, ensuring each ends with newline
        for line in data_ini:
            if not line.endswith("\n"):
                line += "\n"
            f_out.write(line)

        # Add a newline to separate blocks
        f_out.write("\n")

        # Write new PCL commands from data1[0], ensuring each ends with newline
        for line in data1[0]:
            if not line.endswith("\n"):
                line += "\n"
            f_out.write(line)

    # Run current ses file in Patran
    # Define the command as a list of arguments
    patran_command = [r"C:\MSC.Software\Patran_x64\20190\bin\patran.exe", "-sfp", "patran_upd_ini.ses", "-ans", "yes","-b"]
    # Execute the command
    subprocess.run(patran_command, check=True)
    # Read CFD data for skins, TE and tip rib
    filename1 = f"pressure_skin_aoa{a_stall}"
    filename2 = f"pressure_te_aoa{a_stall}"
    filename3 = f"pressure_tip_aoa{a_stall}"
    aero_points_skins, aero_forces_skins, aero_normals_skins  = read_fluent.read_aero_data(filename1)
    aero_points_te, aero_forces_te, aero_normals_te  = read_fluent.read_aero_data(filename2)
    aero_points_tip, aero_forces_tip, aero_normals_tip  = read_fluent.read_aero_data(filename3)
    # Interpolate CFD to structural
    node_forces_skins, grd_skins = cfd_interpolation.conservative_mapping_projection(aero_points_skins, aero_forces_skins, aero_normals_skins, "wing_skins.bdf")
    node_forces_te, grd_te = cfd_interpolation.conservative_mapping_projection(aero_points_te, aero_forces_te, aero_normals_te, "wing_te.bdf")
    node_forces_tip, grd_tip = cfd_interpolation.conservative_mapping_projection(aero_points_tip, aero_forces_tip, aero_normals_tip, "wing_tiprib.bdf")

    # Write loads
    ses_wing = gen_lc_comms.generate_loadcase_commands("wing", grd_skins, node_forces_skins)
    ses_te = gen_lc_comms.generate_loadcase_commands("te", grd_te, node_forces_te)
    ses_tip = gen_lc_comms.generate_loadcase_commands("tip", grd_tip, node_forces_tip)

    # Merge all
    combined_ses = ses_wing[0] + ses_te[0] + ses_tip[0]

    with open("wing_upd.ses", "w") as f:
        for line in combined_ses:
            f.write(line + "\n")
    additional_files = ["sol101.ses", "sol105.ses"]
    # Open main file in append mode
    with open("wing_upd.ses", "a") as outfile:
        for fname in additional_files:
            with open(fname, "r") as infile:
                contents = infile.read()
                outfile.write("\n" + contents)  # Add newline for separation

    # Define the command as a list of arguments
    patran_cmd = [
        r'C:\MSC.Software\Patran_x64\20190\bin\patran.exe',
        '-ifile', 'init_fld.pcl',
        '-skin',
        '-db', 'wing_ini.db',
        '-sfp', 'wing_upd.ses',
        '-ans', 'yes',
        '-b'
    ]

    # Run the command without showing the console output
    result = subprocess.run(patran_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Step 1: Run MSC Nastran in the background
    nastran_cmd = [
        r'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe',
        'wing_linear_static.bdf'
    ]

    subprocess.Popen(nastran_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Step 2: Monitor if 'nastran.exe' is running
    flag = True
    while flag:
        flag = any(proc.name().lower() == 'nastran.exe' for proc in psutil.process_iter(['name']))
        time.sleep(1)  # Wait 1 second before checking again
    # Step 2: Run MSC Nastran in the background
    time.sleep(10)
    nastran_cmd = [
        r'C:\MSC.Software\MSC_Nastran\20190\bin\nastranw.exe',
        'wing_linear_buckling.bdf'
    ]
    subprocess.Popen(nastran_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    flag = True
    while flag:
        flag = any(proc.name().lower() == 'nastran.exe' for proc in psutil.process_iter(['name']))
        time.sleep(1)
    time.sleep(10)
    f_val, g_shells, g_beams, g_buck = collect_output()
    delete_files()
    return f_val, g_shells, g_beams, g_buck
