def createcurve2p(curveid, p1, p2):
    """
    Create a curve between two points (like asm_const_line_2point).

    Parameters:
        curveid (list of int): Current list of curve IDs.
        p1, p2 (int): Point IDs.

    Returns:
        curveid (updated list),
        key1, key2 (PCL command strings)
    """
    new_curveid = curveid[-1] + 1 if curveid else 1
    curveid.append(new_curveid)

    key1 = 'STRING asm_line_2point_created_ids[VIRTUAL]'
    key2 = (
        f'asm_const_line_2point( "{new_curveid}", '
        f'"Point {p1}", "Point {p2}", 0, "", 50., 1, asm_line_2point_created_ids )'
    )

    return curveid, key1, key2

def createcurvexyz(curveid, pointid, point, vector, unq):
    """
    Create a line from a given point in a specified vector direction.

    Parameters:
        curveid (list of int): List of curve IDs.
        pointid (list of int): List of point IDs.
        point (int): Starting point ID.
        vector (list or tuple of 3 floats): Direction vector.
        unq (int): Number of unique points to add (0, 1, or 2).

    Returns:
        curveid (updated list)
        pointid (updated list)
        key1, key2 (PCL command strings)
    """
    new_curveid = curveid[-1] + 1 if curveid else 1
    curveid.append(new_curveid)

    key1 = 'STRING asm_create_line_xyz_created_ids[VIRTUAL]'
    key2 = (
        f'asm_const_line_xyz( "{new_curveid}", "<{vector[0]:.6f} {vector[1]:.6f} {vector[2]:.6f}>", '
        f'"Point {point}", "Coord 0", asm_create_line_xyz_created_ids )'
    )

    # Conditionally update point IDs
    if pointid:
        if unq == 2:
            pointid.extend([pointid[-1] + 1, pointid[-1] + 2])
        elif unq == 1:
            pointid.append(pointid[-1] + 1)
        # else: do nothing
    else:
        if unq == 2:
            pointid.extend([1, 2])
        elif unq == 1:
            pointid.append(1)

    return curveid, pointid, key1, key2



def extractcurvefromsurf(dist, direction, pointid, curveid, surfid, unq, inst):
    """
    Mimics the behavior of the MATLAB extractcurvefromsurf function.

    Parameters:
        dist (float): Normalized spanwise position (0 to 1).
        direction (int): Extraction direction (likely normal to surface).
        pointid (list of int): Current list of point IDs.
        curveid (list of int): Current list of curve IDs.
        surfid (int): Surface ID to extract from.
        unq (int): Determines how many new point IDs to allocate.
        inst (int): Whether it's the first call or not.

    Returns:
        pointid, curveid, key1, key2: updated IDs and PCL command strings.
    """
    # Update curve ID
    if inst == 0:
        new_curveid = curveid[-1] + 1 if curveid else 1
        curveid.append(new_curveid)
    else:
        new_curveid = curveid[-1] + 1
        curveid.append(new_curveid)

    key1 = 'STRING sgm_create_curve_ex_created_ids[VIRTUAL]'
    key2 = f'sgm_const_curve_extract( "{new_curveid}", "Surface {surfid}", {direction}, {dist:.6f}, sgm_create_curve_ex_created_ids )'

    # Update point ID
    if inst == 0:
        new_points = [pointid[-1] + 1, pointid[-1] + 2] if pointid else [1, 2]
        pointid.extend(new_points)
    else:
        if unq == 2:
            new_points = [pointid[-1] + 1, pointid[-1] + 2]
            pointid.extend(new_points)
        elif unq == 1:
            new_points = [pointid[-1] + 1]
            pointid.extend(new_points)
        else:
            pass  # No new points

    return pointid, curveid, key1, key2


def breakcurves(curveid, point, curve1):
    """
    Break a curve at a specific point and generate two new curves.

    Parameters:
        curveid (list of int): Current list of curve IDs.
        point (int): Point ID where the break occurs.
        curve1 (int): ID of the curve to break.

    Returns:
        curveid (updated list),
        key1...key5 (PCL command strings)
    """
    # Add two new curve IDs
    new_curveid1 = curveid[-1] + 1 if curveid else 1
    new_curveid2 = new_curveid1 + 1
    curveid.extend([new_curveid1, new_curveid2])

    # Defensive cleanup: ensure original curve1 not duplicated
    curveid = [cid for cid in curveid if cid != curve1]

    key1 = 'STRING sgm_curve_break_poi_created_ids[VIRTUAL]'
    key2 = (
        f'sgm_edit_curve_break_point( "{new_curveid1}", "Point {point}", "Curve {curve1}", '
        f'TRUE, sgm_curve_break_poi_created_ids )'
    )
    key3 = '$# Question from application SGM'
    key4 = '$#     Do you wish to delete the original curves?'
    key5 = '$? YES 38000217 '

    return curveid, key1, key2, key3, key4, key5


def mergecurves(curveid, curve1, curve2):
    """
    Merge two curves into one, return new curve ID and related PCL commands.

    Parameters:
        curveid (list of int): Current list of curve IDs.
        curve1, curve2 (int): Curve IDs to be merged.

    Returns:
        curveid (updated list)
        key1...key5 (PCL command strings)
    """
    new_curveid = curveid[-1] + 1 if curveid else 1
    curveid.append(new_curveid)

    # Defensive cleanup (MATLAB tried to remove duplicates — you can skip or log if needed)
    curveid = [cid for cid in curveid if cid != curve1 and cid != curve2]

    key1 = 'STRING sgm_edit_curve_merg_created_ids[VIRTUAL]'
    key2 = (
        f'sgm_edit_curve_merge( "{new_curveid}", "Curve {curve1} {curve2}", 1, 0.00050000002, TRUE, '
        f'sgm_edit_curve_merg_created_ids )'
    )
    key3 = '$# Question from application SGM'
    key4 = '$#     Do you wish to delete the original curves?'
    key5 = '$? YES 38000217 '

    return curveid, key1, key2, key3, key4, key5


def createcurveintersectplanesurf(pointid, curveid, planeid, surf, ident):
    if curveid:
        curveid.append(curveid[-1] + 1)
    else:
        curveid = [1]

    if ident == 0:
        pass  # pointid unchanged
    elif ident == 1:
        pointid.append(pointid[-1] + 1)
    else:
        pointid.append(pointid[-1] + 1)
        pointid.append(pointid[-1] + 1)

    key1 = 'STRING sgm_create_curve_in_created_ids[VIRTUAL]'
    key2 = f'sgm_const_curve_intersect( "{curveid[-1]}", 2, "Plane {planeid}", "Surface {surf}", 0.0049999999, 0.00050000001, sgm_create_curve_in_created_ids )'

    return pointid, curveid, key1, key2


def reversecurve(curve):
    key1 = 'STRING sgm_edit_curve_rev_reversed_ids[VIRTUAL]'
    key2 = f'sgm_edit_curve_reverse( TRUE, "Curve {curve}", sgm_edit_curve_rev_reversed_ids)'
    return key1, key2



