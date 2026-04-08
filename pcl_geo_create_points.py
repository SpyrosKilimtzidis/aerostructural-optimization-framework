def extractpointfromcurve(pointid, curveid, dist):
    """
    Extract a point at a specific distance along a curve.

    Parameters:
        pointid (list of int): Current list of point IDs.
        curveid (int): ID of the curve to extract from.
        dist (float): Normalized distance (0 to 1).

    Returns:
        pointid (updated list),
        key1, key2 (PCL command strings)
    """
    new_pointid = pointid[-1] + 1 if pointid else 1
    pointid.append(new_pointid)

    key1 = 'STRING asm_grid_extract_cu_created_ids[VIRTUAL]'
    key2 = (
        f'asm_const_grid_extract_v1( "{new_pointid}", "Curve {curveid}", {dist:.6f}, '
        f'1, asm_grid_extract_cu_created_ids )'
    )

    return pointid, key1, key2


def createpointintersect(pointid, curve1, curve2):
    """
    Create a point at the intersection of two curves.

    Parameters:
        pointid (list of int): List of current point IDs.
        curve1, curve2 (int): Curve IDs to intersect.

    Returns:
        pointid (updated list),
        key1, key2 (PCL command strings)
    """
    new_pointid = pointid[-1] + 1 if pointid else 1
    pointid.append(new_pointid)

    key1 = 'STRING asm_create_grid_int_created_ids[VIRTUAL]'
    key2 = (
        f'asm_const_grid_intersect_v1( "{new_pointid}", "Curve {curve1}", "Curve {curve2}", '
        f'asm_create_grid_int_created_ids )'
    )

    return pointid, key1, key2


def createpprojectcurve(pointid, curveid, p1):
    pointid.append(pointid[-1] + 1)

    key1 = 'STRING asm_create_grid_pro_created_ids[VIRTUAL]'
    key2 = f'asm_const_grid_project_v1( "{pointid[-1]}", "Point {p1}", "Curve {curveid}", 3, FALSE, 1, "<1 0 0>", "Coord 0", asm_create_grid_pro_created_ids )'

    return pointid, key1, key2


def createpointintersectcurve2plane(pointid,curve,plane):
    pointid.append(pointid[-1] + 1)
    key1 = 'STRING asm_create_grid_int_created_ids[VIRTUAL]'
    key2 = f'asm_const_grid_intersect_plo_cv( "{pointid[-1]}", "Curve {curve}", "Plane {plane}", 0., asm_create_grid_int_created_ids )'
    return pointid, key1, key2
