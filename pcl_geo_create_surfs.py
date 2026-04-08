def createsurf3(surfid, curvelist):
    """
    Create a 3-edge surface from a list of 3 curve IDs.

    Parameters:
        surfid (list of int): List of surface IDs.
        curvelist (list of 3 ints): Curve IDs that define the surface edges.

    Returns:
        surfid (updated list)
        key1 through key6 (PCL command strings)
    """
    new_surfid = surfid[-1] + 1 if surfid else 1
    surfid.append(new_surfid)

    key1 = 'STRING sgm_surface_3edge_created_ids[VIRTUAL]'
    key2 = (
        f'sgm_const_surface_3edge( "{new_surfid}", "Curve {curvelist[0]}", '
        f'"Curve {curvelist[1]}", "Curve {curvelist[2]}", sgm_surface_3edge_created_ids )'
    )
    key3 = '$# Question from application SGM'
    key4 = '$#     Warning.  This surface may not be meshable.  Do you wish to create it '
    key5 = '$# anyway?'
    key6 = '$? YES 38000661 '

    return surfid, key1, key2, key3, key4, key5, key6

def createsurf4(surfid, curvelist):
    """
    Create a 4-edge surface from a list of 4 curve IDs.

    Parameters:
        surfid (list of int): Current list of surface IDs.
        curvelist (list of 4 ints): Curve IDs that define the surface edges.

    Returns:
        surfid (updated list)
        key1, key2 (PCL command strings)
    """
    new_surfid = surfid[-1] + 1 if surfid else 1
    surfid.append(new_surfid)

    key1 = 'STRING sgm_surface_4edge_created_ids[VIRTUAL]'
    key2 = (
        f'sgm_const_surface_4edge( "{new_surfid}", '
        f'"Curve {curvelist[0]}", "Curve {curvelist[1]}", '
        f'"Curve {curvelist[2]}", "Curve {curvelist[3]}", '
        f'sgm_surface_4edge_created_ids )'
    )

    return surfid, key1, key2
