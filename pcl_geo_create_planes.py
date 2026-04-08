def createplanepoints(planeid, coordid, p1, axis):
    if planeid:
        planeid.append(planeid[-1] + 1)
    else:
        planeid = [1]

    key1 = 'STRING sgm_create_plane_po_created_ids[VIRTUAL]'
    key2 = f'sgm_const_plane_point_vector( "{planeid[-1]}", "Point {p1}", "Coord {coordid}.{axis}", sgm_create_plane_po_created_ids )'

    return planeid, key1, key2
