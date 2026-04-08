def createcoord(coordid, p1, p2, p3):
    if coordid:
        coordid.append(coordid[-1] + 1)
    else:
        coordid = [1]

    key1 = 'STRING asm_create_cord_axi_created_ids[VIRTUAL]'
    key2 = (
        f'asm_const_coord_axis( "{coordid[-1]}", "XY", "Coord 0", 1, '
        f'"Point {p1}", "Point {p2}", "Point {p3}", asm_create_cord_axi_created_ids )'
    )
    return key1, key2


def createcoordeuler(coordid, p1):
    if coordid:
        coordid.append(coordid[-1] + 1)
    else:
        coordid = [1]

    key1 = 'STRING asm_create_cord_eul_created_ids[VIRTUAL]'
    key2 = f'asm_const_coord_euler( "{coordid[-1]}", 3, 1, 3, 0., 0., 0., "Coord 0", 1, "Point {p1}", asm_create_cord_eul_created_ids )'

    return coordid, key1, key2
