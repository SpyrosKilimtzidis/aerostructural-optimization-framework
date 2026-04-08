import numpy as np
def createpcompprops(laminate, nmat1, lamid, hply, symid):
    # Determine symmetry flag
    if symid == 1:
        symm = 3
    elif symid==0:
        symm = 2
    else:
        symm = 1
# Repeat material ID string
    string1 = ', '.join([f'"{nmat1}"' for _ in range(len(laminate))])

    # Repeat hply values as strings
    string2 = ', '.join([f'"  {hply}"' for _ in range(len(laminate))])

    # Laminate angles
    string3 = ', '.join([f'"  {str(angle)}"' for angle in laminate])

    # Empty fields
    string4 = ', '.join(['""' for _ in range(len(laminate))])

    key1 = (
        f'mat_create_lam3( "{lamid}", "", {symm}, [{string1}], '
        f'[{string2}], [{string3}], [{string4}], {len(laminate)}, "", "Create" )'
    )

    return key1

def createpshellprops(propid, laminateid, coordid, surfs):
    # Convert list or array of surface IDs to a space-separated string
    if isinstance(surfs, (list, tuple, np.ndarray)):
        surf_str = ' '.join(str(s) for s in np.ravel(surfs))
    else:
        surf_str = str(surfs)

    key1 = (
        f'elementprops_create( "{propid}", 51, 25, 35, 1, 3, 20, '
        f'[13, 20, 4037, 4111, 4265, 1005, 5, 8111, 4213], '
        f'[5, 9, 1, 1, 1, 1, 1, 4, 4], '
        f'["m:{laminateid}", "Coord {coordid}", "", "", "", "", "", "", ""], '
        f'"Surface {surf_str}" )'
    )
    return key1

def createcbeamprops(curves, w, h, mat, propname, sectid, orient, tskin, identifier):
    if identifier == 1:
        offset = tskin / 2 + h / 2
    else:
        offset = -(tskin / 2 + h / 2)

    # Convert curves to a space-separated string
    if isinstance(curves, (list, tuple, np.ndarray)):
        curves_str = ' '.join(str(c) for c in np.ravel(curves))
    else:
        curves_str = str(curves)

    # Section creation command
    key1 = f'beam_section_create( "{sectid}", "BAR", ["{w:.6f}", "{h:.6f}"] )'

    # Property creation command
    key2 = (
        f'elementprops_create( "{propname}", 11, 2, 42, 1, 1, 20, '
        f'[39, 13, 6, 4042, 4043, 2047, 2048, 1, 10, 11, 4026, 1026, 4044, 4045, '
        f'4037, 4047, 4048, 4050, 4051, 4053, 4054, 4056, 4057, 8112, 4061, 4303, '
        f'8111, 4403, 4404, 4410, 4411, 8200, 8201, 8202], '
        f'[11, 5, 2, 2, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, '
        f'4, 4, 1, 1, 1, 6, 4, 4, 4], '
        f'["{sectid}", "m:{mat}", "<{orient[0]} {orient[1]} {orient[2]}>", '
        f'"<0 0 {offset:.6f}>", "<0 0 {offset:.6f}>", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", '
        f'"Analysis", "Analysis", "Analysis"], '
        f'"Curve {curves_str}" )'
    )

    return key1, key2


