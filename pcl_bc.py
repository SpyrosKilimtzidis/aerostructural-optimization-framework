def createclamp(surfaces):
    """
    Generates a Patran command string to apply a clamp boundary condition
    to the specified surface IDs.

    Parameters:
    - surfaces: list or array of surface IDs

    Returns:
    - key1: string containing the formatted clamp command
    """
    return (
        'loadsbcs_create2( "clamp", "Displacement", "Nodal", "", "Static", '
        f'["Surface {surfaces[0]}:{surfaces[-1]}"], '
        '"Geometry", "Coord 0", "1.", '
        '["<0,0,0     >", "<0,0,0     >", "<     >", "<     >"], '
        '["", "", "", ""] )'
    )
