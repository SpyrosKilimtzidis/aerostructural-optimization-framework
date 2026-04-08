import numpy as np
def creategroup(grpid, surfs):
    if isinstance(surfs, (list, np.ndarray)):
        surfs_str = ' '.join(str(s) for s in np.ravel(surfs))
    else:
        surfs_str = str(surfs)

    key1 = f'ga_group_create( "{grpid}" )'
    key2 = f'ga_group_entity_add( "{grpid}", "Surface {surfs_str}" )'
    key3 = f'bv_group_reorganize( 1, ["{grpid}"], [0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0] )'

    return key1, key2, key3
