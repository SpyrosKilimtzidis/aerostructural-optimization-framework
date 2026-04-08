def settol(tol):
    key1 = f'pref_geo_set_v1( 0, {tol}, 3 )'
    key2 = f'pref_global_set_v3( TRUE, 3, "{tol}", "SAVE" )'
    key3 = 'pref_env_set_logical( "revert_enabled", FALSE )'
    return [key1, key2, key3]
