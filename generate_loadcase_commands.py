import pcl_loads as pcl_loads
import pcl_fields as pcl_fields

def generate_loadcase_commands(base_name, node_array, force_array, step_size=100):
    """
    Generates Patran SES commands for nodal vector fields and associated forces.

    Parameters:
        base_name (str): Prefix for field and force IDs (e.g., 'wing', 'terib', 'tiprib').
        node_array (np.ndarray): Array of shape (N, 4) where the 4th column contains node IDs.
        force_array (np.ndarray): Array of shape (N, 3) containing force vectors.
        step_size (int): Number of nodes per SES command group.

    Returns:
        data_ses (list): Nested list of SES command strings (like data1[0] format).
        lc_str_list (list): List of force ID strings used (e.g., for group referencing).
    """
    data_ses = [[]]  # Equivalent to data1{1}
    jk = 0
    total_nodes = node_array.shape[0]
    list1 = total_nodes // step_size
    lc_str_list = []

    for ikk in range(list1 + 1):
        if ikk == 0:
            stepini = 0
        elif ikk > 0 and ikk < list1:
            stepini = ikk * step_size
        elif ikk == list1:
            stepini = ikk * step_size
            step_size = total_nodes - stepini

        node_range = slice(stepini, stepini + step_size)
        node_ids = node_array[node_range, 3]
        vectors = force_array[node_range, :]

        field_id = f'lc1_nodalforces_{base_name}_{ikk + 1}'
        force_id = f'lc1_{base_name}_{ikk + 1}'

        # Create SES lines
        data_ses[0].append(pcl_fields.createfieldnodesvector(field_id, node_ids, vectors))
        data_ses[0].append(pcl_loads.createforces_id(force_id, field_id, node_ids))

    return data_ses
