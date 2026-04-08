def createfieldnodesvector(fieldid, nodes, vector):
    nodeid = []
    vecid = []

    for kk in range(len(nodes)):
        if kk < len(nodes) - 1:
            nodeid.append(f'"Node {int(nodes[kk])}",')
            vecid.append(f' "<{vector[kk,0]:.4f},{vector[kk,1]:.4f},{vector[kk,2]:.4f}>", ')
        else:
            nodeid.append(f'"Node {int(nodes[kk])}"')
            vecid.append(f' "<{vector[kk,0]:.4f},{vector[kk,1]:.4f},{vector[kk,2]:.4f}>" ')

    nodeid1 = ''.join(nodeid)
    vecid1 = ''.join(vecid)

    key1 = f'fields_create_dfem( "{fieldid}", "Node", "Vector", {len(nodes)}, [{nodeid1}], [{vecid1}] )'
    return key1
