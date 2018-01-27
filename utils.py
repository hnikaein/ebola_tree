def upgma(names, edit_dist):
    """
    performing upgma algorithm
    :param names: names of elements
    :param edit_dist: matrix of edit distances between elements
    :return: a newick style tree
    """
    names = names[:]
    dist = {}
    size = {}
    depth = {}
    n = len(names)
    for name in names:
        size[name] = 1
        depth[name] = 0
    for i in range(n):
        dist[names[i]] = {}
        for j in range(n):
            dist[names[i]][names[j]] = edit_dist[i][j]
    for k in range(n - 1):
        # find the nearest clusters
        maxx = 10000000
        maxii = (0, 0)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                if dist[names[i]][names[j]] < maxx:
                    maxx = dist[names[i]][names[j]]
                    maxii = (i, j)
        # compute edges and nodes
        maxx /= 2
        maxi = names[maxii[0]], names[maxii[1]]
        edge_0 = maxx - depth[maxi[0]]
        edge_1 = maxx - depth[maxi[1]]
        if edge_1 < edge_0 or (edge_0 == edge_1 and maxi[0] > maxi[1]):
            maxi = maxi[1], maxi[0]
            edge_0, edge_1 = edge_1, edge_0
        # compute new nod attributes
        new_name = "(%s:%d,%s:%d)" % (maxi[0], edge_0, maxi[1], edge_1)
        dist[new_name] = {}
        size[new_name] = size[maxi[0]] + size[maxi[1]]
        depth[new_name] = maxx
        for i in names:
            dist[new_name][i] = (size[maxi[0]] * dist[maxi[0]][i] + size[maxi[1]] * dist[maxi[1]][i]) / size[
                new_name]
            dist[i][new_name] = (size[maxi[0]] * dist[maxi[0]][i] + size[maxi[1]] * dist[maxi[1]][i]) / size[
                new_name]
        names.append(new_name)
        dist[new_name][new_name] = 0
        # delete old nodes
        del dist[maxi[0]]
        del dist[maxi[1]]
        del names[maxii[1]]
        del names[maxii[0]]
    return names[0] + ";"


def neighbor_join(names, edit_dist):
    """
    performing neighbor_join algorithm
    :param names: names of elements
    :param edit_dist: matrix of edit distances between elements
    :return: a newick style tree
    """
    names = names[:]
    dist = {}
    for i in range(len(names)):
        dist[names[i]] = {}
        for j in range(len(names)):
            dist[names[i]][names[j]] = edit_dist[i][j]
    new_name = names[0]
    for k in range(len(names) - 2):
        # find the nearest nodes
        maxx = 0
        maxii = (0, 0)
        d = {}
        for name in names:
            d[name] = sum([dist[name][nname] for nname in names])
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                if d[names[i]] + d[names[j]] - (len(names) - 2) * dist[names[i]][names[j]] > maxx:
                    maxx = d[names[i]] + d[names[j]] - (len(names) - 2) * dist[names[i]][names[j]]
                    maxii = (i, j)
        # compute edges and nodes
        maxi = names[maxii[0]], names[maxii[1]]
        edge_0 = 0.5 * (dist[maxi[0]][maxi[1]] + (d[maxi[0]] - d[maxi[1]]) / (len(names) - 2))
        edge_1 = dist[maxi[0]][maxi[1]] - edge_0
        if edge_1 < edge_0 or (edge_0 == edge_1 and maxi[0] > maxi[1]):
            maxi = maxi[1], maxi[0]
            edge_0, edge_1 = edge_1, edge_0
        # compute new nod attributes
        new_name = "(%s:%d,%s:%d)" % (maxi[0], edge_0, maxi[1], edge_1)
        dist[new_name] = {}
        for i in names:
            dist[new_name][i] = 0.5 * (dist[maxi[0]][i] + dist[maxi[1]][i] - dist[maxi[0]][maxi[1]])
            dist[i][new_name] = dist[new_name][i]
        names.append(new_name)
        dist[new_name][new_name] = 0
        # delete old nodes
        del dist[maxi[0]]
        del dist[maxi[1]]
        del names[maxii[1]]
        del names[maxii[0]]
    # compute last node containing three last nodes
    last_name = names[1] if new_name == names[0] else names[0]
    total = "%s,%s:%d);" % (new_name[:-1], last_name, dist[names[0]][names[1]])
    return total
