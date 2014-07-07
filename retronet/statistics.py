import collections


def count_nodes(graph, root, max_depth):
    """Counts nodes at given depth range.

    Parameters
    ----------
    graph : Networkx DiGraph
        A graph on which counting is done.
    root : a graph node
        A seed node for which counting is done.
    max_depth : integer
        Maximal depth for which nodes are counted.

    Returns
    -------
    counts : dictionary
        A dictionary which keys represents depths (node distances) from the
        root and their values the total number of nodes at a given distance.
    """

    # Usually we measure depth in synthetic steps (no. of reactions separating
    # two chemicals), but here we want to mean it what it usually does --
    # a distance between adjacent nodes.
    max_depth *= 2

    # Initialize counters (values) at given depths (keys).
    counts = {depth: 0 for depth in range(max_depth + 1)}

    discovered = {root}
    queue = collections.deque([(root, 0)])
    while queue:
        child, depth = queue.popleft()
        counts[depth] += 1

        # Do not procceed beyond a given depth.
        if depth == max_depth:
            continue

        # Find node's predecessors.
        parents = set(graph.predecessors(child))
        undiscovered = parents - discovered

        # Enqueue newly discovered nodes.
        queue.extend([(u, depth + 1) for u in undiscovered])
        discovered |= undiscovered
    return counts
