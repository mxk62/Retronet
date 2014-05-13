import collections


def count_nodes(graph, root, max_depth):
    """Counts nodes at given depth.

    Parameters
    ----------
    graph : Networkx DiGraph
        A graph on which counting is done.
    root : a graph node
        A seed node for which counting is done.
    max_depth : integer
        Maximal depth for which counting is performed.

    Returns
    -------
    count : dictionary
        A dictionary which keys represents depths and values number of nodes at
        that depth.
    """

    # Usually we measure depth in synthetic steps (no. of reactions
    # separating two chemicals, but here we want to mean it what it usually
    # does -- a step between adjacent nodes.
    max_depth *= 2

    # Initialize counters (values) at given depths (keys).
    occupation = {depth: 0 for depth in range(max_depth)}

    discovered = {root}
    queue = collections.deque([(root, 0)])
    while queue:
        child, depth = queue.popleft()

        # Do not procceed beyond a given depth.
        if depth == max_depth:
            continue

        # Find node's predecessors.
        parents = set(graph.predecessors(child))
        undiscovered = parents - discovered

        # Update counters.
        occupation[depth] += len(undiscovered)

        # Enqueue newly discovered nodes.
        queue.extend([(u, depth + 1) for u in undiscovered])
        discovered |= undiscovered
    return occupation
