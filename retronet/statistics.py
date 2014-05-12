import collections


def get_stats(graph, max_depth, root):
    """Counts nodes at given depth."""

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
