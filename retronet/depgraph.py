import networkx
import retronet as rn


def create_depgraph(patterns):
    """Arranges molecular patterns in hierarchical manner.

    Parameters
    ----------
    patterns : Pattern sequence
        Sequence of patterns for which hierarchy is to be built.

    Returns
    -------
    graph : Networkx DiGraph
        A directed, acyclic graph representing dependencies between patterns.
    """
    g = networkx.DiGraph()
    while patterns:
        patt = patterns.pop()
        if patt not in g:
            g.add_node(patt, id=patt.smiles)

        duplicates = set()
        for other in patterns:
            if patt.does_contain(other) or other.does_contain(patt):
                d = patt.find_distance(other)
                if abs(d) < rn.EPSILON:
                    duplicates.add(other)
                else:
                    if other not in g:
                        g.add_node(other)
                    edge = (other, patt) if patt.does_contain(other)\
                        else (patt, other)
                    g.add_edge(*edge, weight=d)
        if duplicates:
            g.node[patt]['duplicates'] = [d.smiles for d in duplicates]
        for patt in duplicates:
            patterns.remove(patt)
            if patt in g:
                g.remove_node(patt)
    return prune(g)


def prune(graph):
    """Prune graph from unwanted edges.

    Parameters
    ----------
    graph : Networkx DiGraph
        A graph representing 'dependencies' between patterns.

    Returns
    -------
    graph : Networkx DiGraph
        The graph with edges pruned.
    """
    dependencies = {v: set(u for u, _ in graph.in_edges_iter([v]))
                    for v in graph.nodes()}
    for level, nodes in enumerate(toposort(dependencies)):
        for v in nodes:
            graph.node[v]['lvl'] = level
    for start in graph.nodes():
        for end in graph.successors(start):
            if graph.node[end]['lvl'] - graph.node[start]['lvl'] > 1:
                graph.remove_edge(start, end)
    for v in graph.nodes():
        del graph.node[v]['lvl']
    return graph


def toposort(data):
    """Performs topological sorting on data.

    Parameters
    ----------
    data : dictionary
        A dictionary whose keys are items and values are sets of dependent
        items.

    Returns
    -------
    ordering : set
        A list of sets in topological order. The first set consists of items
        with no predecessors, each subsequent set consists of items that
        are successors of items in the preceding sets.
    """

    # Ignore self dependencies.
    for k, v in data.items():
        v.discard(k)

    # Find items which do not depend on anything and add empty dependencies
    # where needed.
    extra_items = reduce(set.union, data.values()) - set(data.keys())
    data.update({item: set() for item in extra_items})

    while True:
        # Find nodes having no dependencies.
        ordered = set(item for item, dep in data.items() if not dep)
        if not ordered:
            break
        yield ordered

        # Remove nodes having no dependencies and update the list of
        # adjacent nodes of the remaining ones accordingly.
        data = {item: (dep - ordered) for item, dep in data.items()
                if item not in ordered}
    if data:
        raise ValueError('Cycle(s) exist among input data.')
