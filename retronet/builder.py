import itertools
import networkx


def populate(tasks, results, smiles, transforms, depth=2):
    """Builds a retrosynthetic network around a given chemical.

    Function use breath-first search approach to generate retrosynthetic
    network around a specified target compound.

    Parameters
    ----------
    tasks : Queue
        A 'ventilator', distributes tasks among workers.
    results : Queue
        A 'sink', collects results from workers.
    smiles : string
        A target compound encoded in Daylight's SMILES notation.
    transforms : sequence of Transforms
        A sequence of all available retrosynthetic transforms which can be
        used to generate reactions.
    depth : integer (optional)
        Maximal number of synthetics steps, defaults to 2.

    Returns
    -------
    g : Networkx's DiGraph
        A directed, biparite graph representing the chemical network.
    """

    # Divide available transforms into batches.
    batch_size = 3000
    batches = partition(transforms.keys(), batch_size)

    # Initialize graph representing the network.
    graph = networkx.DiGraph()

    rxn_counter = itertools.count()
    lvl_counter = itertools.count()

    current_depth = lvl_counter.next()
    current_nodes = {(smiles, None)}
    while current_nodes:
        next_nodes = set()

        for chem_smi, rxn_id in current_nodes:

            # Add a node to the graph if it is not already there.
            if chem_smi not in graph:
                graph.add_node(chem_smi, bipartite=0)

                # Find out valid retrosynthetic reactions using the
                # available transforms, providing we didn't reach the
                # terminal depth.
                if current_depth < depth:

                    # Distribute transforms to perform among the workers.
                    for batch in batches:
                        msg = {'smiles': chem_smi, 'trans_ids': batch}
                        tasks.put(msg)

                    # Collect results from different workers.
                    reactant_sets = set()
                    for i in range(len(batches)):
                        msg = results.get()
                        if msg['results'] is not None:
                            reactant_sets.update(
                                frozenset(smis) for smis in msg['results'])

                    # Add reaction associated with each reactant set and
                    # enqueue chemicals from those sets.
                    for reactants in reactant_sets:
                        rxn_no = rxn_counter.next()
                        graph.add_node(rxn_no, bipartite=1)
                        graph.add_edge(rxn_no, chem_smi)
                        next_nodes.update((smi, rxn_no) for smi in reactants)

            # Finally, update outgoing edges if there is a successor
            # (a reaction which uses the chemical as a reactant).
            if rxn_id is not None:
                graph.add_edge(chem_smi, rxn_id)

        current_nodes = set(next_nodes)
        current_depth = lvl_counter.next()

    return graph


def partition(l, n):
    """Partitions a list into nonoverlapping sublists of a given length.

    Parameters
    ----------
    l : list
        A list to be partitioned.
    n : integer
        Requested chunk size.

    Returns
    -------
    p : list
        A sequence of list representing result of the partitioning.
    """
    return [l[i:i+n] for i in range(0, len(l), n)]
