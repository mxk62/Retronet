import itertools
import networkx
import time


def populate(connections, smiles, transforms, depth=2):
    """Builds a retrosynthetic chemical network around a given chemical.

    Parameters
    ----------
    connections : list of Connections
        A sequence of connections to available workers.
    smiles : string
        A SMILES representing target molecule.
    transforms : sequence of Transform
        A sequence of retrosynthetic transforms which should be used to
        generated reactions.
    depth : integer (optional)
        Maximal number of synthetics steps, defaults to 2.

    Returns
    -------
    g : Networkx DiGraph
        A directed, biparite graph representing the chemical network.
    """
    batch_size = 100
    batches = partition(transforms.keys(), batch_size)

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
                # available transforms.
                if current_depth < depth:

                    # Distribute transforms to perform among the workers.
                    for no, batch in enumerate(batches):
                        msg = {'smiles': chem_smi, 'trans_ids': batch}
                        connections[no % len(connections)].send(msg)

                        print 'Ventilator: sent \'%s\' to %d.' % (msg, no % len(connections))
                        time.sleep(1)

                    # Collect results from different workers.
                    reactant_sets = set()
                    for i in range(len(batches)):
                        msg = connections[i % len(connections)].recv()

                        print 'Ventilator: received \'%s\' from %d.' % (msg, no % len(connections))
                        time.sleep(1)

                        if msg['results'] != 'NONE':
                            reactant_sets.update(
                                frozenset(smis) for smis in msg['results'])

                    # Add reaction associated with each reactant set and
                    # enqueue chemicals from those sets.
                    for reactants in reactant_sets:
                        rxn_no = rxn_counter.next()
                        graph.add_node(rxn_no, bipartite=1)
                        graph.add_edge(rxn_no, chem_smi)
                        next_nodes.update((smi, rxn_no) for smi in reactants)

            # Update outgoing edges if there is a successor (a reaction
            # which uses it as a reactant).
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
