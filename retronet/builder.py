import itertools
import json
import networkx
import zmq


def populate(smiles, transforms, depth=2):
    """Builds a retrosynthetic chemical network around a given chemical.

    Parameters
    ----------
    smiles : string
        A SMILES representing the target chemical.
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
    context = zmq.Context()

    # Set up a channel to communicate with workers.
    socket = context.socket(zmq.PUSH)
    socket.bind('tcp://*:5555')

    # Set up a channel to gather results
    receiver = context.socket(zmq.PULL)
    receiver.bind('tcp://*:5556')

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
                    for batch in batches:
                        message = {'smiles': chem_smi, 'trans_ids': batch}
                        socket.send_json(json.dumps(message))

                    # Collect results from different workers.
                    reactant_sets = set()
                    for i in range(len(batches)):
                        message = json.loads(receiver.recv_json())
                        if message['results'] != 'NONE':
                            reactant_sets.update(
                                frozenset(smis) for smis in message['results'])

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


def partition(list, n):
    return [list[i:i+n] for i in range(0, len(list), n)]
