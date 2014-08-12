import collections
import itertools
import networkx
from retronet import Chemical


def populate_serial(smiles, transforms, depth=2):
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

                if current_depth < depth:
                    # Find out valid retrosynthetic reactions using the
                    # available transforms.
                    chem = Chemical(chem_smi)
                    reactant_sets = apply(chem, transforms)

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
    graph = networkx.DiGraph()

    rxn_count = itertools.count()

    processed = set()
    queue = collections.deque([(smiles, None, depth)])
    while queue:
        chem_smi, rxn_id, depth = queue.popleft()

        # Add a node to the graph if it is not already there.
        if chem_smi not in graph:
            graph.add_node(chem_smi, bipartite=0)

        # Update its outgoing edges if there is a successor (a reaction which
        # uses it as a reactant).
        if rxn_id is not None:
            graph.add_edge(chem_smi, rxn_id)

        # Do NOT proceed further if node was already processed or is at
        # maximal depth.
        if chem_smi in processed or depth == 0:
            continue

        # Find out valid retrosynthetic reactions using the available
        # transforms.
        chem = Chemical(chem_smi)
        reactant_sets = apply(chem, transforms)

        # Add reaction associated with each reactant set and enqueue chemicals
        # from those sets.
        for reactants in reactant_sets:
            rxn_no = rxn_count.next()
            graph.add_node(rxn_no, bipartite=1)
            graph.add_edge(rxn_no, chem_smi)
            queue.extend([(smi, rxn_no, depth - 1) for smi in reactants])

        processed.add(chem_smi)

    return graph


def apply(chem, transforms):
    """Applies available transform to a given chemical.

    Parameters
    ----------
    chem : Chemical
        A representation of a given chemical compound.
    transforms : sequence of Transforms
        List of availalbe transforms.

    Returns
    -------
    reactants : set
        Possible reactant sets obtained after applying all compatible
        transforms.
    """
    reactant_sets = set()
    for t in transforms:
        try:
            current_sets = chem.make_retrostep(t)
        except RuntimeError:
            continue
        reactant_sets |= current_sets
    return reactant_sets


def recreate(seed_smiles, db, depth=2):
    """Recreates a fragment of existing network around a given chemical.

    Parameters
    ----------
    smiles : string
        A SMILES representing initial target molecule.
    db : mongo database
        A mongo database storing the existing chemical network.
    depth : integer (optional)
        Maximal number of synthetics steps, defaults to 2.

    Returns
    -------
    g : Networkx DiGraph
        A directed, biparite graph representing the chemical network.
    """
    graph = networkx.DiGraph()

    processed = set()
    queue = collections.deque([(seed_smiles, None, depth)])
    while queue:
        chem_smi, rxn_id, depth = queue.popleft()

        # Add a node to the graph if it is not already there.
        if chem_smi not in graph:
            graph.add_node(chem_smi, bipartite=0)

        # Update its outgoing edges if there is a successor (a reaction which
        # uses it as a reactant).
        if rxn_id is not None:
            graph.add_edge(chem_smi, rxn_id)

        # Do NOT processed further if node is as already processed or is at
        # maximal depth.
        if chem_smi in processed or depth == 0:
            continue

        # Find chemical predecessors (reactions which produce it).
        rxn_ids = db.chemicals.find_one({'smiles': chem_smi})['reactions'][
            'produced']

        rxn_query = {'_id': {'$in': rxn_ids}}
        for rxn_rec in db.reactions.find(rxn_query):
            rxn_id = rxn_rec['_id']
            graph.add_node(rxn_id, bipartite=1)
            graph.add_edge(rxn_id, chem_smi)

            chem_query = {'_id': {'$in': rxn_rec['reactants']}}
            reactants = set(chem_rec['smiles']
                            for chem_rec in db.chemicals.find(chem_query))

            queue.extend([(smi, rxn_id, depth - 1) for smi in reactants])

        processed.add(chem_smi)

    return graph
