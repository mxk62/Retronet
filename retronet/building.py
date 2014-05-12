import collections
import itertools
import networkx as nx
import retronet as rn


def populate(smiles, transforms, depth=2):
    """Builds a new retrosynthetic network around a given chemical.

    Parameters
    ----------
    smiles : string
        A SMILES representing initial target molecule.
    transforms : Transform sequence
        A sequence of retrosynthetic transforms which should be used to
        generated reactions.
    depth : integer (optional)
        Maximal number of synthetics steps, defaults to 2.

    Returns
    -------
    g : Networkx DiGraph
        A directed, biparite graph representing the chemical network.
    """
    graph = nx.DiGraph()

    rxn_count = itertools.count()

    processed = set()
    queue = collections.deque([(smiles, None, depth)])
    while queue:
        chem_smi, rxn_id, depth = queue.popleft()

        # Add a node to the graph if it is not already there.
        if not graph.has_node(chem_smi):
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
        chem = rn.Chemical(chem_smi)
        reactant_sets = set()
        for t in transforms:
            reactant_sets.update(chem.make_retrostep(t))

        # Add reaction associated with each reactant set and enqueue chemicals
        # from product sets.
        for reactants in reactant_sets:
            rxn_no = rxn_count.next()
            graph.add_node(rxn_no, bipartite=1)
            graph.add_edge(rxn_no, chem_smi)
            queue.extend([(smi, rxn_no, depth - 1) for smi in reactants])

        processed.add(chem_smi)

    return graph


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
    graph = nx.DiGraph()

    processed = set()
    queue = collections.deque([(seed_smiles, None, depth)])
    while queue:
        chem_smi, rxn_id, depth = queue.popleft()

        # Add a node to the graph if it is not already there.
        if not graph.has_node(chem_smi):
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
