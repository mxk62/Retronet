import argparse
import collections
import itertools
import networkx as nx
import sys
import time
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
from retronet.transform import Transform
from retronet.chemical import Chemical


def populate(graph, max_depth, seed_smiles, transforms):
    """Builds a new retrosynthetic network around a given chemical."""

    rxn_count = itertools.count()

    processed = set()
    queue = collections.deque([(seed_smiles, None, 0)])
    while queue:
        chem_smi, rxn_id, depth = queue.popleft()

        # Add a node to the graph if it is not already there.
        if not graph.has_node(chem_smi):
            graph.add_node(chem_smi, bipartite=0)

        # Update its outgoing edges if there is a successor (a reaction which
        # uses it as a reactant).
        if rxn_id is not None:
            graph.add_edge(chem_smi, rxn_id)

        # Do NOT processed further if node is at maximal depth.
        if chem_smi in processed or depth == max_depth:
            continue

        # Find out valid retrosynthetic reactions using the available
        # transforms.
        chem = Chemical(chem_smi)
        reactant_sets = set()
        for t in transforms:
            reactant_sets.update(chem.make_retrostep(t))

        # Add reaction associated with each reactant set and enqueue chemicals
        # from product sets.
        for reactants in reactant_sets:
            rxn_no = rxn_count.next()
            graph.add_node(rxn_no, bipartite=1)
            graph.add_edge(rxn_no, chem_smi)

            queue.extend([(smi, rxn_no, depth + 1) for smi in reactants])

        processed.add(chem_smi)


def recreate(graph, max_depth, seed_smiles, db):
    """Recreates a fragment of existing chemical network."""

    processed = set()
    queue = collections.deque([(seed_smiles, None, 0)])
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
        if chem_smi in processed or depth == max_depth:
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

            queue.extend([(smi, rxn_id, depth + 1) for smi in reactants])

        processed.add(chem_smi)


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


# Parse command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument('--ip', type=str, default='127.0.0.1',
                    help='Database server ip.')
parser.add_argument('-d', '--depth', type=int, default=3,
                    help='Maximal network depth.')
parser.add_argument('file', type=str,
                    help='File with chemicals SMILES.')
args = parser.parse_args()

# Connect to the mongo database.
try:
    client = MongoClient(args.ip)
except ConnectionFailure:
    sys.exit(1)
db = client.data

# Acquire all valid, single-retron transforms.
transforms = []
for idx, rec in enumerate(db.retro.find(), start=1):
    try:
        t = Transform(rec['reaction_smarts'].encode('ascii'))
    except (KeyError, ValueError):
        continue
    if len(t.retrons) == 1:
        transforms.append(t)

# Read the SMILES from an external file.
smiles = []
with open(args.file) as f:
    for line in f.readlines():
        smiles.append(line.strip())

# Build the network, either using existing reaction database or by applying
# retrosynthetic transforms.
network = nx.DiGraph()
for idx, smi in enumerate(smiles):
    # Time each build individually.
    start = time.clock()
    populate(network, args.depth, smi, transforms)
    #recreate(network, args.depth, smi, db)
    eta = time.clock() - start

    # Gather statistics and write them down to a file.
    stats = get_stats(network, args.depth, smi)
    line = '\t'.join([str(stats[key]) for key in sorted(stats)])
    with open('occupation.dat', 'a') as f:
        f.write('{0}\t{1}\t{2}\n'.format(idx, eta, line))

    # Write graph for debugging purposes
    #nx.write_gexf(network, 'testgraph.gexf')

    # Discard current network entirely.
    network.clear()
